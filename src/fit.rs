//! 基准变换的拟合：由（线性输入, 参考输出）配对求解。
//!
//! # 配对的两端分别从哪来
//!
//! ```text
//! 输入的线性端：解码 NEF → 相机矩阵 → ProPhoto 线性（本管线的工作空间）
//! 参考的线性端：读参考 TIF → 按内嵌 ICC 线性化 → ProPhoto 线性
//! ```
//!
//! 两端都在 ProPhoto 线性上，所以它们之间的映射就是**基准渲染本身**——不含任何
//! 色域转换的干扰。
//!
//! # 为什么要拆成「一维色调曲线 + 三维色彩查找表」
//!
//! 基准渲染同时包含非线性色调与色彩重映射：纯矩阵表达不了色调曲线，纯 3D LUT 覆盖
//! 不了无限动态范围。拆开后各自可拟合、可单独诊断。
//!
//! # 拟合的第一步是"看"，不是"解"
//!
//! 先把配对关系画出来（分箱中位数），看它是不是一条平滑单调的曲线。**如果不是**，
//! 说明参考样本里混进了非基准的因素（某张被调整过、ADL 没关掉、配对错了），这时
//! 求解只会把问题一起拟合进去。所以 [`diagnose`] 先于 [`fit_curve`]。

use crate::error::{Error, Result};
use crate::icc::ColorSpacePlan;

/// 一张 16 位 RGB 图像。
#[derive(Debug, Clone)]
pub struct Image16 {
    pub width: usize,
    pub height: usize,
    /// 行主序，每像素 3 个分量。
    pub pixels: Vec<u16>,
}

impl Image16 {
    #[inline]
    pub fn at(&self, x: usize, y: usize) -> Option<[u16; 3]> {
        if x >= self.width || y >= self.height {
            return None;
        }
        let i = (y * self.width + x) * 3;
        Some([self.pixels[i], self.pixels[i + 1], self.pixels[i + 2]])
    }
}

/// 读一张**未压缩、单平面、16 位**的 RGB TIFF。
///
/// 只支持这一种——参考导出的格式由尼康工坊决定且固定。遇到别的组合会明确报错，
/// 而不是猜着读：猜错的代价是整批标定数据静默偏移。
pub fn read_rgb16(data: &[u8]) -> Result<Image16> {
    use crate::tiff::Tiff;
    let t = Tiff::locate(data)?;
    let ifd = t.ifd0_offset()?;

    // 按条目自身的类型码取值——TIFF 里同一个标签的宽度由写文件的实现决定，
    // 一律按 LONG 读会在 SHORT 上炸（`Truncated { LONG, need 4, have 2 }`）。
    fn as_u32(t: &crate::tiff::Tiff<'_>, e: &crate::tiff::Entry) -> Result<u32> {
        match e.typ {
            3 => Ok(t.u16_value(e)? as u32),
            4 => t.u32_value(e),
            other => Err(Error::BadPayload {
                what: "参考 TIFF",
                detail: format!("标签 {:#06x} 的类型码为 {other}，只支持 SHORT 与 LONG", e.tag),
            }),
        }
    }

    let get = |tag: u16, what: &'static str| -> Result<u32> {
        let e = t.find(ifd, tag)?.ok_or_else(|| Error::BadPayload {
            what: "参考 TIFF",
            detail: format!("缺少{what}标签（{tag:#06x}）"),
        })?;
        as_u32(&t, &e)
    };

    let width = get(0x0100, "宽度")? as usize;
    let height = get(0x0101, "高度")? as usize;

    // 位深：可能是单值，也可能是每通道一个
    let bits_entry = t.find(ifd, 0x0102)?.ok_or_else(|| Error::BadPayload {
        what: "参考 TIFF",
        detail: "缺少位深标签".into(),
    })?;
    let bits = as_u32(&t, &bits_entry)?;
    if bits != 16 {
        return Err(Error::BadPayload {
            what: "参考 TIFF",
            detail: format!("位深为 {bits}，本读取器只支持 16 位"),
        });
    }
    let spp = get(0x0115, "通道数")?;
    if spp != 3 {
        return Err(Error::BadPayload {
            what: "参考 TIFF",
            detail: format!("通道数为 {spp}，本读取器只支持 3"),
        });
    }
    let compression = match t.find(ifd, 0x0103)? {
        Some(e) => as_u32(&t, &e)?,
        None => 1,
    };
    if compression != 1 {
        return Err(Error::BadPayload {
            what: "参考 TIFF",
            detail: format!("压缩方式为 {compression}，本读取器只支持未压缩"),
        });
    }
    let planar = match t.find(ifd, 0x011C)? {
        Some(e) => as_u32(&t, &e)?,
        None => 1,
    };
    if planar != 1 {
        return Err(Error::BadPayload {
            what: "参考 TIFF",
            detail: "只支持单平面（chunky）排列".into(),
        });
    }

    // 条带
    let offs_e = t.find(ifd, 0x0111)?.ok_or_else(|| Error::BadPayload {
        what: "参考 TIFF",
        detail: "缺少条带偏移标签".into(),
    })?;
    let counts_e = t.find(ifd, 0x0117)?.ok_or_else(|| Error::BadPayload {
        what: "参考 TIFF",
        detail: "缺少条带字节数标签".into(),
    })?;
    let offsets: Vec<u32> = read_u32_list(&t, &offs_e)?;
    let counts: Vec<u32> = read_u32_list(&t, &counts_e)?;
    if offsets.len() != counts.len() {
        return Err(Error::BadPayload {
            what: "参考 TIFF",
            detail: format!("条带偏移 {} 个而字节数 {} 个，不匹配", offsets.len(), counts.len()),
        });
    }

    let need = width * height * 3 * 2;
    let mut pixels = Vec::with_capacity(need / 2);
    for (o, c) in offsets.iter().zip(counts.iter()) {
        let s = *o as usize;
        let e = s + *c as usize;
        if e > data.len() {
            return Err(Error::BadPayload {
                what: "参考 TIFF",
                detail: format!("条带 {s}..{e} 越出文件（长度 {}）", data.len()),
            });
        }
        for ch in data[s..e].chunks_exact(2) {
            pixels.push(u16::from_le_bytes([ch[0], ch[1]]));
        }
    }

    if pixels.len() < width * height * 3 {
        return Err(Error::BadPayload {
            what: "参考 TIFF",
            detail: format!(
                "像素数据不足：需 {} 个分量，实得 {}",
                width * height * 3,
                pixels.len()
            ),
        });
    }
    pixels.truncate(width * height * 3);

    Ok(Image16 { width, height, pixels })
}

fn read_u32_list(t: &crate::tiff::Tiff<'_>, e: &crate::tiff::Entry) -> Result<Vec<u32>> {
    // 条带标签的类型既可能是 LONG 也可能是 SHORT——由写文件的实现决定，
    // 不能假定。早前按 LONG 读，遇到 SHORT 就报 `Truncated { LONG, need 4, have 2 }`。
    let n = e.count as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    let unit = match e.typ {
        3 => 2usize, // SHORT
        4 => 4usize, // LONG
        other => {
            return Err(Error::BadPayload {
                what: "参考 TIFF",
                detail: format!("条带标签的类型码为 {other}，只支持 SHORT 与 LONG"),
            })
        }
    };
    let raw = t.raw();
    if n * unit <= 4 {
        // 值直接放在条目的值字段里
        let field = e.offset + 8;
        if field + 4 > raw.len() {
            return Err(Error::BadPayload { what: "参考 TIFF", detail: "条带值字段越界".into() });
        }
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let p = field + i * unit;
            out.push(if unit == 2 {
                u16::from_le_bytes([raw[p], raw[p + 1]]) as u32
            } else {
                u32::from_le_bytes([raw[p], raw[p + 1], raw[p + 2], raw[p + 3]])
            });
        }
        return Ok(out);
    }
    let off = t.deref_offset(e)?;
    if off + n * unit > raw.len() {
        return Err(Error::BadPayload { what: "参考 TIFF", detail: "条带数组越界".into() });
    }
    Ok((0..n)
        .map(|i| {
            let p = off + i * unit;
            if unit == 2 {
                u16::from_le_bytes([raw[p], raw[p + 1]]) as u32
            } else {
                u32::from_le_bytes([raw[p], raw[p + 1], raw[p + 2], raw[p + 3]])
            }
        })
        .collect())
}

/// 一对样本：我们的线性值 与 参考的线性值，都在 ProPhoto 线性。
#[derive(Debug, Clone, Copy)]
pub struct Sample {
    pub ours: [f32; 3],
    pub theirs: [f32; 3],
}

/// 把参考 TIF 转成工作空间的线性值。
pub fn reference_to_working(img: &Image16, plan: &ColorSpacePlan) -> Vec<[f32; 3]> {
    let inv = 1.0 / 65535.0;
    img.pixels
        .chunks_exact(3)
        .map(|p| {
            plan.to_working_space([
                p[0] as f32 * inv,
                p[1] as f32 * inv,
                p[2] as f32 * inv,
            ])
        })
        .collect()
}

/// 分箱诊断的结论。
#[derive(Debug, Clone)]
pub struct CurveDiagnosis {
    /// 每个箱的中位数：`(输入中位, 输出中位, 样本数)`。
    pub bins: Vec<(f32, f32, usize)>,
    /// 相邻箱输出中位数的**逆序**次数——平滑单调的曲线应为 0。
    pub inversions: usize,
    /// 曲线在相邻箱之间的斜率范围（仅统计样本充足的箱）。
    pub slope_range: Option<(f32, f32)>,
}

impl CurveDiagnosis {
    /// 大致是否像一条干净的单调整曲线。
    pub fn looks_clean(&self) -> bool {
        self.inversions == 0 && self.bins.iter().filter(|b| b.2 > 100).count() >= 4
    }
}

/// 按亮度分箱，给出「输入 → 输出」的中位数映射。
///
/// 这是拟合前的"看一眼"：若曲线有逆序或斜率剧烈跳变，说明样本里混进了非基准因素，
/// 此时求解只会把问题拟合进去。
pub fn diagnose(samples: &[Sample], bins: usize) -> CurveDiagnosis {
    let lum = |v: [f32; 3]| 0.2126 * v[0] + 0.7152 * v[1] + 0.0722 * v[2];
    let mut buckets: Vec<Vec<f32>> = vec![Vec::new(); bins];
    for s in samples {
        let x = lum(s.ours).clamp(0.0, 1.0);
        let b = ((x * bins as f32) as usize).min(bins - 1);
        buckets[b].push(lum(s.theirs));
    }

    let mut out = Vec::new();
    for (i, v) in buckets.iter_mut().enumerate() {
        if v.is_empty() {
            continue;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = v[v.len() / 2];
        let lo = i as f32 / bins as f32;
        out.push((lo, med, v.len()));
    }

    let mut inversions = 0;
    for w in out.windows(2) {
        if w[1].1 < w[0].1 - 1e-6 {
            inversions += 1;
        }
    }

    let slopes: Vec<f32> = out
        .windows(2)
        .filter(|w| w[0].2 > 100 && w[1].2 > 100)
        .map(|w| (w[1].1 - w[0].1) / (w[1].0 - w[0].0).max(1e-6))
        .filter(|s| s.is_finite())
        .collect();
    let slope_range = if slopes.is_empty() {
        None
    } else {
        Some((
            slopes.iter().cloned().fold(f32::INFINITY, f32::min),
            slopes.iter().cloned().fold(f32::NEG_INFINITY, f32::max),
        ))
    };

    CurveDiagnosis { bins: out, inversions, slope_range }
}
