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

// ---------------------------------------------------------------------------
// 感知均匀的编码域
// ---------------------------------------------------------------------------

/// 线性 → sRGB 编码（0..1 → 0..1）。
///
/// # 为什么曲线要在这个域上拟合
///
/// 第一版按**线性**域等距分箱：128 箱覆盖 0..1，箱宽 0.0078。实测 ΔE00 中位数
/// 2.097、P95 5.944，**未达标**，且最大误差全在深阴影（输入 0.02–0.04、参考
/// 0.09–0.13，差 4~5 倍）。
///
/// 根因是定义域：线性 0.023 对应 L*≈18、线性 0.119 对应 L*≈42——同一个箱宽在阴影处
/// 的感知跨度比在高光处大一个数量级。**用线性等距的网格去拟合感知上高度不均匀的
/// 映射，阴影段必然欠拟合。**
///
/// 编码域里等距的箱对应的感知跨度大致相当，阴影与高光因此得到相称的分辨率。
pub fn encode_srgb(x: f32) -> f32 {
    let x = x.clamp(0.0, 1.0);
    if x <= 0.003_130_8 {
        12.92 * x
    } else {
        1.055 * x.powf(1.0 / 2.4) - 0.055
    }
}

/// sRGB 编码 → 线性。
pub fn decode_srgb(x: f32) -> f32 {
    let x = x.clamp(0.0, 1.0);
    if x <= 0.040_45 {
        x / 12.92
    } else {
        ((x + 0.055) / 1.055).powf(2.4)
    }
}

/// 曲线的前向施加：**编码 → 查曲线 → 解码**。
///
/// 曲线的定义域与值域都是 0..1 的**编码值**，而管线其余部分走线性，因此这里必须做
/// 两次转换。少了任何一次，曲线都会作用在错误的定义域上——而那种错误不会报错，
/// 只会让画面整体偏亮或偏暗。
pub fn forward_curve(curve: &[f32], v: [f32; 3]) -> [f32; 3] {
    if curve.len() < 2 {
        return v;
    }
    [
        decode_srgb(apply_curve(curve, encode_srgb(v[0]))),
        decode_srgb(apply_curve(curve, encode_srgb(v[1]))),
        decode_srgb(apply_curve(curve, encode_srgb(v[2]))),
    ]
}

/// 拟合一维色调曲线。
///
/// # 做法与取舍
///
/// 按输入值分箱、取箱内输出值的**中位数**——中位数而非均值，是因为它能挡住少量
/// 错配或极端像素。空箱用两侧邻居线性插值补上，最后**强制单调不减**。
///
/// 强制单调的理由：基准渲染本身是单调的，非单调的拟合结果只可能来自噪声。放它过去
/// 会让画面出现局部反转的"色带"，而这种缺陷很难在缩略图上被发现。
///
/// 返回 `bins + 1` 个采样点，定义域与值域都是 0..1 的**编码值**。
pub fn fit_curve(samples: &[Sample], bins: usize) -> Vec<f32> {
    let bins = bins.max(2);
    let mut buckets: Vec<Vec<f32>> = vec![Vec::new(); bins];
    for s in samples {
        for k in 0..3 {
            let x = encode_srgb(s.ours[k]);
            let b = ((x * bins as f32) as usize).min(bins - 1);
            buckets[b].push(encode_srgb(s.theirs[k]));
        }
    }

    // 各箱中位数；空箱记 None
    let mut med: Vec<Option<f32>> = buckets
        .iter_mut()
        .map(|v| {
            if v.is_empty() {
                return None;
            }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            Some(v[v.len() / 2])
        })
        .collect();

    // 空箱用两侧已知点线性插值；两端则用最近点外推
    let known: Vec<usize> = (0..bins).filter(|i| med[*i].is_some()).collect();
    if known.is_empty() {
        // 没有任何样本——返回恒等，调用方应据此判断拟合无效
        return (0..=bins).map(|i| i as f32 / bins as f32).collect();
    }
    for i in 0..bins {
        if med[i].is_some() {
            continue;
        }
        let lo = known.iter().rev().find(|k| **k < i).copied();
        let hi = known.iter().find(|k| **k > i).copied();
        med[i] = match (lo, hi) {
            (Some(l), Some(h)) => {
                let f = (i - l) as f32 / (h - l) as f32;
                Some(med[l].unwrap() * (1.0 - f) + med[h].unwrap() * f)
            }
            (Some(l), None) => med[l],
            (None, Some(h)) => med[h],
            (None, None) => Some(i as f32 / bins as f32),
        };
    }

    // 采样成 bins+1 个点，并强制单调不减
    let mut out: Vec<f32> = (0..=bins)
        .map(|i| {
            if i == bins {
                med[bins - 1].unwrap()
            } else {
                med[i].unwrap()
            }
        })
        .collect();
    for i in 1..out.len() {
        if out[i] < out[i - 1] {
            out[i] = out[i - 1];
        }
    }
    out
}

/// LUT 网格坐标：**编码域**等距。
///
/// # 为什么网格索引必须走编码域
///
/// 曲线早就搬到编码域了，但 LUT 曾经仍按**线性值**分格（`x * (edge-1)`）。后果实测：
/// 17³ = 4913 个格子**只有 220 个有样本，95.52% 是空的**；暗角格 (0,0,0) 吃掉
/// **34.35%** 的样本，其 ΔE00 中位 4.70（全局 2.16 倍），占坏尾部的 **72.29%**。
///
/// 更关键的是分阶段归因（同一评估集，暗部最暗 10%）：
///
/// ```text
/// 恒等 ΔE00 1.61  →  仅曲线 1.87  →  曲线 + LUT 5.56
/// 中位 Y   0.0068  →      0.0075  →           0.0186      参考 0.0072
/// ```
///
/// 即**原始解码与曲线在深阴影里都是对的，是 LUT 把深阴影抬亮了 2.6 倍**——因为整个
/// 暗部被压进同一个格子，而节点值取的是**格内样本的中位数**，那些更亮的样本把节点
/// 抬高，插值再把抬高的节点摊回整个暗部。
///
/// 索引搬到编码域后，暗部样本分散到多个格子，每格内的取值范围大幅收窄，中位数才
/// 真正代表该节点。
fn lut_coord(x: f32, edge: usize) -> f32 {
    encode_srgb(x) * (edge - 1) as f32
}

/// 在采样点上对曲线求值（线性插值）。
pub fn apply_curve(curve: &[f32], x: f32) -> f32 {
    if curve.len() < 2 {
        return x;
    }
    let x = x.clamp(0.0, 1.0);
    let pos = x * (curve.len() - 1) as f32;
    let i = pos.floor() as usize;
    let j = (i + 1).min(curve.len() - 1);
    let f = pos - i as f32;
    curve[i] * (1.0 - f) + curve[j] * f
}

/// 拟合三维色彩查找表。
///
/// # 它在补什么
///
/// 一维曲线只处理**明暗**，补不了通道之间的差异。LUT 接收曲线之后的 RGB，输出参考
/// 的 RGB——因此它承担的是色彩重映射那一半。
///
/// # 空网格怎么填
///
/// 画面分布不均，很多网格没有样本（尤其高光端）。空网格先由**邻域均值**迭代填充，
/// 剩余仍为空的用最近的有效网格填充。这与"留 0"不同：留 0 会在那些区域产生**黑色
/// 斑块**，而它们恰好落在人眼敏感的高光与暗部。
///
/// 返回 `edge³ × 3` 个值，行主序，R 变化最快。
pub fn fit_lut(samples: &[Sample], curve: &[f32], edge: usize) -> Vec<f32> {
    let edge = edge.max(2);
    let n = edge * edge * edge;
    let mut sum = vec![[0f64; 3]; n];
    let mut cnt = vec![0u32; n];

    let idx = |v: [f32; 3]| -> usize {
        let q = |x: f32| (lut_coord(x, edge).round() as usize).min(edge - 1);
        (q(v[0]) * edge + q(v[1])) * edge + q(v[2])
    };

    for s in samples {
        let c = forward_curve(curve, s.ours);
        let i = idx(c);
        for (k, t) in s.theirs.iter().enumerate() {
            sum[i][k] += *t as f64;
        }
        cnt[i] += 1;
    }

    // 已知网格取均值
    let mut lut = vec![[f32::NAN; 3]; n];
    for i in 0..n {
        if cnt[i] > 0 {
            for k in 0..3 {
                lut[i][k] = (sum[i][k] / cnt[i] as f64) as f32;
            }
        }
    }

    // 邻域迭代填充
    for _ in 0..edge * 2 {
        let mut filled = 0;
        let mut next = lut.clone();
        for x in 0..edge {
            for y in 0..edge {
                for z in 0..edge {
                    let i = (x * edge + y) * edge + z;
                    if !lut[i][0].is_nan() {
                        continue;
                    }
                    let mut acc = [0f32; 3];
                    let mut c = 0u32;
                    for (dx, dy, dz) in [
                        (-1i32, 0i32, 0i32),
                        (1, 0, 0),
                        (0, -1, 0),
                        (0, 1, 0),
                        (0, 0, -1),
                        (0, 0, 1),
                    ] {
                        let (nx, ny, nz) = (x as i32 + dx, y as i32 + dy, z as i32 + dz);
                        if nx < 0 || ny < 0 || nz < 0 {
                            continue;
                        }
                        let (nx, ny, nz) = (nx as usize, ny as usize, nz as usize);
                        if nx >= edge || ny >= edge || nz >= edge {
                            continue;
                        }
                        let j = (nx * edge + ny) * edge + nz;
                        if lut[j][0].is_nan() {
                            continue;
                        }
                        for k in 0..3 {
                            acc[k] += lut[j][k];
                        }
                        c += 1;
                    }
                    if c > 0 {
                        for k in 0..3 {
                            next[i][k] = acc[k] / c as f32;
                        }
                        filled += 1;
                    }
                }
            }
        }
        lut = next;
        if filled == 0 {
            break;
        }
    }

    // 仍有空的（极端情况）：取全体已知网格的均值
    let mut avg = [0f32; 3];
    let mut m = 0u32;
    for v in lut.iter() {
        if !v[0].is_nan() {
            for k in 0..3 {
                avg[k] += v[k];
            }
            m += 1;
        }
    }
    let fallback = if m > 0 {
        [avg[0] / m as f32, avg[1] / m as f32, avg[2] / m as f32]
    } else {
        [0.0, 0.0, 0.0]
    };

    let mut out = Vec::with_capacity(n * 3);
    for v in lut {
        for (k, x) in v.iter().enumerate() {
            out.push(if x.is_nan() { fallback[k] } else { *x });
        }
    }
    out
}

/// 在 LUT 上做三线性插值求值。
pub fn apply_lut(lut: &[f32], edge: usize, v: [f32; 3]) -> [f32; 3] {
    if edge < 2 || lut.len() < edge * edge * edge * 3 {
        return v;
    }
    let f = |x: f32| -> (usize, usize, f32) {
        let p = lut_coord(x, edge);
        let i = (p.floor() as usize).min(edge - 2);
        (i, i + 1, p - i as f32)
    };
    let (r0, r1, fr) = f(v[0]);
    let (g0, g1, fg) = f(v[1]);
    let (b0, b1, fb) = f(v[2]);
    let at = |r: usize, g: usize, b: usize, k: usize| lut[((r * edge + g) * edge + b) * 3 + k];

    let mut out = [0f32; 3];
    for (k, slot) in out.iter_mut().enumerate() {
        let c00 = at(r0, g0, b0, k) * (1.0 - fb) + at(r0, g0, b1, k) * fb;
        let c01 = at(r0, g1, b0, k) * (1.0 - fb) + at(r0, g1, b1, k) * fb;
        let c10 = at(r1, g0, b0, k) * (1.0 - fb) + at(r1, g0, b1, k) * fb;
        let c11 = at(r1, g1, b0, k) * (1.0 - fb) + at(r1, g1, b1, k) * fb;
        let c0 = c00 * (1.0 - fg) + c01 * fg;
        let c1 = c10 * (1.0 - fg) + c11 * fg;
        *slot = c0 * (1.0 - fr) + c1 * fr;
    }
    out
}

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
