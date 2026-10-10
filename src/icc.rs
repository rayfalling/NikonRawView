//! ICC 配置文件的解析与转换。
//!
//! # 为什么必须尊重内嵌 ICC，不能假定 sRGB
//!
//! 参考导出是**尼康工坊写的 TIFF**，内嵌一份 `Nikon sRGB 4.0.0.3002`。实测它的
//! 色料矩阵与 sRGB 相同，但**三通道的 TRC 是各 4096 点的采样曲线，不是 sRGB 那条
//! 解析式**。差异就藏在曲线里——当成普通 sRGB 处理会在暗部与高光累积偏差，而这
//! 份偏差会被后续的标定**误当成基准渲染的一部分拟合进去**。
//!
//! 因此这里的做法是：曲线与矩阵**一律从文件里读**，不写死任何一条已知曲线。
//!
//! # 只支持矩阵/TRC 型
//!
//! 本模块只认 `mntr` 类别 + `XYZ ` PCS + `rXYZ/gXYZ/bXYZ` + `rTRC/gTRC/bTRC` 的
//! 组合（ICC 里最常见的一类）。遇到基于 LUT 的 profile 会**明确报错**，而不是
//! 悄悄退回 sRGB——那正是本模块要防的事。

use crate::error::{Error, Result};
use crate::mat3::{self, Mat3};

/// 解析出来的 ICC 配置文件。
#[derive(Debug, Clone, PartialEq)]
pub struct IccProfile {
    /// `desc` 标签里的可读名称。
    pub description: String,
    /// 线性 RGB → PCS XYZ（D50）的色料矩阵。
    pub colorants: Mat3,
    /// 源白点（ICC 头里的 `wtpt`）。
    pub white: [f32; 3],
    /// 三通道的色调曲线（已归一化到 0..1 的定义域与值域）。
    pub trc: [Vec<f32>; 3],
}

const TAG_DESC: u32 = 0x6465_7363; // 'desc'
const TAG_WTPT: u32 = 0x7774_7074; // 'wtpt'
const TAG_RXYZ: u32 = 0x7258_595A; // 'rXYZ'
const TAG_GXYZ: u32 = 0x6758_595A; // 'gXYZ'
const TAG_BXYZ: u32 = 0x6258_595A; // 'bXYZ'
const TAG_RTRC: u32 = 0x7254_5243; // 'rTRC'
const TAG_GTRC: u32 = 0x6754_5243; // 'gTRC'
const TAG_BTRC: u32 = 0x6254_5243; // 'bTRC'

impl IccProfile {
    /// 解析一份 ICC 配置文件。
    pub fn parse(d: &[u8]) -> Result<Self> {
        if d.len() < 132 {
            return Err(Error::BadPayload { what: "ICC 配置", detail: format!("长度仅 {}", d.len()) });
        }
        if &d[36..40] != b"acsp" {
            return Err(Error::BadPayload { what: "ICC 配置", detail: "magic 不是 acsp".into() });
        }
        if &d[16..20] != b"RGB " {
            return Err(Error::BadPayload {
                what: "ICC 配置",
                detail: format!("数据色彩空间不是 RGB，而是 {:?}", String::from_utf8_lossy(&d[16..20])),
            });
        }
        if &d[20..24] != b"XYZ " {
            return Err(Error::BadPayload {
                what: "ICC 配置",
                detail: format!("PCS 不是 XYZ，而是 {:?}", String::from_utf8_lossy(&d[20..24])),
            });
        }

        let count = be32(d, 128)? as usize;
        let mut tags = std::collections::BTreeMap::new();
        for i in 0..count {
            let p = 132 + i * 12;
            if p + 12 > d.len() {
                return Err(Error::BadPayload { what: "ICC 配置", detail: "标签表越界".into() });
            }
            let sig = be32(d, p)?;
            let off = be32(d, p + 4)? as usize;
            let size = be32(d, p + 8)? as usize;
            if off.saturating_add(size) > d.len() {
                return Err(Error::BadPayload {
                    what: "ICC 配置",
                    detail: format!("标签 {sig:#010x} 越界"),
                });
            }
            tags.insert(sig, (off, size));
        }

        let need = |sig: u32, what: &'static str| -> Result<(usize, usize)> {
            tags.get(&sig).copied().ok_or_else(|| Error::BadPayload {
                what: "ICC 配置",
                detail: format!("缺少{what}标签——本模块只支持矩阵/TRC 型 profile，不退回 sRGB"),
            })
        };

        let description = tags
            .get(&TAG_DESC)
            .and_then(|(o, s)| read_desc(d, *o, *s).ok())
            .unwrap_or_else(|| "（无描述）".into());

        let colorants: Mat3 = [
            read_xyz(d, need(TAG_RXYZ, " rXYZ")?)?,
            read_xyz(d, need(TAG_GXYZ, " gXYZ")?)?,
            read_xyz(d, need(TAG_BXYZ, " bXYZ")?)?,
        ];
        // 注意：色料矩阵的行是 r/g/b 的 XYZ，构矩阵时按 **列** 放（XYZ = M · rgb）
        let m: Mat3 = [
            [colorants[0][0], colorants[1][0], colorants[2][0]],
            [colorants[0][1], colorants[1][1], colorants[2][1]],
            [colorants[0][2], colorants[1][2], colorants[2][2]],
        ];

        let white = read_xyz(d, need(TAG_WTPT, " wtpt")?)?;

        let trc = [
            read_curv(d, need(TAG_RTRC, " rTRC")?)?,
            read_curv(d, need(TAG_GTRC, " gTRC")?)?,
            read_curv(d, need(TAG_BTRC, " bTRC")?)?,
        ];

        Ok(Self { description, colorants: m, white, trc })
    }

    /// 去掉曲线的编码，得到线性光值。
    pub fn linearize(&self, rgb: [f32; 3]) -> [f32; 3] {
        [
            apply_curve(&self.trc[0], rgb[0]),
            apply_curve(&self.trc[1], rgb[1]),
            apply_curve(&self.trc[2], rgb[2]),
        ]
    }

    /// 把该 profile 描述的颜色转换到 **ProPhoto 线性（D50）**——本项目的工作空间。
    pub fn to_working_space(&self, rgb: [f32; 3]) -> [f32; 3] {
        let lin = self.linearize(rgb);
        let xyz = mat3::mul_vec(self.colorants, lin);
        mat3::mul_vec(mat3::XYZ_TO_PROPHOTO, xyz)
    }
}

/// 参考导出的色彩空间处理方案。
///
/// 存在的意义是**让"假定"必须被写下来**：没有内嵌 ICC 时不能默默按 sRGB 处理，
/// 而要把这个假定连同理由带进后续的报告里——否则一个静默的假定会一路混进标定，
/// 直到最后也没人知道它存在。
#[derive(Debug, Clone, PartialEq)]
pub enum ColorSpacePlan {
    /// 按文件内嵌的 ICC 转换。
    Embedded(Box<IccProfile>),
    /// 文件没有内嵌 ICC，按 sRGB 假定处理——**这是一个假定，不是事实**。
    AssumedSrgb { reason: String },
}

impl ColorSpacePlan {
    /// 该方案的可读说明，供报告与日志使用。
    pub fn describe(&self) -> String {
        match self {
            ColorSpacePlan::Embedded(p) => format!("按内嵌 ICC「{}」转换", p.description),
            ColorSpacePlan::AssumedSrgb { reason } => {
                format!("**假定 sRGB**（{reason}）——该文件未内嵌 ICC，此结果依赖假定")
            }
        }
    }

    /// 是否为假定。
    pub fn is_assumption(&self) -> bool {
        matches!(self, ColorSpacePlan::AssumedSrgb { .. })
    }

    /// 把某色彩空间下的 RGB 转到工作空间（ProPhoto 线性）。
    ///
    /// 假定 sRGB 的分支走的是标准 sRGB 解析式——那正是"假定"二字的含义。
    pub fn to_working_space(&self, rgb: [f32; 3]) -> [f32; 3] {
        match self {
            ColorSpacePlan::Embedded(p) => p.to_working_space(rgb),
            ColorSpacePlan::AssumedSrgb { .. } => {
                let lin = rgb.map(|x| {
                    let x = x.clamp(0.0, 1.0);
                    if x <= 0.04045 {
                        x / 12.92
                    } else {
                        ((x + 0.055) / 1.055).powf(2.4)
                    }
                });
                mat3::mul_vec(
                    mat3::XYZ_TO_PROPHOTO,
                    mat3::mul_vec(SRGB_TO_XYZ_D50, lin),
                )
            }
        }
    }
}

/// sRGB 的线性 RGB → PCS XYZ（D50）矩阵。仅用于「缺 ICC 时假定 sRGB」这一条路径。
const SRGB_TO_XYZ_D50: Mat3 = [
    [0.436_07, 0.385_07, 0.143_05],
    [0.222_50, 0.716_87, 0.060_61],
    [0.013_92, 0.097_06, 0.713_99],
];

/// 为一份参考导出决定色彩空间处理方案。
///
/// 有内嵌 ICC 就用它；没有则**明确记为假定**并把理由带出来。
pub fn plan_for(tiff: &[u8]) -> Result<ColorSpacePlan> {
    match embedded_profile(tiff)? {
        Some(raw) => Ok(ColorSpacePlan::Embedded(Box::new(IccProfile::parse(&raw)?))),
        None => Ok(ColorSpacePlan::AssumedSrgb {
            reason: "文件未内嵌 ICC 配置".into(),
        }),
    }
}

/// 从 TIFF 数据里取出内嵌的 ICC 配置（IFD0 标签 `0x8773`）。
///
/// 返回 `None` 表示该文件没有内嵌配置——调用方**必须明确报告并标注假定**，
/// 不得默默按 sRGB 处理。
pub fn embedded_profile(tiff: &[u8]) -> Result<Option<Vec<u8>>> {
    use crate::tiff::Tiff;
    let t = Tiff::locate(tiff)?;
    let ifd = t.ifd0_offset()?;
    match t.find(ifd, 0x8773)? {
        Some(e) => Ok(Some(t.bytes(&e)?.to_vec())),
        None => Ok(None),
    }
}

fn be32(d: &[u8], at: usize) -> Result<u32> {
    let b = d.get(at..at + 4).ok_or_else(|| Error::BadPayload {
        what: "ICC 配置",
        detail: format!("读取 4 字节时越界（偏移 {at}）"),
    })?;
    Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
}

fn s15f16(d: &[u8], at: usize) -> Result<f32> {
    let v = be32(d, at)? as i32;
    Ok(v as f32 / 65536.0)
}

/// 读一个 `XYZ ` 类型标签。
fn read_xyz(d: &[u8], (o, s): (usize, usize)) -> Result<[f32; 3]> {
    if s < 20 || &d[o..o + 4] != b"XYZ " {
        return Err(Error::BadPayload {
            what: "ICC 配置",
            detail: format!("XYZ 标签格式不符（长度 {s}）"),
        });
    }
    Ok([s15f16(d, o + 8)?, s15f16(d, o + 12)?, s15f16(d, o + 16)?])
}

/// 读一个 `curv` 类型标签，归一化为 0..1 → 0..1 的采样表。
fn read_curv(d: &[u8], (o, s): (usize, usize)) -> Result<Vec<f32>> {
    if s < 12 || &d[o..o + 4] != b"curv" {
        return Err(Error::BadPayload { what: "ICC 配置", detail: "TRC 标签格式不符".into() });
    }
    let n = be32(d, o + 8)? as usize;
    if n == 0 {
        // 恒等
        return Ok(vec![0.0, 1.0]);
    }
    if n == 1 {
        let g = u16::from_be_bytes([d[o + 12], d[o + 13]]) as f32 / 256.0;
        // 采样成表，后续统一按表插值
        return Ok((0..=1024).map(|i| (i as f32 / 1024.0).powf(g)).collect());
    }
    if o + 12 + n * 2 > d.len() {
        return Err(Error::BadPayload { what: "ICC 配置", detail: "TRC 采样越界".into() });
    }
    Ok((0..n)
        .map(|i| {
            let v = u16::from_be_bytes([d[o + 12 + i * 2], d[o + 13 + i * 2]]);
            v as f32 / 65535.0
        })
        .collect())
}

/// 在采样表上做线性插值。
fn apply_curve(t: &[f32], x: f32) -> f32 {
    if t.len() < 2 {
        return x;
    }
    let x = x.clamp(0.0, 1.0);
    let pos = x * (t.len() - 1) as f32;
    let i = pos.floor() as usize;
    let j = (i + 1).min(t.len() - 1);
    let f = pos - i as f32;
    t[i] * (1.0 - f) + t[j] * f
}

/// 读 `desc`（ICC v2 的 `desc` 或 v4 的 `mluc`）。
fn read_desc(d: &[u8], o: usize, s: usize) -> Result<String> {
    if s < 12 {
        return Err(Error::BadPayload { what: "ICC 配置", detail: "desc 过短".into() });
    }
    let typ = &d[o..o + 4];
    if typ == b"desc" {
        let n = be32(d, o + 8)? as usize;
        let end = (o + 12 + n).min(d.len());
        return Ok(String::from_utf8_lossy(&d[o + 12..end]).trim_end_matches('\0').to_string());
    }
    if typ == b"mluc" {
        let cnt = be32(d, o + 8)? as usize;
        if cnt > 0 {
            let len = be32(d, o + 16)? as usize;
            let off = be32(d, o + 20)? as usize;
            let st = o + off;
            if st + len <= d.len() {
                let u16s: Vec<u16> = d[st..st + len]
                    .chunks_exact(2)
                    .map(|c| u16::from_be_bytes([c[0], c[1]]))
                    .collect();
                return Ok(String::from_utf16_lossy(&u16s));
            }
        }
    }
    Err(Error::BadPayload {
        what: "ICC 配置",
        detail: format!("不认识的 desc 类型 {:?}", String::from_utf8_lossy(typ)),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("simple/DSC_0001.TIF")
    }

    /// 从 TIFF 里取出内嵌 ICC。
    fn real_icc() -> Option<Vec<u8>> {
        let p = sample_path();
        let d = std::fs::read(&p).ok()?;
        embedded_profile(&d).ok().flatten()
    }

    #[test]
    fn parses_the_real_nikon_profile() {
        let Some(icc) = real_icc() else {
            eprintln!("跳过：simple/DSC_0001.TIF 或其内嵌 ICC 不存在");
            return;
        };
        let p = IccProfile::parse(&icc).expect("应能解析");
        eprintln!("desc = {}", p.description);
        eprintln!("white = {:?}", p.white);
        for row in &p.colorants {
            eprintln!("  colorant row [{:.5} {:.5} {:.5}]", row[0], row[1], row[2]);
        }
        eprintln!("TRC 长度 = {:?}", p.trc.iter().map(Vec::len).collect::<Vec<_>>());

        assert!(p.description.contains("Nikon"), "描述应含 Nikon：{}", p.description);
        assert!(
            p.trc.iter().all(|t| t.len() > 100),
            "TRC 应是采样曲线而非单点 gamma"
        );
        // 白点接近 D65（sRGB 系）
        assert!((p.white[1] - 1.0).abs() < 1e-3, "白点 Y 应归一化：{:?}", p.white);
    }

    #[test]
    fn nikon_srgb_profile_matches_analytic_srgb() {
        // 实测结论：`Nikon sRGB 4.0.0.3002` **就是**标准 sRGB。
        //
        // 它的色料矩阵与 sRGB 相同，4096 点采样曲线与 sRGB 解析式之差仅
        // 1.5e-5，而 16 位量化步长是 1/65535 ≈ 1.53e-5——差异完全来自采样量化。
        //
        // 这条断言的意义不是"两者应当相同"，而是**记录这个事实**：它与
        // `trc_is_read_from_file_not_assumed` 一起，把"读文件"与"假定 sRGB"
        // 这两种实现区分开——前者对任意 profile 都成立，后者只在这一个上凑巧成立。
        let Some(icc) = real_icc() else {
            eprintln!("跳过：无内嵌 ICC");
            return;
        };
        let p = IccProfile::parse(&icc).unwrap();
        let srgb = |x: f32| {
            if x <= 0.04045 {
                x / 12.92
            } else {
                ((x + 0.055) / 1.055).powf(2.4)
            }
        };

        let mut worst = 0f32;
        for i in 0..=1000 {
            let x = i as f32 / 1000.0;
            worst = worst.max((p.linearize([x, x, x])[0] - srgb(x)).abs());
        }
        eprintln!("内嵌 TRC 与 sRGB 解析式的最大差异 = {worst:.8}（量化步长 1/65535 = {:.8}）", 1.0 / 65535.0);
        assert!(
            worst < 2.0 / 65535.0,
            "差异 {worst:.8} 超出量化范围——要么实现读错了曲线，要么这份 profile 确实不是 sRGB"
        );
    }

    #[test]
    fn trc_is_read_from_file_not_assumed() {
        // 这条才真正锁死「必须读文件里的曲线」：造一份 gamma=1.8 的 profile，
        // 若实现退回 sRGB 解析式，结果就会与它一致，断言随即失败。
        let icc = synthetic_profile(1.8);
        let p = IccProfile::parse(&icc).expect("应能解析合成 profile");
        let srgb = |x: f32| {
            if x <= 0.04045 {
                x / 12.92
            } else {
                ((x + 0.055) / 1.055).powf(2.4)
            }
        };
        let x = 0.5f32;
        let ours = p.linearize([x, x, x])[0];
        let gamma18 = 0.5f32.powf(1.8);
        eprintln!("gamma1.8 profile: 本实现 {ours:.6}，gamma1.8 期望 {gamma18:.6}，sRGB {:.6}", srgb(x));
        assert!((ours - gamma18).abs() < 1e-3, "应复现文件里的 gamma 1.8");
        assert!(
            (ours - srgb(x)).abs() > 0.05,
            "与 sRGB 差异过小——实现可能没读曲线"
        );
    }

    /// 造一份最小的矩阵/TRC 型 ICC，三通道都用 gamma 曲线。
    fn synthetic_profile(gamma: f32) -> Vec<u8> {
        let mut d = vec![0u8; 512];
        d[8..12].copy_from_slice(&0x0220_0000u32.to_be_bytes());
        d[16..20].copy_from_slice(b"RGB ");
        d[20..24].copy_from_slice(b"XYZ ");
        d[36..40].copy_from_slice(b"acsp");

        let col = [
            (0.43607f32, 0.22250, 0.01392),
            (0.38507, 0.71687, 0.09706),
            (0.14305, 0.06061, 0.71399),
        ];
        // 数据区必须从标签表之后开始：表头 4 字节 + 7 个 12 字节条目 = 216`n        let mut off = 128 + 4 + 7 * 12;
        let put_xyz = |d: &mut Vec<u8>, o: usize, v: (f32, f32, f32)| {
            d[o..o + 4].copy_from_slice(b"XYZ ");
            d[o + 4..o + 8].copy_from_slice(&0u32.to_be_bytes());
            for (i, x) in [v.0, v.1, v.2].iter().enumerate() {
                let raw = (x * 65536.0).round() as i32;
                d[o + 8 + i * 4..o + 12 + i * 4].copy_from_slice(&raw.to_be_bytes());
            }
        };
        let put_curv = |d: &mut Vec<u8>, o: usize, g: f32| {
            d[o..o + 4].copy_from_slice(b"curv");
            d[o + 4..o + 8].copy_from_slice(&0u32.to_be_bytes());
            d[o + 8..o + 12].copy_from_slice(&1u32.to_be_bytes());
            let g16 = (g * 256.0).round() as u16;
            d[o + 12..o + 14].copy_from_slice(&g16.to_be_bytes());
        };

        let entry_off = 128 + 4;
        let sigs: [(u32, usize); 4] = [
            (TAG_RXYZ, 20),
            (TAG_GXYZ, 20),
            (TAG_BXYZ, 20),
            (TAG_WTPT, 20),
        ];
        d[128..132].copy_from_slice(&7u32.to_be_bytes());
        // 数据区必须从标签表**之后**开始：表头 4 字节 + 7 个 12 字节条目 = 216。
        // 早前算成 188，导致写标签条目时把色料数据覆盖掉了。
        let mut off = 128 + 4 + 7 * 12;
        for (i, (sig, size)) in sigs.iter().enumerate() {
            let p = entry_off + i * 12;
            d[p..p + 4].copy_from_slice(&sig.to_be_bytes());
            d[p + 4..p + 8].copy_from_slice(&(off as u32).to_be_bytes());
            d[p + 8..p + 12].copy_from_slice(&(*size as u32).to_be_bytes());
            let v = match *sig {
                TAG_RXYZ => col[0],
                TAG_GXYZ => col[1],
                TAG_BXYZ => col[2],
                _ => (0.9642, 1.0, 0.8249),
            };
            put_xyz(&mut d, off, v);
            off += size;
        }
        for (i, sig) in [TAG_RTRC, TAG_GTRC, TAG_BTRC].iter().enumerate() {
            let p = entry_off + (4 + i) * 12;
            d[p..p + 4].copy_from_slice(&sig.to_be_bytes());
            d[p + 4..p + 8].copy_from_slice(&(off as u32).to_be_bytes());
            d[p + 8..p + 12].copy_from_slice(&14u32.to_be_bytes());
            put_curv(&mut d, off, gamma);
            off += 14;
        }
        d.truncate(off);
        d
    }

    #[test]
    fn missing_icc_becomes_an_explicit_assumption() {
        // 造一份没有 0x8773 标签的最小 TIFF
        let mut d = vec![0u8; 32];
        d[0..2].copy_from_slice(b"II");
        d[2..4].copy_from_slice(&42u16.to_le_bytes());
        d[4..8].copy_from_slice(&8u32.to_le_bytes());
        d[8..10].copy_from_slice(&0u16.to_le_bytes()); // IFD0 条目数 0
        let plan = plan_for(&d).expect("无 ICC 不应报错，而应记为假定");
        assert!(plan.is_assumption(), "应标注为假定");
        let desc = plan.describe();
        eprintln!("{desc}");
        assert!(desc.contains("假定"), "说明中必须写明是假定：{desc}");
        assert!(desc.contains("未内嵌 ICC"));
    }

    #[test]
    fn real_export_uses_embedded_profile_not_an_assumption() {
        let p = sample_path();
        let Ok(d) = std::fs::read(&p) else {
            eprintln!("跳过：无参考导出样本");
            return;
        };
        let plan = plan_for(&d).expect("应能决定方案");
        eprintln!("{}", plan.describe());
        assert!(!plan.is_assumption(), "真实导出内嵌了 ICC，不应走假定分支");
        match plan {
            ColorSpacePlan::Embedded(prof) => {
                assert!(prof.description.contains("Nikon"), "{}", prof.description)
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn rejects_non_matrix_profiles() {
        let mut d = vec![0u8; 200];
        d[36..40].copy_from_slice(b"acsp");
        d[16..20].copy_from_slice(b"RGB ");
        d[20..24].copy_from_slice(b"XYZ ");
        d[128..132].copy_from_slice(&0u32.to_be_bytes());
        // 缺 rXYZ → 必须明确报错，不得退回 sRGB
        let e = IccProfile::parse(&d);
        assert!(matches!(e, Err(Error::BadPayload { .. })), "{e:?}");
        assert!(e.unwrap_err().to_string().contains("不退回 sRGB"));
    }

    #[test]
    fn rejects_lut_based_pcs() {
        let mut d = vec![0u8; 200];
        d[36..40].copy_from_slice(b"acsp");
        d[16..20].copy_from_slice(b"RGB ");
        d[20..24].copy_from_slice(b"Lab ");
        let e = IccProfile::parse(&d);
        assert!(e.is_err());
    }

    #[test]
    fn rejects_bad_magic_and_truncation() {
        assert!(IccProfile::parse(&[0u8; 8]).is_err());
        let mut d = vec![0u8; 300];
        d[36..40].copy_from_slice(b"junk");
        assert!(IccProfile::parse(&d).is_err());
    }

    #[test]
    fn curve_is_monotonic_and_anchored() {
        let t: Vec<f32> = (0..=100).map(|i| (i as f32 / 100.0).powf(2.2)).collect();
        assert!(apply_curve(&t, 0.0).abs() < 1e-6);
        assert!((apply_curve(&t, 1.0) - 1.0).abs() < 1e-6);
        assert!(apply_curve(&t, 0.5) < 0.5, "gamma 2.2 在中点应低于线性");
        // 夹取
        assert!(apply_curve(&t, -1.0).abs() < 1e-6);
        assert!((apply_curve(&t, 2.0) - 1.0).abs() < 1e-6);
    }
}
