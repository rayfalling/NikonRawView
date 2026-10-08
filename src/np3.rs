//! NP3 / NCP 配方容器解析。
//!
//! # 容器结构（实测确认）
//!
//! ```text
//! 0x00  "NCP\0"                       magic（NP3 与 NCP 同族）
//! 0x04  u32BE 版本（0x00000100）
//! 0x08  u32BE 4
//! 0x0C  4 字节家族码（"0310" 或 "0300"）
//! 0x10  chunk 序列：[u32BE tag][u32BE len][payload]，至 tag=0 && len=0 结束
//! ```
//!
//! # 为什么按 tag 而非固定偏移
//!
//! 公开的 NP3 布局实现基于 392 字节固定模板，其曲线块固定在 `0x1CC`。实测
//! 本机 `0300` 代文件的曲线块落在 `0xF4` 与 `0x13A`——按固定偏移会**读到
//! 错误数据且不报错**，静默产出看似合理的错误曲线。因此本实现一律按 tag 定位。
//!
//! 标量同样按 tag 读取：其数值是 tag 的函数，而非文件中的固定位置。
//! 实测 `06 00` ↔ 锐化（0x52）、`07 00` ↔ 清晰度（0x5C），与公开布局一致。

use std::collections::BTreeMap;

use crate::error::{Error, Result};

pub const MAGIC: &[u8; 4] = b"NCP\0";

/// 具名 chunk 标签。
pub mod tag {
    /// 曲线块。
    pub const CURVE: u32 = 0x0000_0002;
    /// 配方显示名（20 字节）。
    pub const NAME: u32 = 0x0000_0200;
    /// 基准色彩编码（2 字节）。
    pub const BASE: u32 = 0x0000_0300;
    /// 语义未定（2 字节）。
    pub const UNKNOWN_0400: u32 = 0x0000_0400;
    /// 色相混合器（28 字节 = 8 段 × 3 值）。
    pub const BLENDER: u32 = 0x0000_1F00;
    /// 校色（20 字节）。
    pub const GRADING: u32 = 0x0000_2000;
    /// 首个标量 tag。
    pub const SCALAR_FIRST: u32 = 0x0000_0500;
    /// 末个标量 tag。
    pub const SCALAR_LAST: u32 = 0x0000_1E00;
}

/// 标量字段名。tag 与字段的对应由实测确定。
pub fn scalar_name(t: u32) -> Option<&'static str> {
    Some(match t {
        0x0000_0500 => "unknown_0500",
        0x0000_0600 => "sharpening",
        0x0000_0700 => "clarity",
        0x0000_1600 => "mid_range_sharpening",
        0x0000_1900 => "contrast",
        0x0000_1A00 => "highlights",
        0x0000_1B00 => "shadows",
        0x0000_1C00 => "white_level",
        0x0000_1D00 => "black_level",
        0x0000_1E00 => "saturation",
        t if (tag::SCALAR_FIRST..=tag::SCALAR_LAST).contains(&t) => "scalar",
        _ => return None,
    })
}

/// 一个原始 chunk，未识别者也保留。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chunk {
    pub tag: u32,
    pub offset: usize,
    pub payload: Vec<u8>,
}

/// 色相混合器的一段。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Band {
    pub hue: i8,
    pub chroma: i8,
    pub brightness: i8,
}

/// 混合器的 8 段，顺序固定。
pub const BAND_ORDER: [&str; 8] = ["Red", "Orange", "Yellow", "Green", "Cyan", "Blue", "Purple", "Magenta"];

/// 校色的一个区。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Zone {
    pub hue: u16,
    pub chroma: i8,
    pub brightness: i8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Grading {
    pub highlights: Zone,
    pub mid_tone: Zone,
    pub shadows: Zone,
    pub blending: u8,
    pub balance: u8,
}

/// 色调曲线块的内容。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Curve {
    /// 控制点数量（载荷 +8）。
    pub point_count: u8,
    pub points: Vec<(u8, u8)>,
    /// 257 级 16 位大端 LUT。
    pub lut: Vec<u16>,
    /// LUT 在曲线载荷内的偏移——两代模板在此处相差 2 字节，故不写死。
    pub lut_offset: usize,
}

impl Curve {
    /// LUT 是否严格单调递增。
    pub fn is_monotonic(&self) -> bool {
        self.lut.windows(2).all(|w| w[0] < w[1])
    }
}

/// 解析结果。
#[derive(Debug, Clone)]
pub struct Container {
    pub version: u32,
    /// 家族码，`"0310"` 或 `"0300"`。
    pub family: String,
    pub name: String,
    pub base_code: Option<u16>,
    /// tag → **值字节**，按 tag 升序。
    ///
    /// 标量载荷为 2 字节 `[值][步进/标志]`，值取**首字节**（实测步进字节恒为
    /// `0x04`）。原始两字节仍保留在 `chunks` 中。
    pub scalars: BTreeMap<u32, u8>,
    pub blender: Option<[Band; 8]>,
    pub grading: Option<Grading>,
    pub curve: Option<Curve>,
    /// 全部 chunk（含未识别者），供前向兼容与诊断。
    pub chunks: Vec<Chunk>,
}

/// 标量的解码方式。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scale {
    /// `(b − 0x80) / 4`，四分之一步进。用于锐化、清晰度、中频锐化。
    Quarter,
    /// `b − 0x80`。用于对比度、亮度、饱和度等。
    Linear,
}

impl Container {
    /// 是否含自定义色调曲线。
    pub fn has_custom_curve(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| !is_flat(c))
    }

    /// 取某个标量标签的原始值字节。
    pub fn scalar_raw(&self, tag: u32) -> Option<u8> {
        self.scalars.get(&tag).copied()
    }

    /// 取某个标量标签按指定方式解码后的值。
    pub fn scalar(&self, tag: u32, scale: Scale) -> Option<f32> {
        self.scalar_raw(tag).map(|b| match scale {
            Scale::Quarter => (b as f32 - 128.0) / 4.0,
            Scale::Linear => b as f32 - 128.0,
        })
    }

    /// 锐化（四分之一步进）。
    pub fn sharpening(&self) -> Option<f32> {
        self.scalar(0x0000_0600, Scale::Quarter)
    }
    /// 清晰度（四分之一步进）。
    pub fn clarity(&self) -> Option<f32> {
        self.scalar(0x0000_0700, Scale::Quarter)
    }
    /// 中频锐化（四分之一步进）。
    pub fn mid_range_sharpening(&self) -> Option<f32> {
        self.scalar(0x0000_1600, Scale::Quarter)
    }
    /// 对比度。
    pub fn contrast(&self) -> Option<f32> {
        self.scalar(0x0000_1900, Scale::Linear)
    }
    /// 饱和度。
    pub fn saturation(&self) -> Option<f32> {
        self.scalar(0x0000_1E00, Scale::Linear)
    }
}

fn is_flat(c: &Curve) -> bool {
    // 恒等曲线：lut[i] ≈ i * 32767 / 256
    c.lut.iter().enumerate().all(|(i, &v)| {
        let expect = ((i as u32 * 32767) / 256) as i32;
        (v as i32 - expect).abs() <= 1
    })
}

/// 遍历 chunk 序列。
///
/// 结束符是**单个 u32 = 0**（实测：834 字节的 `0300` 文件在曲线块后只剩 4 字节，
/// 978 字节的 `0310` 文件同样只剩 4 字节）。因此先读 tag，为 0 即停止——按
/// "8 字节全零"处理会把结束符读成缺斤少两的块头而误报截断。
pub fn chunks(data: &[u8]) -> Result<Vec<Chunk>> {
    if data.len() < 0x10 {
        return Err(Error::Truncated { what: "容器头", need: 0x10, have: data.len() });
    }
    if &data[0..4] != MAGIC {
        return Err(Error::NotRecipeContainer([data[0], data[1], data[2], data[3]]));
    }
    let mut out = Vec::new();
    let mut p = 0x10usize;
    // 防跑飞：每个 chunk 至少 8 字节，故上限为剩余长度 / 8
    let guard = data.len() / 8 + 1;
    for _ in 0..guard {
        let tag_raw = data.get(p..p + 4).ok_or(Error::Truncated {
            what: "chunk 标签",
            need: 4,
            have: data.len().saturating_sub(p),
        })?;
        let t = u32::from_be_bytes([tag_raw[0], tag_raw[1], tag_raw[2], tag_raw[3]]);
        if t == 0 {
            return Ok(out);
        }
        let len_raw = data.get(p + 4..p + 8).ok_or(Error::Truncated {
            what: "chunk 长度",
            need: 4,
            have: data.len().saturating_sub(p + 4),
        })?;
        let len = u32::from_be_bytes([len_raw[0], len_raw[1], len_raw[2], len_raw[3]]) as usize;
        let start = p + 8;
        let payload = data.get(start..start + len).ok_or(Error::Truncated {
            what: "chunk 载荷",
            need: len,
            have: data.len().saturating_sub(start),
        })?;
        out.push(Chunk { tag: t, offset: p, payload: payload.to_vec() });
        p = start + len;
    }
    Err(Error::Truncated { what: "chunk 终止符", need: 4, have: 0 })
}

fn fixed_ascii(b: &[u8]) -> String {
    let end = b.iter().position(|&c| c == 0).unwrap_or(b.len());
    String::from_utf8_lossy(&b[..end]).trim().to_string()
}

fn u16_be(b: &[u8]) -> u16 {
    u16::from_be_bytes([b[0], b[1]])
}

/// 在曲线载荷内定位 257 级 LUT。
///
/// 不写死偏移：`0310` 代在载荷 +64，`0300` 代在 +62。改为**先试 +64，再扫描**，
/// 两者都以「257 个值严格递增、起点接近 0、终点接近满量程」为接受条件，因此
/// 不会静默接受错误数据。
fn locate_lut(payload: &[u8]) -> Option<(usize, Vec<u16>)> {
    const N: usize = 257;
    let accept = |off: usize| -> Option<Vec<u16>> {
        let bytes = payload.get(off..off + N * 2)?;
        let v: Vec<u16> = (0..N).map(|i| u16_be(&bytes[i * 2..i * 2 + 2])).collect();
        let inc = v.windows(2).all(|w| w[0] < w[1]);
        let starts_low = v[0] < 0x1000;
        let ends_high = v[N - 1] > 0x7000;
        (inc && starts_low && ends_high).then_some(v)
    };
    if let Some(v) = accept(64) {
        return Some((64, v));
    }
    let hi = payload.len().checked_sub(N * 2)?;
    (0..=hi).find_map(|off| accept(off).map(|v| (off, v)))
}

fn parse_curve(payload: &[u8]) -> Option<Curve> {
    let (lut_offset, lut) = locate_lut(payload)?;
    let point_count = *payload.get(8)?;
    let mut points = Vec::new();
    for i in 0..point_count as usize {
        let p = 9 + i * 2;
        if p + 2 > payload.len() {
            break;
        }
        points.push((payload[p], payload[p + 1]));
    }
    Some(Curve { point_count, points, lut, lut_offset })
}

/// 有符号单字节字段的解码：`b − 0x80`。
///
/// **不能直接 `as i8`**——中性值在文件里存为 `0x80`，强转会得到 −128 而不是 0。
fn s8(b: u8) -> i8 {
    (b as i16 - 128) as i8
}

fn parse_blender(p: &[u8]) -> Option<[Band; 8]> {
    if p.len() < 24 {
        return None;
    }
    let mut out = [Band { hue: 0, chroma: 0, brightness: 0 }; 8];
    for (i, b) in out.iter_mut().enumerate() {
        let o = i * 3;
        *b = Band { hue: s8(p[o]), chroma: s8(p[o + 1]), brightness: s8(p[o + 2]) };
    }
    Some(out)
}

fn parse_zone(p: &[u8], o: usize) -> Zone {
    let hue = (((p[o] & 0x0F) as u16) << 8) | p[o + 1] as u16;
    Zone { hue, chroma: s8(p[o + 2]), brightness: s8(p[o + 3]) }
}

fn parse_grading(p: &[u8]) -> Option<Grading> {
    if p.len() < 20 {
        return None;
    }
    Some(Grading {
        highlights: parse_zone(p, 0),
        mid_tone: parse_zone(p, 4),
        shadows: parse_zone(p, 8),
        blending: p[16],
        balance: p[18],
    })
}

/// 解析一份配方容器。
pub fn parse(data: &[u8]) -> Result<Container> {
    let version = if data.len() >= 8 {
        u32::from_be_bytes([data[4], data[5], data[6], data[7]])
    } else {
        0
    };
    let family = if data.len() >= 0x10 { fixed_ascii(&data[0x0C..0x10]) } else { String::new() };

    let all = chunks(data)?;
    let mut name = String::new();
    let mut base_code = None;
    let mut scalars = BTreeMap::new();
    let mut blender = None;
    let mut grading = None;
    let mut curve = None;

    for c in &all {
        match c.tag {
            tag::NAME => name = fixed_ascii(&c.payload),
            tag::BASE if c.payload.len() >= 2 => base_code = Some(u16_be(&c.payload)),
            tag::CURVE => curve = parse_curve(&c.payload),
            tag::BLENDER => blender = parse_blender(&c.payload),
            tag::GRADING => grading = parse_grading(&c.payload),
            t if (tag::SCALAR_FIRST..=tag::SCALAR_LAST).contains(&t)
                && c.payload.len() >= 2 =>
            {
                scalars.insert(t, c.payload[0]);
            }
            _ => {}
        }
    }

    Ok(Container { version, family, name, base_code, scalars, blender, grading, curve, chunks: all })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 合成一份 `0310` 代容器：name + base + 若干标量 + 曲线块。
    fn synth(family: &[u8; 4], name: &str, base: u16, lut_start: usize) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&0x0000_0100u32.to_be_bytes());
        out.extend_from_slice(&4u32.to_be_bytes());
        out.extend_from_slice(family);

        let push = |out: &mut Vec<u8>, t: u32, payload: &[u8]| {
            out.extend_from_slice(&t.to_be_bytes());
            out.extend_from_slice(&(payload.len() as u32).to_be_bytes());
            out.extend_from_slice(payload);
        };
        let mut nm = vec![0u8; 20];
        nm[..name.len()].copy_from_slice(name.as_bytes());
        push(&mut out, tag::NAME, &nm);
        push(&mut out, tag::BASE, &base.to_be_bytes());
        push(&mut out, tag::UNKNOWN_0400, &[0, 0]);
        push(&mut out, 0x0000_0600, &[0x82, 0x04]); // sharpening +0.5
        push(&mut out, 0x0000_0700, &[0x82, 0x04]); // clarity +0.5
        push(&mut out, 0x0000_1900, &[0x01, 0x04]); // contrast = sentinel

        // 曲线载荷：header + 257×u16BE，LUT 起始位置可调
        let mut cp = vec![0u8; 578];
        cp[0] = b'I';
        cp[1] = b'0';
        cp[3] = 0xFF;
        cp[8] = 5; // 控制点数
        for i in 0..5usize {
            cp[9 + i * 2] = (i * 50) as u8;
            cp[10 + i * 2] = (i * 60) as u8;
        }
        for i in 0..257usize {
            let v = ((i as u32 * 32767) / 256) as u16;
            cp[lut_start + i * 2] = (v >> 8) as u8;
            cp[lut_start + i * 2 + 1] = v as u8;
        }
        push(&mut out, tag::CURVE, &cp);
        out.extend_from_slice(&0u32.to_be_bytes());
        out.extend_from_slice(&0u32.to_be_bytes());
        out
    }

    #[test]
    fn parses_0310_template() {
        let d = synth(b"0310", "LINKS-Nature", 0x03C2, 64);
        let c = parse(&d).unwrap();
        assert_eq!(c.family, "0310");
        assert_eq!(c.name, "LINKS-Nature");
        assert_eq!(c.base_code, Some(0x03C2));
        assert_eq!(c.scalars.get(&0x0000_0600), Some(&0x82));
        assert_eq!(c.sharpening(), Some(0.5));
        assert_eq!(c.clarity(), Some(0.5));
        assert_eq!(scalar_name(0x0000_0600), Some("sharpening"));
        assert_eq!(scalar_name(0x0000_0700), Some("clarity"));
        let cv = c.curve.as_ref().unwrap();
        assert_eq!(cv.point_count, 5);
        assert_eq!(cv.points.len(), 5);
        assert_eq!(cv.lut.len(), 257);
        assert_eq!(cv.lut_offset, 64);
        assert!(cv.is_monotonic());
        // 合成的是恒等曲线
        assert!(!c.has_custom_curve());
    }

    /// `0300` 代的曲线块载荷内 LUT 偏移比 `0310` 少 2 字节——写死偏移必错。
    #[test]
    fn parses_0300_template_with_shifted_lut() {
        let d = synth(b"0300", "Kodak Ektar Green", 0x00C3, 62);
        let c = parse(&d).unwrap();
        assert_eq!(c.family, "0300");
        assert_eq!(c.base_code, Some(0x00C3));
        let cv = c.curve.as_ref().unwrap();
        assert_eq!(cv.lut_offset, 62, "应回退扫描到 62");
        assert!(cv.is_monotonic());
    }

    /// 反例：若把 LUT 读在 64 而实际在 62，得到的数据必须被拒绝而不是静默接受。
    #[test]
    fn wrong_lut_offset_is_rejected_not_silently_accepted() {
        let d = synth(b"0300", "X", 0x0000, 62);
        let c = parse(&d).unwrap();
        let cv = c.curve.as_ref().unwrap();
        assert_ne!(cv.lut_offset, 64);
        // 直接按 64 读会越界或非单调，locate_lut 的接受条件会挡住
        let payload = &c.chunks.iter().find(|ch| ch.tag == tag::CURVE).unwrap().payload;
        let bad: Vec<u16> = (0..257)
            .map(|i| u16_be(&payload[64 + i * 2..64 + i * 2 + 2]))
            .collect();
        let looks_valid = bad.windows(2).all(|w| w[0] < w[1]) && bad[0] < 0x1000 && bad[256] > 0x7000;
        assert!(!looks_valid, "偏移错误时必须无法通过接受条件");
    }

    #[test]
    fn rejects_non_container() {
        let d = vec![0xAAu8; 64];
        assert!(matches!(parse(&d), Err(Error::NotRecipeContainer(_))));
    }

    #[test]
    fn retains_unknown_chunks() {
        let d = synth(b"0310", "X", 0x0000, 64);
        let c = parse(&d).unwrap();
        assert!(c.chunks.iter().any(|ch| ch.tag == tag::UNKNOWN_0400));
    }

    /// 中性值在文件里存为 `0x80`，解码后必须是 0——`as i8` 强转会得到 −128。
    #[test]
    fn neutral_bytes_decode_to_zero_not_minus_128() {
        let mut p = vec![0u8; 28];
        for b in p.iter_mut() {
            *b = 0x80;
        }
        let bands = parse_blender(&p).unwrap();
        for b in &bands {
            assert_eq!((b.hue, b.chroma, b.brightness), (0, 0, 0));
        }
        let mut g = vec![0u8; 20];
        for b in g.iter_mut() {
            *b = 0x80;
        }
        let gr = parse_grading(&g).unwrap();
        assert_eq!(gr.highlights.chroma, 0);
        assert_eq!(gr.shadows.brightness, 0);
    }

    /// 有符号解码的端点。
    #[test]
    fn signed_decode_endpoints() {
        assert_eq!(s8(0x80), 0);
        assert_eq!(s8(0x00), -128);
        assert_eq!(s8(0xFF), 127);
        assert_eq!(s8(0xE4), 100); // 0xE4 = 228 → 100
    }

    #[test]
    fn truncated_container_errors() {
        let d = &synth(b"0310", "X", 0, 64)[..20];
        assert!(parse(d).is_err());
    }

    /// 结束符是单个 u32=0。实测真实文件在末块之后恰好只剩 4 字节，
    /// 按"8 字节全零"处理会误报截断——这条用例把它锁死。
    #[test]
    fn terminator_is_a_single_u32_zero() {
        let full = synth(b"0310", "X", 0x0000, 64);
        // 合成器写的是 8 字节结束符；砍掉后 4 字节，模拟真实文件
        let four = &full[..full.len() - 4];
        let c = parse(four).expect("4 字节结束符必须被接受");
        assert_eq!(c.name, "X");
        assert!(c.curve.is_some());
    }

    /// 末块之后不足 4 字节必须报截断，而不是静默当成正常结束。
    #[test]
    fn short_tail_is_truncation() {
        let full = synth(b"0310", "X", 0x0000, 64);
        let short = &full[..full.len() - 5];
        assert!(matches!(parse(short), Err(Error::Truncated { .. })));
    }
}
