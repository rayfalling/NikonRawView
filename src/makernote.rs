//! Nikon MakerNote 定位。
//!
//! Nikon 的 MakerNote 自带一个内嵌 TIFF 头（`Nikon\0` + 2 字节版本 + 2 字节填充，
//! 其后是 `II`/`MM` + 42 + IFD 偏移）。MakerNote 内部所有条目的偏移都相对那个
//! 内嵌 TIFF 头，因此这里复用 [`Tiff::at`]，把基址设为 `makernote + 10`。

use crate::error::Result;
use crate::tiff::Tiff;

/// `Nikon\0` 头长度 + 版本 2 字节 + 填充 2 字节。
pub const NIKON_MN_HEADER: usize = 10;

/// 返回指向 MakerNote IFD 的 TIFF 视图（基址为内嵌 TIFF 头）。
pub fn nikon_maker_note(data: &[u8]) -> Result<Tiff<'_>> {
    let outer = Tiff::locate(data)?;
    let mn = outer.maker_note_offset()?;
    Tiff::at(data, mn + NIKON_MN_HEADER)
}

/// Nikon MakerNote 的 Active D-Lighting 标签。
pub const TAG_ACTIVE_D_LIGHTING: u16 = 0x0022;

/// Active D-Lighting 的档位。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveDLighting {
    Off,
    Low,
    Normal,
    High,
    ExtraHigh,
    /// 自动。
    Auto,
    /// 未记录的取值。
    Other(u16),
}

impl ActiveDLighting {
    pub fn from_raw(v: u16) -> Self {
        match v {
            0 => Self::Off,
            1 => Self::Low,
            3 => Self::Normal,
            5 => Self::High,
            7 => Self::ExtraHigh,
            0xFFFF => Self::Auto,
            other => Self::Other(other),
        }
    }

    pub fn is_off(&self) -> bool {
        matches!(self, Self::Off)
    }

    pub fn name(&self) -> String {
        match self {
            Self::Off => "关闭".into(),
            Self::Low => "低".into(),
            Self::Normal => "标准".into(),
            Self::High => "高".into(),
            Self::ExtraHigh => "超高".into(),
            Self::Auto => "自动".into(),
            Self::Other(v) => format!("未知({v})"),
        }
    }
}

/// 读出一个 MakerNote 标签的原始字节。标签不存在时返回 `None`。
pub fn tag_bytes(data: &[u8], tag: u16) -> Result<Option<Vec<u8>>> {
    let mn = nikon_maker_note(data)?;
    let ifd = mn.ifd0_offset()?;
    match mn.find(ifd, tag)? {
        Some(e) => Ok(Some(mn.bytes(&e)?.to_vec())),
        None => Ok(None),
    }
}

/// 读出一个单值 SHORT 标签。
pub fn tag_u16(data: &[u8], tag: u16) -> Result<Option<u16>> {
    let mn = nikon_maker_note(data)?;
    let ifd = mn.ifd0_offset()?;
    match mn.find(ifd, tag)? {
        Some(e) => Ok(Some(mn.u16_value(&e)?)),
        None => Ok(None),
    }
}

/// 读出 Active D-Lighting 档位。
pub fn active_d_lighting(data: &[u8]) -> Result<Option<ActiveDLighting>> {
    Ok(tag_u16(data, TAG_ACTIVE_D_LIGHTING)?.map(ActiveDLighting::from_raw))
}

/// Nikon MakerNote 的 **BlackLevel** 标签。
///
/// 这是 Z8 黑电平的**唯一来源**（LibRaw `src/metadata/nikon.cpp:578-585`）。
/// `adobe_coeff` 表项里 Z8 的 `t_black` 是 0，不参与；12-bit 缩放分支也不触发
/// （Z8 是 14 位）。实测四个样本都是 `[1008, 1008, 1008, 1008]`。
pub const TAG_BLACK_LEVEL: u16 = 0x003d;

/// 黑电平的四通道值，以及按 LibRaw 语义归算出的公共部分。
///
/// # 为什么要复刻 LibRaw 的归算
///
/// LibRaw 读 tag `0x003d` 后并不直接当四通道用，而是先取**四者最小值**作为公共黑
/// 电平 `black`，逐通道相对它做差（`nikon.cpp:580-584`）：
///
/// ```text
/// black = min(cblack[0..3])
/// cblack[c] -= black        // 剩下的部分才是逐通道差异
/// ```
///
/// 而 `RGGB_2_RGBG(q) = q ^ (q >> 1)`（`internal/defines.h:171`）把文件序
/// `[R, G1, B, G2]` 映射到内部序 `[R, G1, G2, B]`。实测 Z8 四个值相同，所以
/// `black = 1008`、逐通道差全为 0——**等效于统一扣 1008**。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlackLevel {
    /// 文件序的四通道原始值（RGGB 顺序）。
    pub raw: [u16; 4],
    /// 公共黑电平（四者最小值）——LibRaw 的 `black`。
    pub black: u16,
    /// 逐通道相对公共部分之差——LibRaw 的 `cblack[0..3]`（内部序 R, G1, G2, B）。
    pub per_channel: [u16; 4],
}

impl BlackLevel {
    /// 四通道是否完全一致（此时逐通道差全为 0，等效于统一扣一个标量）。
    pub fn is_uniform(&self) -> bool {
        self.per_channel.iter().all(|v| *v == 0)
    }

    /// 可读描述，供诊断与日志。
    pub fn describe(&self) -> String {
        if self.is_uniform() {
            format!("黑电平 {}（四通道一致：{:?}）", self.black, self.raw)
        } else {
            format!(
                "黑电平 公共 {} + 逐通道 {:?}（文件序 {:?}）",
                self.black, self.per_channel, self.raw
            )
        }
    }
}

/// 读出 Nikon MakerNote 的黑电平（tag `0x003d`）。
///
/// 返回 `None` 表示该文件没有这个标签——**不要当成 0**，那会让暗部留下未知的
/// pedestal。调用方应把"读不到"与"读到 0"区分开。
pub fn black_level(data: &[u8]) -> Result<Option<BlackLevel>> {
    let mn = nikon_maker_note(data)?;
    let ifd = mn.ifd0_offset()?;
    let Some(e) = mn.find(ifd, TAG_BLACK_LEVEL)? else {
        return Ok(None);
    };
    // 必须是 4 个 SHORT —— 类型或个数不符就不要猜
    if e.typ != 3 || e.count != 4 {
        return Ok(None);
    }
    let b = mn.bytes(&e)?;
    if b.len() < 8 {
        return Ok(None);
    }
    let order = mn.order;
    let mut raw = [0u16; 4];
    for (i, slot) in raw.iter_mut().enumerate() {
        *slot = order.u16(&b[i * 2..i * 2 + 2]);
    }

    // 复刻 LibRaw 语义：公共部分取最小值，逐通道减掉它
    let black = raw.iter().copied().min().unwrap_or(0);
    // RGGB_2_RGBG(q) = q ^ (q >> 1)：文件序 [R,G1,B,G2] → 内部序 [R,G1,G2,B]
    let mut per_channel = [0u16; 4];
    for (c, v) in raw.iter().enumerate() {
        per_channel[c ^ (c >> 1)] = v.saturating_sub(black);
    }
    Ok(Some(BlackLevel { raw, black, per_channel }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rggb_to_rgbg_mapping_is_the_libraw_one() {
        // RGGB_2_RGBG(q) = q ^ (q >> 1)
        let m: Vec<usize> = (0..4).map(|q| q ^ (q >> 1)).collect();
        assert_eq!(m, vec![0, 1, 3, 2], "文件序 [R,G1,B,G2] 应映射到内部序 [R,G1,G2,B]");
    }

    #[test]
    fn uniform_black_level_is_recognised() {
        let b = BlackLevel {
            raw: [1008, 1008, 1008, 1008],
            black: 1008,
            per_channel: [0, 0, 0, 0],
        };
        assert!(b.is_uniform());
        assert!(b.describe().contains("1008"));
    }

    #[test]
    fn non_uniform_black_level_reports_per_channel() {
        // 造一个逐通道不同的例子，验证归算与描述
        let raw = [1008u16, 1010, 1012, 1010];
        let black = *raw.iter().min().unwrap();
        assert_eq!(black, 1008);
        let mut per = [0u16; 4];
        for c in 0..4usize {
            per[c ^ (c >> 1)] = raw[c] - black;
        }
        let b = BlackLevel { raw, black, per_channel: per };
        assert!(!b.is_uniform());
        let d = b.describe();
        assert!(d.contains("公共 1008"), "{d}");
        assert!(d.contains("逐通道"), "{d}");
    }
}
