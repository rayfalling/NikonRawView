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
