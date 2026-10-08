//! 从 NEF/JPEG 读出 Picture Control 的身份信息。
//!
//! 载荷位于 Nikon MakerNote 的 tag `0x0023`，按 ExifTool 的 `PictureControl3`
//! 布局解码。实测该载荷在 11926 个文件上**恒为 108 字节**。
//!
//! # 为什么不用字节模式匹配
//!
//! 载荷之前存在 8 字节前导头，其中间一个 u16 会随文件变化（实测 `0x01` 与
//! `0x03`）。以固定字节串定位会把该字段当成常量，从而把变化误判为"另一种
//! 格式"。本实现**只经 IFD 链读 tag**，那个字段从构造上就不会被触及。

use crate::error::{Error, Result};
use crate::makernote::nikon_maker_note;

/// 配方身份标签。
pub const TAG_PICTURE_CONTROL: u16 = 0x0023;

/// PictureControl3 各字段相对载荷起点的偏移。
mod off {
    pub const VERSION: usize = 0;
    pub const VERSION_LEN: usize = 8;
    pub const NAME: usize = 8;
    pub const NAME_LEN: usize = 20;
    pub const BASE: usize = 28;
    pub const BASE_LEN: usize = 20;
    pub const ADJUST: usize = 54;
    pub const QUICK_ADJUST: usize = 55;
    pub const SHARPNESS: usize = 57;
    pub const MID_RANGE_SHARPNESS: usize = 59;
    pub const CLARITY: usize = 61;
    pub const CONTRAST: usize = 63;
    pub const BRIGHTNESS: usize = 65;
    pub const SATURATION: usize = 67;
    pub const HUE: usize = 69;
    pub const FILTER_EFFECT: usize = 71;
    pub const TONING_EFFECT: usize = 72;
    pub const TONING_SATURATION: usize = 73;
    /// PictureControl3 具名字段的末尾（不含后续未定名字节）。
    pub const NAMED_END: usize = 74;
}

/// 调整模式（载荷 @54）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Adjust {
    DefaultSettings,
    QuickAdjust,
    FullControl,
    Other(u8),
}

impl Adjust {
    fn from_u8(v: u8) -> Self {
        match v {
            0 => Adjust::DefaultSettings,
            1 => Adjust::QuickAdjust,
            2 => Adjust::FullControl,
            other => Adjust::Other(other),
        }
    }
}

/// 曲线哨兵值：存在自定义色调曲线时，对比度与亮度会被覆写为该值。
pub const CURVE_SENTINEL: u8 = 0x01;

/// 一份配方身份。
///
/// 不含 `Eq`：锐化/清晰度等以四分之一步进表示，用 `f32` 承载更直观。
#[derive(Debug, Clone, PartialEq)]
pub struct Identity {
    /// 载荷前 8 字节（实测恒为 `03100310`，与配方无关，仅作记录）。
    pub version: String,
    /// 自定义配方名；内置配方时与 `base` 相同。
    pub name: String,
    /// 基准色彩名（如 `NEUTRAL`、`FLEXIBLE COLOR`）。
    pub base: String,
    pub adjust: Adjust,
    pub quick_adjust: u8,
    /// 锐化，按 `(b − 0x80) / 4` 解码。
    pub sharpness: f32,
    pub mid_range_sharpness: f32,
    pub clarity: f32,
    /// 对比度，按 `b − 0x80` 解码。
    pub contrast: i16,
    pub brightness: i16,
    pub saturation: i16,
    pub hue: i16,
    pub filter_effect: u8,
    pub toning_effect: u8,
    pub toning_saturation: u8,
    /// 载荷原始字节，用于交叉校验与诊断。
    pub raw: Vec<u8>,
}

impl Identity {
    /// 对比度与亮度是否同时为哨兵值。
    ///
    /// **这是提示，不是判据。** 实测存在反例：某配方确含 257 级非恒等曲线，
    /// 其哨兵却为假。因此调用方 MUST NOT 据此跳过配方库解析。
    pub fn curve_sentinel(&self) -> bool {
        self.raw.get(off::CONTRAST) == Some(&CURVE_SENTINEL)
            && self.raw.get(off::BRIGHTNESS) == Some(&CURVE_SENTINEL)
    }

    /// 是否为自定义配方（名称与基准不同）。
    pub fn is_custom(&self) -> bool {
        !self.name.eq_ignore_ascii_case(&self.base)
    }
}

fn fixed_ascii(b: &[u8]) -> String {
    let end = b.iter().position(|&c| c == 0).unwrap_or(b.len());
    String::from_utf8_lossy(&b[..end]).trim().to_string()
}

fn quarter(b: u8) -> f32 {
    (b as f32 - 128.0) / 4.0
}

fn linear(b: u8) -> i16 {
    b as i16 - 0x80
}

/// 从载荷字节解码身份。载荷长度不足具名区时返回截断错误。
pub fn decode(raw: &[u8]) -> Result<Identity> {
    if raw.len() < off::NAMED_END {
        return Err(Error::Truncated {
            what: "PictureControl3 载荷",
            need: off::NAMED_END,
            have: raw.len(),
        });
    }
    let at = |o: usize| raw[o];
    Ok(Identity {
        version: fixed_ascii(&raw[off::VERSION..off::VERSION + off::VERSION_LEN]),
        name: fixed_ascii(&raw[off::NAME..off::NAME + off::NAME_LEN]),
        base: fixed_ascii(&raw[off::BASE..off::BASE + off::BASE_LEN]),
        adjust: Adjust::from_u8(at(off::ADJUST)),
        quick_adjust: at(off::QUICK_ADJUST),
        sharpness: quarter(at(off::SHARPNESS)),
        mid_range_sharpness: quarter(at(off::MID_RANGE_SHARPNESS)),
        clarity: quarter(at(off::CLARITY)),
        contrast: linear(at(off::CONTRAST)),
        brightness: linear(at(off::BRIGHTNESS)),
        saturation: linear(at(off::SATURATION)),
        hue: linear(at(off::HUE)),
        filter_effect: at(off::FILTER_EFFECT),
        toning_effect: at(off::TONING_EFFECT),
        toning_saturation: at(off::TONING_SATURATION),
        raw: raw.to_vec(),
    })
}

/// 直接从文件字节读出配方身份。
pub fn read(data: &[u8]) -> Result<Identity> {
    let mn = nikon_maker_note(data)?;
    let ifd = mn.ifd0_offset()?;
    let entry = mn
        .find(ifd, TAG_PICTURE_CONTROL)?
        .ok_or(Error::NoPictureControlTag)?;
    let raw = mn.bytes(&entry)?;
    decode(raw)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 构造一份符合 PictureControl3 布局的合成载荷。
    fn synth(name: &str, base: &str, contrast: u8, brightness: u8) -> Vec<u8> {
        let mut b = vec![0u8; 108];
        b[..8].copy_from_slice(b"03100310");
        b[off::NAME..off::NAME + name.len()].copy_from_slice(name.as_bytes());
        b[off::BASE..off::BASE + base.len()].copy_from_slice(base.as_bytes());
        b[off::ADJUST] = 1;
        b[off::QUICK_ADJUST] = 0xFF;
        b[off::SHARPNESS] = 0x82; // +0.5
        b[off::MID_RANGE_SHARPNESS] = 0x80; // 0
        b[off::CLARITY] = 0x82; // +0.5
        b[off::CONTRAST] = contrast;
        b[off::BRIGHTNESS] = brightness;
        b[off::SATURATION] = 0x85; // +5
        b[off::HUE] = 0x78; // -8
        b[off::FILTER_EFFECT] = 0xFF;
        b[off::TONING_EFFECT] = 0xFF;
        b[off::TONING_SATURATION] = 0xFF;
        b
    }

    #[test]
    fn decodes_custom_recipe() {
        let raw = synth("LINKS-Nature", "NEUTRAL", 0x01, 0x01);
        let id = decode(&raw).unwrap();
        assert_eq!(id.name, "LINKS-Nature");
        assert_eq!(id.base, "NEUTRAL");
        assert!(id.is_custom());
        assert_eq!(id.adjust, Adjust::QuickAdjust);
        assert_eq!(id.sharpness, 0.5);
        assert_eq!(id.mid_range_sharpness, 0.0);
        assert_eq!(id.clarity, 0.5);
        assert_eq!(id.saturation, 5);
        assert_eq!(id.hue, -8);
        assert!(id.curve_sentinel());
    }

    #[test]
    fn decodes_builtin_recipe() {
        let raw = synth("PORTRAIT", "PORTRAIT", 0x80, 0x80);
        let id = decode(&raw).unwrap();
        assert!(!id.is_custom());
        assert!(!id.curve_sentinel());
        assert_eq!(id.contrast, 0);
        assert_eq!(id.brightness, 0);
    }

    #[test]
    fn sentinel_requires_both_fields() {
        let mut raw = synth("X", "NEUTRAL", 0x01, 0x01);
        raw[off::BRIGHTNESS] = 0x80;
        let id = decode(&raw).unwrap();
        assert!(!id.curve_sentinel(), "仅一个字段为哨兵时不得判为真");
    }

    #[test]
    fn truncated_payload_is_an_error() {
        let raw = vec![0u8; 40];
        assert!(matches!(decode(&raw), Err(Error::Truncated { .. })));
    }
}
