//! 解析失败的类型化错误。
//!
//! 解析的是逆向得来的二进制格式，因此"明确的失败"比"猜测后继续"重要：
//! 每个变体都对应一个可诊断的具体原因，调用方据此决定是跳过、报错还是降级。

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// 数据不足以完成解析（截断）。
    Truncated { what: &'static str, need: usize, have: usize },
    /// 未找到 TIFF 头。
    NotTiff,
    /// 未找到 Exif APP1 段（JPEG）。
    NoExifSegment,
    /// 未找到 MakerNote 或其头部不是 `Nikon\0`。
    NoMakerNote,
    /// MakerNote 中不含配方标签（tag 0x0023）。
    NoPictureControlTag,
    /// 字节序标记既不是 `II` 也不是 `MM`。
    BadByteOrder([u8; 2]),
    /// 容器 magic 不是 `NCP\0`。
    NotRecipeContainer([u8; 4]),
    /// 配方名在配方库中不存在。
    RecipeNotFound(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Truncated { what, need, have } => {
                write!(f, "{what} 数据不足：需要 {need} 字节，实际 {have} 字节")
            }
            Error::NotTiff => write!(f, "未找到 TIFF 头"),
            Error::NoExifSegment => write!(f, "JPEG 中未找到 Exif APP1 段"),
            Error::NoMakerNote => write!(f, "未找到 Nikon MakerNote"),
            Error::NoPictureControlTag => write!(f, "MakerNote 中不含配方标签 0x0023"),
            Error::BadByteOrder(b) => write!(f, "非法字节序标记：{b:02X?}"),
            Error::NotRecipeContainer(m) => write!(f, "非配方容器，magic = {m:02X?}"),
            Error::RecipeNotFound(n) => write!(f, "配方库中未找到配方：{n}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;
