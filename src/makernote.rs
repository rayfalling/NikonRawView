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
