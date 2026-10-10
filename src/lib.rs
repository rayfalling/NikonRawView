//! 解析尼康 Picture Control 配方。
//!
//! 本 crate 覆盖两件事——**从照片读配方身份**、**从配方文件读配方本体**：
//!
//! - [`picture_control`]：NEF/JPEG 的 MakerNote 中 tag `0x0023` 的解码
//! - [`np3`]：NP3/NCP 配方容器的解析（两代模板）
//! - [`library`]：本机配方来源的发现与按名解析
//! - [`tiff`]：TIFF/IFD 解析，NEF 与 JPEG 共用（以基址参数区分）
//!
//! # 为什么是两段式
//!
//! 实测确认 **NEF 装不下完整配方**：相机写入的载荷恒为 108 字节，而一条
//! 257 级的色调曲线就需要 514 字节。因此身份在照片里，本体在配方文件里。

pub mod calibrate;
pub mod camera;
pub mod color;
pub mod deltae;
pub mod error;
pub mod fit;
pub mod icc;
pub mod library;
pub mod libraw;
pub mod makernote;
pub mod mat3;
pub mod np3;
pub mod picture_control;
pub mod render;
pub mod tiff;
pub mod transform;

pub use error::{Error, Result};
pub use picture_control::{Identity, Adjust};
