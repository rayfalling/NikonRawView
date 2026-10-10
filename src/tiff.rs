//! TIFF / IFD 解析，NEF 与 JPEG 共用一套实现。
//!
//! 两者结构相同，唯一差别是 **TIFF 头的位置**：
//! - NEF：TIFF 头就在文件偏移 0
//! - JPEG：TIFF 头在 APP1 段内 `Exif\0\0` 之后
//!
//! 所有 IFD 条目中的偏移都相对于 TIFF 头，因此解析器把该位置作为 `base`
//! 携带。忽略这一点正是早前 JPEG 侧全部解析失败的根因——本模块的
//! `locate()` 与测试一起把这条规则锁死。

use crate::error::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteOrder {
    Little,
    Big,
}

impl ByteOrder {
    pub fn u16(self, b: &[u8]) -> u16 {
        let a = [b[0], b[1]];
        match self {
            ByteOrder::Little => u16::from_le_bytes(a),
            ByteOrder::Big => u16::from_be_bytes(a),
        }
    }
    pub fn u32(self, b: &[u8]) -> u32 {
        let a = [b[0], b[1], b[2], b[3]];
        match self {
            ByteOrder::Little => u32::from_le_bytes(a),
            ByteOrder::Big => u32::from_be_bytes(a),
        }
    }
}

/// 一个 12 字节的 IFD 条目。
#[derive(Debug, Clone, Copy)]
pub struct Entry {
    pub tag: u16,
    /// TIFF 类型码（3 = SHORT，4 = LONG，7 = UNDEFINED …）。
    pub typ: u16,
    pub count: u32,
    /// 条目自身的绝对偏移（12 字节块起始）。
    pub offset: usize,
}

/// 每个 TIFF 类型码的单元字节数；未知类型返回 0。
fn unit_size(typ: u16) -> usize {
    match typ {
        1 | 2 | 6 | 7 => 1, // BYTE / ASCII / SBYTE / UNDEFINED
        3 | 8 => 2,         // SHORT / SSHORT
        4 | 9 | 11 => 4,    // LONG / SLONG / FLOAT
        5 | 10 | 12 => 8,   // RATIONAL / SRATIONAL / DOUBLE
        _ => 0,
    }
}

/// 解析后的 TIFF 视图：持有原始数据与 TIFF 头基址。
pub struct Tiff<'a> {
    data: &'a [u8],
    /// TIFF 头在 `data` 中的绝对位置。NEF 为 0，JPEG 为 APP1 内 Exif 之后。
    base: usize,
    pub order: ByteOrder,
}

impl<'a> Tiff<'a> {
    /// 定位 TIFF 头：先看文件开头，再找 JPEG 的 Exif APP1 段。
    pub fn locate(data: &'a [u8]) -> Result<Self> {
        if data.len() >= 8 && (data.starts_with(b"II") || data.starts_with(b"MM")) {
            return Self::at(data, 0);
        }
        // JPEG：在 `Exif\0\0` 之后
        let needle = b"Exif\x00\x00";
        let pos = data
            .windows(needle.len())
            .position(|w| w == needle)
            .ok_or(Error::NoExifSegment)?;
        Self::at(data, pos + needle.len())
    }

    /// 在指定基址构造；校验字节序标记与首个 IFD 偏移。
    pub fn at(data: &'a [u8], base: usize) -> Result<Self> {
        let head = data
            .get(base..base + 8)
            .ok_or(Error::Truncated { what: "TIFF 头", need: 8, have: data.len().saturating_sub(base) })?;
        let order = match &head[..2] {
            b"II" => ByteOrder::Little,
            b"MM" => ByteOrder::Big,
            other => return Err(Error::BadByteOrder([other[0], other[1]])),
        };
        let magic = order.u16(&head[2..4]);
        if magic != 42 {
            return Err(Error::NotTiff);
        }
        Ok(Self { data, base, order })
    }

    pub fn base(&self) -> usize {
        self.base
    }

    /// 底层字节。偏移量与 [`Entry`] 里给出的一样，是**绝对**偏移。
    ///
    /// 供需要自行批量读取的调用方使用（如按条带读 16 位像素）——那条路径上逐条目
    /// 走 [`Self::bytes`] 会把每 2 字节都过一次边界检查。
    pub fn raw(&self) -> &'a [u8] {
        self.data
    }

    /// IFD0 的绝对偏移。
    pub fn ifd0_offset(&self) -> Result<usize> {
        let p = self.base + 4;
        let raw = self
            .data
            .get(p..p + 4)
            .ok_or(Error::Truncated { what: "IFD0 指针", need: 4, have: 0 })?;
        Ok(self.base + self.order.u32(raw) as usize)
    }

    /// 读出一个 IFD 的全部条目。
    pub fn entries(&self, ifd_offset: usize) -> Result<Vec<Entry>> {
        let cnt_raw = self.data.get(ifd_offset..ifd_offset + 2).ok_or(Error::Truncated {
            what: "IFD 条目数",
            need: 2,
            have: self.data.len().saturating_sub(ifd_offset),
        })?;
        let n = self.order.u16(cnt_raw) as usize;
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let p = ifd_offset + 2 + i * 12;
            let blk = self.data.get(p..p + 12).ok_or(Error::Truncated {
                what: "IFD 条目",
                need: 12,
                have: self.data.len().saturating_sub(p),
            })?;
            out.push(Entry {
                tag: self.order.u16(&blk[0..2]),
                typ: self.order.u16(&blk[2..4]),
                count: self.order.u32(&blk[4..8]),
                offset: p,
            });
        }
        Ok(out)
    }

    /// 条目数据的位置与长度。数据量 ≤ 4 字节时内联在条目内，否则是一个相对基址的偏移。
    fn data_range(&self, e: &Entry) -> Result<(usize, usize)> {
        let unit = unit_size(e.typ);
        let len = unit.saturating_mul(e.count as usize);
        if unit == 0 {
            return Ok((e.offset + 8, 0));
        }
        if len <= 4 {
            return Ok((e.offset + 8, len));
        }
        let p = e.offset + 8;
        let raw = self
            .data
            .get(p..p + 4)
            .ok_or(Error::Truncated { what: "条目偏移", need: 4, have: 0 })?;
        let abs = self.base + self.order.u32(raw) as usize;
        Ok((abs, len))
    }

    /// 条目原始字节。
    pub fn bytes(&self, e: &Entry) -> Result<&'a [u8]> {
        let (p, len) = self.data_range(e)?;
        if len == 0 {
            return Ok(&[]);
        }
        self.data.get(p..p + len).ok_or(Error::Truncated {
            what: "条目数据",
            need: len,
            have: self.data.len().saturating_sub(p),
        })
    }

    /// 单值 SHORT。
    pub fn u16_value(&self, e: &Entry) -> Result<u16> {
        let b = self.bytes(e)?;
        if b.len() < 2 {
            return Err(Error::Truncated { what: "SHORT", need: 2, have: b.len() });
        }
        Ok(self.order.u16(&b[..2]))
    }

    /// 单值 LONG（亦用于承接指针）。
    pub fn u32_value(&self, e: &Entry) -> Result<u32> {
        let b = self.bytes(e)?;
        if b.len() < 4 {
            return Err(Error::Truncated { what: "LONG", need: 4, have: b.len() });
        }
        Ok(self.order.u32(&b[..4]))
    }

    /// 按 tag 查找条目（返回首个匹配）。
    pub fn find(&self, ifd_offset: usize, tag: u16) -> Result<Option<Entry>> {
        Ok(self.entries(ifd_offset)?.into_iter().find(|e| e.tag == tag))
    }

    /// 条目**值字段**（12 字节条目中的最后 4 字节）本身，不做任何解引用。
    ///
    /// 指针类标签（如 ExifIFD、MakerNote）必须用它——这些标签的 count 往往大于 4，
    /// 若走 [`Self::bytes`] 会"顺着偏移把被指向的数据取回来"，再把那批数据的头 4
    /// 字节当成偏移，得到完全错误的地址。
    pub fn value_field(&self, e: &Entry) -> Result<&'a [u8]> {
        let p = e.offset + 8;
        self.data
            .get(p..p + 4)
            .ok_or(Error::Truncated { what: "条目值字段", need: 4, have: 0 })
    }

    /// 值字段按 LONG 解读（用于指针）。
    pub fn value_field_u32(&self, e: &Entry) -> Result<u32> {
        Ok(self.order.u32(self.value_field(e)?))
    }

    /// 某个条目所指向的**相对**偏移换算成绝对位置。
    pub fn deref_offset(&self, e: &Entry) -> Result<usize> {
        Ok(self.base + self.value_field_u32(e)? as usize)
    }

    /// MakerNote 的绝对偏移（tag 0x927C，位于 ExifIFD）。
    pub fn maker_note_offset(&self) -> Result<usize> {
        let ifd0 = self.ifd0_offset()?;
        let exif = self
            .find(ifd0, 0x8769)?
            .ok_or(Error::NoMakerNote)?;
        let exif_abs = self.deref_offset(&exif)?;
        let mn = self
            .find(exif_abs, 0x927C)?
            .ok_or(Error::NoMakerNote)?;
        let abs = self.deref_offset(&mn)?;
        // Nikon MakerNote 自带 "Nikon\0" 头
        match self.data.get(abs..abs + 6) {
            Some(h) if h == b"Nikon\x00" => Ok(abs),
            _ => Err(Error::NoMakerNote),
        }
    }
}
