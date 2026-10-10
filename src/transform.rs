//! 基准变换数据的文件格式。
//!
//! # 为什么要单独一份数据文件
//!
//! 基准配方的渲染变换不在任何文件里（已实测：NX Studio 的 `PicCon*.bin` 是携带配方
//! 记录的模板 TIFF，**不含 257 级曲线**），只能从参考导出拟合。拟合结果必须能脱离
//! 参考样张独立分发——**日常渲染不应依赖参考导出文件的存在**。
//!
//! # 为什么是「一维色调曲线 + 三维色彩查找表」
//!
//! 基准渲染同时包含**非线性色调**与**色彩重映射**两部分：纯矩阵表达不了色调曲线，
//! 纯 3D LUT 覆盖不了无限动态范围。拆成两部分后各自可拟合、可单独诊断。
//!
//! # 格式
//!
//! 小端、定长字段 + 变长载荷，全部按 4 字节对齐读取：
//!
//! ```text
//! 0x00  "NBT\0"          magic
//! 0x04  u32              版本
//! 0x08  u16              基准编码
//! 0x0A  u16              保留（必须为 0）
//! 0x0C  u32 + bytes      基准名（UTF-8）
//!       u32 + bytes      来源清单（UTF-8，换行分隔）
//!       u32 + f32×n      色调曲线
//!       u32              3D LUT 边长
//!       u32 + f32×n×3    3D LUT（行主序，R 变化最快）
//! ```
//!
//! 版本不匹配时**明确报错**，不做兼容猜测——猜错的代价是一整批照片被悄悄染错色。

use crate::error::{Error, Result};
use crate::render::Transform;

pub const MAGIC: &[u8; 4] = b"NBT\0";
/// 当前格式版本。任何不兼容改动都必须递增。
pub const VERSION: u32 = 1;

/// 一份基准变换。
#[derive(Debug, Clone, PartialEq)]
pub struct BaseTransform {
    /// 基准编码（来自配方身份的 base code）。**选择依据是它，不是名称。**
    pub base_code: u16,
    /// 基准名，仅用于展示与诊断。
    pub name: String,
    /// 拟合所用的参考文件清单。
    pub sources: Vec<String>,
    /// 一维色调曲线；空表示恒等。
    pub curve: Vec<f32>,
    /// 3D LUT 的边长（每轴采样数）。
    pub lut_edge: usize,
    /// 3D LUT，长度 = edge³ × 3，行主序，R 变化最快。
    pub lut: Vec<f32>,
}

impl BaseTransform {
    pub fn new(base_code: u16, name: impl Into<String>) -> Self {
        Self {
            base_code,
            name: name.into(),
            sources: Vec::new(),
            curve: Vec::new(),
            lut_edge: 0,
            lut: Vec::new(),
        }
    }

    /// 由拟合结果构造（见 [`crate::fit`]）。
    pub fn from_fit(
        base_code: u16,
        name: impl Into<String>,
        sources: Vec<String>,
        curve: Vec<f32>,
        lut_edge: usize,
        lut: Vec<f32>,
    ) -> Self {
        Self { base_code, name: name.into(), sources, curve, lut_edge, lut }
    }

    /// 施加到**工作空间的线性值**上：先曲线，后 LUT。
    ///
    /// 顺序不是随便定的——曲线负责明暗、LUT 是在曲线之后拟合出来的，它接收的正是
    /// 曲线之后的 RGB。反过来施加会让 LUT 落在它没被拟合过的定义域上。
    pub fn apply(&self, v: [f32; 3]) -> [f32; 3] {
        let c = crate::fit::forward_curve(&self.curve, v);
        if self.lut_edge < 2 || self.lut.is_empty() {
            c
        } else {
            crate::fit::apply_lut(&self.lut, self.lut_edge, c)
        }
    }

    /// 是否为恒等变换（没有曲线也没有 LUT）。
    pub fn is_identity(&self) -> bool {
        self.curve.is_empty() && self.lut.is_empty()
    }

    /// 转成可施加的变换。
    pub fn to_transform(&self) -> Transform {
        let mut t = Transform::default();
        if !self.curve.is_empty() {
            t.curve = Some(self.curve.clone());
        }
        t
    }

    /// 序列化。
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(64 + self.curve.len() * 4 + self.lut.len() * 4);
        out.extend_from_slice(MAGIC);
        put_u32(&mut out, VERSION);
        put_u16(&mut out, self.base_code);
        put_u16(&mut out, 0); // 保留位，固定为 0
        put_str(&mut out, &self.name);
        put_str(&mut out, &self.sources.join("\n"));
        put_u32(&mut out, self.curve.len() as u32);
        for v in &self.curve {
            put_f32(&mut out, *v);
        }
        put_u32(&mut out, self.lut_edge as u32);
        put_u32(&mut out, (self.lut.len() / 3) as u32);
        for v in &self.lut {
            put_f32(&mut out, *v);
        }
        out
    }

    /// 反序列化。版本不符或内容截断都会给出明确错误。
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut r = Reader { d: data, p: 0 };
        let magic = r.take(4, "magic")?;
        if magic != MAGIC {
            return Err(Error::BadPayload {
                what: "基准变换数据",
                detail: format!("magic 应为 {MAGIC:02X?}，实际为 {magic:02X?}"),
            });
        }
        let version = r.u32("版本")?;
        if version != VERSION {
            return Err(Error::UnsupportedVersion {
                what: "基准变换数据",
                found: version,
                expected: VERSION,
            });
        }
        let base_code = r.u16("基准编码")?;
        let reserved = r.u16("保留位")?;
        if reserved != 0 {
            return Err(Error::BadPayload {
                what: "基准变换数据",
                detail: format!("保留位应为 0，实际为 {reserved}"),
            });
        }
        let name = r.string("基准名")?;
        let sources_raw = r.string("来源清单")?;
        let curve_len = r.u32("曲线长度")? as usize;
        let mut curve = Vec::with_capacity(curve_len);
        for _ in 0..curve_len {
            curve.push(r.f32("曲线采样")?);
        }
        let lut_edge = r.u32("LUT 边长")? as usize;
        let lut_n = r.u32("LUT 项数")? as usize;
        if lut_edge > 0 && lut_n != lut_edge * lut_edge * lut_edge {
            return Err(Error::BadPayload {
                what: "基准变换数据",
                detail: format!("LUT 项数 {lut_n} 与边长 {lut_edge} 不符（应为 {}）", lut_edge.pow(3)),
            });
        }
        let mut lut = Vec::with_capacity(lut_n * 3);
        for _ in 0..lut_n * 3 {
            lut.push(r.f32("LUT 分量")?);
        }
        if r.p != data.len() {
            return Err(Error::BadPayload {
                what: "基准变换数据",
                detail: format!("末尾有 {} 字节多余数据", data.len() - r.p),
            });
        }

        Ok(Self {
            base_code,
            name,
            sources: sources_raw
                .split('\n')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(String::from)
                .collect(),
            curve,
            lut_edge,
            lut,
        })
    }
}

fn put_u16(o: &mut Vec<u8>, v: u16) {
    o.extend_from_slice(&v.to_le_bytes());
}
fn put_u32(o: &mut Vec<u8>, v: u32) {
    o.extend_from_slice(&v.to_le_bytes());
}
fn put_f32(o: &mut Vec<u8>, v: f32) {
    o.extend_from_slice(&v.to_le_bytes());
}
fn put_str(o: &mut Vec<u8>, s: &str) {
    put_u32(o, s.len() as u32);
    o.extend_from_slice(s.as_bytes());
}

struct Reader<'a> {
    d: &'a [u8],
    p: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize, what: &'static str) -> Result<&'a [u8]> {
        let s = self.d.get(self.p..self.p + n).ok_or_else(|| Error::BadPayload {
            what: "基准变换数据",
            detail: format!("读取{what}时数据不足：需要 {n} 字节，剩余 {}", self.d.len().saturating_sub(self.p)),
        })?;
        self.p += n;
        Ok(s)
    }
    fn u16(&mut self, what: &'static str) -> Result<u16> {
        let b = self.take(2, what)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }
    fn u32(&mut self, what: &'static str) -> Result<u32> {
        let b = self.take(4, what)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn f32(&mut self, what: &'static str) -> Result<f32> {
        let b = self.take(4, what)?;
        Ok(f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn string(&mut self, what: &'static str) -> Result<String> {
        let n = self.u32(what)? as usize;
        let b = self.take(n, what)?;
        String::from_utf8(b.to_vec()).map_err(|e| Error::BadPayload {
            what: "基准变换数据",
            detail: format!("{what} 不是合法 UTF-8：{e}"),
        })
    }
}

/// 把基准变换写入文件。
pub fn write(path: &std::path::Path, t: &BaseTransform) -> Result<()> {
    std::fs::write(path, t.to_bytes())
        .map_err(|e| Error::InvalidInput(format!("写入 {} 失败：{e}", path.display())))
}

/// 从文件读出基准变换。
pub fn read(path: &std::path::Path) -> Result<BaseTransform> {
    let d = std::fs::read(path)
        .map_err(|e| Error::InvalidInput(format!("读取 {} 失败：{e}", path.display())))?;
    BaseTransform::from_bytes(&d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> BaseTransform {
        let mut t = BaseTransform::new(0x03C2, "NEUTRAL");
        t.sources = vec!["simple/DSC_4143.NEF".into(), "Z:/摄影/Nikon/Z8/2025.12.20-12.27 冰岛/DSC_1001.NEF".into()];
        t.curve = (0..257).map(|i| (i as f32 / 256.0).powf(0.9)).collect();
        t.lut_edge = 3;
        t.lut = (0..27).flat_map(|i| [i as f32 / 27.0; 3]).collect();
        t
    }

    #[test]
    fn roundtrip_is_field_identical() {
        let t = sample();
        let back = BaseTransform::from_bytes(&t.to_bytes()).expect("应能读回");
        assert_eq!(back, t, "写入后读回应逐字段一致");
    }

    #[test]
    fn roundtrip_handles_empty_payloads() {
        let t = BaseTransform::new(0, "");
        let back = BaseTransform::from_bytes(&t.to_bytes()).unwrap();
        assert_eq!(back, t);
        assert!(back.is_identity());
    }

    #[test]
    fn roundtrip_survives_non_ascii_paths() {
        // 本项目的照片库路径普遍含中文，来源清单必须原样往返
        let t = sample();
        let back = BaseTransform::from_bytes(&t.to_bytes()).unwrap();
        assert_eq!(back.sources, t.sources);
        assert!(back.sources[1].contains("冰岛"));
    }

    #[test]
    fn version_mismatch_is_explicit() {
        let mut b = sample().to_bytes();
        b[4..8].copy_from_slice(&999u32.to_le_bytes());
        match BaseTransform::from_bytes(&b) {
            Err(Error::UnsupportedVersion { found, expected, .. }) => {
                assert_eq!((found, expected), (999, VERSION));
            }
            other => panic!("应报版本不受支持，得到 {other:?}"),
        }
    }

    #[test]
    fn bad_magic_is_rejected() {
        let mut b = sample().to_bytes();
        b[0] = b'X';
        assert!(matches!(BaseTransform::from_bytes(&b), Err(Error::BadPayload { .. })));
    }

    #[test]
    fn truncated_payload_is_rejected() {
        let full = sample().to_bytes();
        for cut in [4, 8, 12, 20, full.len() / 2, full.len() - 1] {
            let e = BaseTransform::from_bytes(&full[..cut]);
            assert!(e.is_err(), "截断到 {cut} 字节时应报错");
        }
    }

    #[test]
    fn trailing_garbage_is_rejected() {
        let mut b = sample().to_bytes();
        b.push(0);
        assert!(matches!(BaseTransform::from_bytes(&b), Err(Error::BadPayload { .. })));
    }

    #[test]
    fn nonzero_reserved_is_rejected() {
        let mut b = sample().to_bytes();
        b[10..12].copy_from_slice(&7u16.to_le_bytes());
        assert!(matches!(BaseTransform::from_bytes(&b), Err(Error::BadPayload { .. })));
    }

    #[test]
    fn inconsistent_lut_size_is_rejected() {
        let mut b = sample().to_bytes();
        // 把 LUT 项数改小，使 edge³ 与之不符
        let p = b.len() - (27 * 3 * 4) - 4;
        b[p..p + 4].copy_from_slice(&5u32.to_le_bytes());
        assert!(matches!(BaseTransform::from_bytes(&b), Err(Error::BadPayload { .. })));
    }

    #[test]
    fn file_roundtrip() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(format!("base-transform-test-{}.nbt", std::process::id()));
        let t = sample();
        write(&p, &t).expect("应能写入");
        let back = read(&p).expect("应能读出");
        let _ = std::fs::remove_file(&p);
        assert_eq!(back, t);
    }

    #[test]
    fn missing_file_is_an_explicit_error() {
        let p = std::path::Path::new("no-such-transform.nbt");
        assert!(matches!(read(p), Err(Error::InvalidInput(_))));
    }
}
