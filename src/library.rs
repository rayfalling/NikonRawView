//! 本机配方来源的发现与按名解析。
//!
//! 两个来源：
//! - 注册表 `HKCU\Software\Nikon\Common\PictureControl\NP3\CustomCurves`
//! - NX Studio 边车 `.nksc`（明文 XMP，`nikon::PictureControl` 滤镜内以 Base64 内嵌配方）
//!
//! 配方名一律取自**载荷内部**，而非注册表值名或文件名——实测注册表值名是数字
//! ID（如 `50331649`），文件内才是可读名。

use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::np3::{self, Container};

/// 一份配方的来源。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Source {
    /// Windows 注册表，记其值名。
    Registry { value_name: String },
    /// NX Studio 边车文件，记其内的字段名。
    Sidecar { path: PathBuf, field: String },
    /// 独立配方文件。
    File { path: PathBuf },
}

impl Source {
    pub fn label(&self) -> String {
        match self {
            Source::Registry { value_name } => format!("registry:{value_name}"),
            Source::Sidecar { path, field } => {
                format!("sidecar:{}#{field}", path.file_name().unwrap_or_default().to_string_lossy())
            }
            Source::File { path } => format!("file:{}", path.to_string_lossy()),
        }
    }
}

/// 一份已解析的配方。
#[derive(Debug, Clone)]
pub struct Recipe {
    pub name: String,
    pub source: Source,
    pub container: Container,
}

impl Recipe {
    /// 家族码，代次标识。
    pub fn family(&self) -> &str {
        &self.container.family
    }
}

/// 配方库：去重后的若干配方。
#[derive(Debug, Default)]
pub struct Library {
    pub recipes: Vec<Recipe>,
}

impl Library {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, r: Recipe) {
        self.recipes.push(r);
    }

    /// 按名解析（忽略大小写与首尾空白，因名称在文件内是定长 NUL 填充字段）。
    pub fn find(&self, name: &str) -> Vec<&Recipe> {
        let want = name.trim();
        self.recipes
            .iter()
            .filter(|r| r.name.eq_ignore_ascii_case(want))
            .collect()
    }

    /// 去重后的配方名集合。
    pub fn names(&self) -> Vec<String> {
        let mut v: Vec<String> = Vec::new();
        for r in &self.recipes {
            if !v.iter().any(|x| x.eq_ignore_ascii_case(&r.name)) {
                v.push(r.name.clone());
            }
        }
        v.sort();
        v
    }

    /// 按名解析，未命中时返回错误并附可用名清单。
    pub fn require(&self, name: &str) -> Result<&Recipe> {
        match self.find(name).first() {
            Some(r) => Ok(r),
            None if self.recipes.is_empty() => Err(Error::RecipeNotFound(name.to_string())),
            None => {
                let avail = self.names().join("、");
                Err(Error::RecipeNotFound(format!("{name}（本机可用：{avail}）")))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Base64（自实现，避免引入依赖）
// ---------------------------------------------------------------------------

/// 标准 Base64 解码。
///
/// 只容忍空白，遇到其他非字母表字符即停止——若放宽为"跳过一切未知字符"，
/// 未反转义的实体（如 `&#xA;`）中的 `A`、`x` 恰好是合法 Base64 字符，
/// 会被静默当成数据注入。因此调用方 MUST 先反转义 XML。
pub fn b64_decode(s: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len() / 4 * 3);
    let mut buf: u32 = 0;
    let mut bits: u32 = 0;
    for ch in s.bytes() {
        let v = match ch {
            b'A'..=b'Z' => ch - b'A',
            b'a'..=b'z' => ch - b'a' + 26,
            b'0'..=b'9' => ch - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' => break,
            b' ' | b'\t' | b'\r' | b'\n' => continue,
            _ => break,
        } as u32;
        buf = (buf << 6) | v;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((buf >> bits) as u8);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// .nksc 边车
// ---------------------------------------------------------------------------

/// 反转义 XML 实体并去掉空白，便于随后提取标签内容。
fn unescape(s: &str) -> String {
    let mut t = s
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&");
    // 去掉数字实体（如 &#xA;）——它们出现在 Base64 内部作换行
    let mut out = String::with_capacity(t.len());
    let b = t.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'&' && i + 2 < b.len() && b[i + 1] == b'#' {
            if let Some(semi) = t[i..].find(';') {
                i += semi + 1;
                continue;
            }
        }
        out.push(b[i] as char);
        i += 1;
    }
    t = out;
    t
}

/// 从 XMP 文本中取出某个标签的全部内容。
fn tag_contents<'a>(xml: &'a str, tag: &str) -> Vec<&'a str> {
    let open = format!("<{tag}>");
    let close = format!("</{tag}>");
    let mut out = Vec::new();
    let mut from = 0usize;
    while let Some(a) = xml[from..].find(&open) {
        let s = from + a + open.len();
        match xml[s..].find(&close) {
            Some(b) => {
                out.push(&xml[s..s + b]);
                from = s + b + close.len();
            }
            None => break,
        }
    }
    out
}

/// 解析 `.nksc` 边车，取出其中内嵌的配方载荷。
///
/// 边车是明文 XMP，`id="nikon::PictureControl"` 的滤镜内含 `ExportData` 与
/// `CustomData` 两个 Base64 字段。**边车记录的是 NX Studio 的编辑状态，未必等于
/// 拍摄时所用配方**，因此身份仍以照片的 MakerNote 为权威，边车只作本体来源。
pub fn parse_sidecar(text: &str, path: &Path) -> Vec<Recipe> {
    let xml = unescape(text);
    let mut out = Vec::new();
    for (tag, field) in [("ExportData", "ExportData"), ("CustomData", "CustomData")] {
        for content in tag_contents(&xml, tag) {
            let bytes = b64_decode(content);
            if bytes.len() < 4 || &bytes[0..4] != np3::MAGIC {
                continue;
            }
            if let Ok(c) = np3::parse(&bytes) {
                let name = c.name.clone();
                out.push(Recipe {
                    name,
                    source: Source::Sidecar { path: path.to_path_buf(), field: field.to_string() },
                    container: c,
                });
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Windows 注册表
// ---------------------------------------------------------------------------

/// 读取注册表配方来源。
///
/// 通过 `reg.exe` 读取以保持零依赖；沙箱与受限环境下失败时返回空列表而非报错，
/// 由调用方决定是否视为致命。
#[cfg(windows)]
pub fn scan_registry() -> Vec<Recipe> {
    use std::process::Command;
    const KEY: &str = r"HKCU\Software\Nikon\Common\PictureControl\NP3\CustomCurves";
    let out = match Command::new("reg").args(["query", KEY, "/s"]).output() {
        Ok(o) if o.status.success() => o,
        _ => return Vec::new(),
    };
    let text = String::from_utf8_lossy(&out.stdout);
    parse_reg_output(&text)
}

/// 非 Windows 平台无此来源。
#[cfg(not(windows))]
pub fn scan_registry() -> Vec<Recipe> {
    Vec::new()
}

/// 解析 `reg query /s` 的输出。
///
/// 每行形如 `    50331649    REG_BINARY    4E435000...`。
pub fn parse_reg_output(text: &str) -> Vec<Recipe> {
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 || !f[1].eq_ignore_ascii_case("REG_BINARY") {
            continue;
        }
        let hex: String = f[2..].join("");
        let Some(bytes) = hex_decode(&hex) else { continue };
        if bytes.len() < 4 || &bytes[0..4] != np3::MAGIC {
            continue;
        }
        if let Ok(c) = np3::parse(&bytes) {
            out.push(Recipe {
                name: c.name.clone(),
                source: Source::Registry { value_name: f[0].to_string() },
                container: c,
            });
        }
    }
    out
}

fn hex_decode(s: &str) -> Option<Vec<u8>> {
    if s.len() % 2 != 0 {
        return None;
    }
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(s.len() / 2);
    let d = |c: u8| -> Option<u8> {
        match c {
            b'0'..=b'9' => Some(c - b'0'),
            b'a'..=b'f' => Some(c - b'a' + 10),
            b'A'..=b'F' => Some(c - b'A' + 10),
            _ => None,
        }
    };
    for i in (0..b.len()).step_by(2) {
        out.push((d(b[i])? << 4) | d(b[i + 1])?);
    }
    Some(out)
}

// ---------------------------------------------------------------------------
// 目录扫描
// ---------------------------------------------------------------------------

/// 递归扫描目录下的 `.nksc` 边车。
pub fn scan_sidecars(root: &Path) -> Vec<Recipe> {
    let mut out = Vec::new();
    collect(root, &mut out, 0);
    out
}

fn collect(dir: &Path, out: &mut Vec<Recipe>, depth: usize) {
    if depth > 4 {
        return;
    }
    let Ok(rd) = std::fs::read_dir(dir) else { return };
    for e in rd.flatten() {
        let p = e.path();
        if p.is_dir() {
            collect(&p, out, depth + 1);
        } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case("nksc")) {
            // 用 lossy 读取：边车是 XMP 文本，但不应因个别非 UTF-8 字节整份丢弃
            if let Ok(bytes) = std::fs::read(&p) {
                let text = String::from_utf8_lossy(&bytes);
                out.extend(parse_sidecar(text.as_ref(), &p));
            }
        }
    }
}

/// 汇总全部本机来源。
pub fn scan_all(sidecar_roots: &[PathBuf]) -> Library {
    let mut lib = Library::new();
    for r in scan_registry() {
        lib.push(r);
    }
    for root in sidecar_roots {
        for r in scan_sidecars(root) {
            lib.push(r);
        }
    }
    lib
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_roundtrip_known_vector() {
        // "NCP\0" 的 Base64 前缀
        assert_eq!(&b64_decode("TkNQAAAA")[..4], b"NCP\0");
    }

    #[test]
    fn base64_tolerates_whitespace_only() {
        assert_eq!(b64_decode("TkNQAAAA"), b64_decode("Tk NQ\r\nAA\tAA"));
    }

    /// 未反转义的实体必须**不能**被当作数据吞掉——它含合法 Base64 字符。
    /// 这正是解析边车前必须先 `unescape` 的原因。
    #[test]
    fn base64_stops_at_unescaped_entities() {
        let with_entity = b64_decode("TkNQ&#xA;AAAA");
        let clean = b64_decode("TkNQAAAA");
        assert_ne!(with_entity, clean, "实体不得被静默吞掉");
        // 反转义后应恢复
        assert_eq!(b64_decode(&unescape("TkNQ&#xA;AAAA")), clean);
    }

    #[test]
    fn hex_decode_basic() {
        assert_eq!(hex_decode("4E435000").unwrap(), b"NCP\0");
        assert_eq!(hex_decode("4e4350").unwrap(), b"NCP");
        assert!(hex_decode("4E4").is_none());
    }

    #[test]
    fn parses_reg_query_output() {
        // 一段真实结果的缩略：值名是数字 ID，配方名在载荷内部
        let mut payload = Vec::new();
        payload.extend_from_slice(b"NCP\0");
        payload.extend_from_slice(&0x0000_0100u32.to_be_bytes());
        payload.extend_from_slice(&4u32.to_be_bytes());
        payload.extend_from_slice(b"0310");
        payload.extend_from_slice(&0x0000_0200u32.to_be_bytes());
        payload.extend_from_slice(&20u32.to_be_bytes());
        let mut nm = vec![0u8; 20];
        nm[..5].copy_from_slice(b"TestX");
        payload.extend_from_slice(&nm);
        payload.extend_from_slice(&0u32.to_be_bytes());
        payload.extend_from_slice(&0u32.to_be_bytes());
        let hex: String = payload.iter().map(|b| format!("{b:02X}")).collect();
        let text = format!(
            "HKEY_CURRENT_USER\\Software\\Nikon\\Common\\PictureControl\\NP3\\CustomCurves\r\n    \
             50331649    REG_BINARY    {hex}\r\n"
        );
        let got = parse_reg_output(&text);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].name, "TestX");
        assert_eq!(got[0].source, Source::Registry { value_name: "50331649".into() });
    }

    #[test]
    fn sidecar_extracts_export_and_custom() {
        // 构造一个最小边车：同样的载荷放在 ExportData 与 CustomData
        let mut payload = Vec::new();
        payload.extend_from_slice(b"NCP\0");
        payload.extend_from_slice(&0x0000_0100u32.to_be_bytes());
        payload.extend_from_slice(&4u32.to_be_bytes());
        payload.extend_from_slice(b"0310");
        payload.extend_from_slice(&0x0000_0200u32.to_be_bytes());
        payload.extend_from_slice(&20u32.to_be_bytes());
        let mut nm = vec![0u8; 20];
        nm[..4].copy_from_slice(b"Side");
        payload.extend_from_slice(&nm);
        payload.extend_from_slice(&0u32.to_be_bytes());
        payload.extend_from_slice(&0u32.to_be_bytes());
        let b64 = b64_encode(&payload);
        let text = format!(
            "<x>&lt;ExportData&gt;{b64}&lt;/ExportData&gt;&lt;Export&gt;\
             &lt;CustomData&gt;{b64}&lt;/CustomData&gt;</x>"
        );
        let got = parse_sidecar(&text, Path::new("a.nksc"));
        assert_eq!(got.len(), 2);
        assert!(got.iter().all(|r| r.name == "Side"));
    }

    #[test]
    fn library_find_is_case_insensitive() {
        let mut lib = Library::new();
        lib.push(Recipe {
            name: "LINKS-Nature".into(),
            source: Source::File { path: PathBuf::from("x") },
            container: np3::parse(&{
                let mut p = Vec::new();
                p.extend_from_slice(b"NCP\0");
                p.extend_from_slice(&1u32.to_be_bytes());
                p.extend_from_slice(&4u32.to_be_bytes());
                p.extend_from_slice(b"0310");
                p.extend_from_slice(&0u32.to_be_bytes());
                p.extend_from_slice(&0u32.to_be_bytes());
                p
            })
            .unwrap(),
        });
        assert_eq!(lib.find("links-nature").len(), 1);
        assert_eq!(lib.find("  LINKS-NAture ").len(), 1);
        assert!(lib.find("nope").is_empty());
        assert!(lib.require("nope").is_err());
    }

    /// 测试用编码器。
    pub(crate) fn b64_encode(data: &[u8]) -> String {
        const T: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut s = String::new();
        for c in data.chunks(3) {
            let b = [c[0], *c.get(1).unwrap_or(&0), *c.get(2).unwrap_or(&0)];
            let n = ((b[0] as u32) << 16) | ((b[1] as u32) << 8) | b[2] as u32;
            s.push(T[(n >> 18) as usize & 63] as char);
            s.push(T[(n >> 12) as usize & 63] as char);
            s.push(if c.len() > 1 { T[(n >> 6) as usize & 63] as char } else { '=' });
            s.push(if c.len() > 2 { T[n as usize & 63] as char } else { '=' });
        }
        s
    }
}
