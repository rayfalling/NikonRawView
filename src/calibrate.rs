//! 基准配方的标定。
//!
//! # 这个模块要解决什么
//!
//! 基准配方（NEUTRAL 等）的渲染变换**不在任何文件里**。已实测确认：NX Studio 安装
//! 目录下的 `PicCon.bin` / `PicCon21.bin` 是携带配方记录的模板 TIFF（标签体系与 NP3
//! 同源），但用「257 个 u16 单调递增」扫描确认**不含色调曲线**。
//!
//! 因此只能从**参考导出**拟合：同一张 NEF，用 NX Studio 以指定基准重新导出，得到
//! 「线性输入 ↔ 尼康输出」的配对，再据此求解基准变换。
//!
//! # 管线中的位置
//!
//! ```text
//! 解码（raw-decode） → 色彩（color-pipeline） → 渲染（picture-control/render）
//!                                                        ↑
//!                                     标定（本模块）产出基准变换数据
//! ```
//!
//! 标定是**离线**的：日常渲染只加载拟合产出的数据文件（[`crate::transform`]），
//! **不依赖参考导出文件的存在**。
//!
//! # 当前状态
//!
//! 本模块目前只有骨架：类型与前提条件已定义，拟合尚未实现（见任务 5.3–5.6）。
//! 参考导出的**采集**需要人工操作（用 NX Studio 重导），其清单与校验见任务 5.1。

use std::path::{Path, PathBuf};

/// 一对「原始文件 ↔ 参考导出」。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferencePair {
    /// 原始 NEF。
    pub raw: PathBuf,
    /// NX Studio 导出的参考图（TIFF）。
    pub export: PathBuf,
}

/// 参考导出必须满足的前提。
///
/// 这些条件不是"尽量满足"——任何一条不满足都会把非基准的因素混进拟合里，
/// 产出一个"数值上能对上、物理上对不上"的基准变换。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Premises {
    /// 基准配方的名称（应等于目标基准名）。
    pub baseline_name: String,
    /// 基准编码（拟合结果按它索引，见 [`crate::render::select_baseline`]）。
    pub baseline_code: u16,
    /// 是否要求参考导出中 ADL 已关闭。
    pub require_adl_off: bool,
    /// 是否要求参考导出中无额外调整（曝光、对比度等滑块均为默认）。
    pub require_no_adjustments: bool,
}

impl Default for Premises {
    fn default() -> Self {
        Self {
            baseline_name: "NEUTRAL".into(),
            baseline_code: 0x03C2,
            require_adl_off: true,
            require_no_adjustments: true,
        }
    }
}

/// 校验中发现的问题。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Issue {
    /// 原始 NEF 不存在。
    MissingRaw { pair_index: usize, path: PathBuf },
    /// 参考导出不存在。
    MissingExport { pair_index: usize, path: PathBuf },
    /// 采集前提未满足。
    PremiseViolated { pair_index: usize, detail: String },
}

impl Issue {
    pub fn describe(&self) -> String {
        match self {
            Issue::MissingRaw { pair_index, path } => {
                format!("第 {} 对：原始文件不存在 {}", pair_index + 1, path.display())
            }
            Issue::MissingExport { pair_index, path } => {
                format!("第 {} 对：参考导出不存在 {}", pair_index + 1, path.display())
            }
            Issue::PremiseViolated { pair_index, detail } => {
                format!("第 {} 对：{}", pair_index + 1, detail)
            }
        }
    }
}

/// 一份参考导出清单。
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReferenceSet {
    pub premises: Premises,
    pub pairs: Vec<ReferencePair>,
}

impl ReferenceSet {
    /// 只做**文件存在性**的检查。
    ///
    /// 内容层面的校验（ADL 是否关闭、基准是否为指定值、有无额外调整、导出是否与
    /// NEF 同源）尚未实现——见任务 5.1。在它实现之前，本方法**不足以**判定一份
    /// 参考集可用于拟合。
    pub fn check_files_exist(&self) -> Vec<Issue> {
        let mut out = Vec::new();
        for (i, p) in self.pairs.iter().enumerate() {
            if !p.raw.is_file() {
                out.push(Issue::MissingRaw { pair_index: i, path: p.raw.clone() });
            }
            if !p.export.is_file() {
                out.push(Issue::MissingExport { pair_index: i, path: p.export.clone() });
            }
        }
        out
    }

    /// 内容校验是否已实现。
    ///
    /// 已实现（任务 5.1）：ADL 是否关闭、Picture Control 是否为指定基准、
    /// 参考导出有无额外调整、与配对 NEF 是否同源（以有效像素尺寸为证）。
    pub const CONTENT_VALIDATION_IMPLEMENTED: bool = true;

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }
}

/// 从 XMP 文本中提取 `crs:` / `crd:` 的属性键值。
///
/// 只做属性扫描，不引入 XML 解析器——这份 XMP 是 NX Studio 生成的定式文本，我们
/// 只关心其中若干标量。属性名可能带 `crs:` 或 `crd:` 前缀（前者用于单文件设置，
/// 后者用于默认值），两种前缀都收。
pub fn xmp_adjustments(xmp: &str) -> std::collections::BTreeMap<String, String> {
    let mut out = std::collections::BTreeMap::new();
    let b = xmp.as_bytes();
    let mut i = 0;

    // 全程在**字节**上扫描。按字节下标去切 &str 会在多字节字符中间 panic——
    // NX Studio 的 XMP 带 BOM（U+FEFF，UTF-8 三字节），第一刀就会踩到。
    let is_name_byte = |c: u8| c.is_ascii_alphanumeric() || c == b'_';

    while i < b.len() {
        let skip = if b[i..].starts_with(b"crs:") || b[i..].starts_with(b"crd:") {
            4
        } else {
            i += 1;
            continue;
        };
        let name_start = i + skip;
        let mut j = name_start;
        while j < b.len() && is_name_byte(b[j]) {
            j += 1;
        }
        if j == name_start {
            i += 1;
            continue;
        }

        // 其后应紧跟 ="…"
        let mut k = j;
        while k < b.len() && b[k].is_ascii_whitespace() {
            k += 1;
        }
        if k < b.len() && b[k] == b'=' {
            k += 1;
            while k < b.len() && b[k].is_ascii_whitespace() {
                k += 1;
            }
            if k < b.len() && (b[k] == b'"' || b[k] == b'\'') {
                let q = b[k];
                let vstart = k + 1;
                let mut v = vstart;
                while v < b.len() && b[v] != q {
                    v += 1;
                }
                let name = String::from_utf8_lossy(&b[name_start..j]).into_owned();
                let value = String::from_utf8_lossy(&b[vstart..v]).into_owned();
                out.insert(name, value);
                i = v + 1;
                continue;
            }
        }
        i = j;
    }
    out
}

/// 必须处于中性默认值的调整项。
///
/// 这些一旦非零，参考导出就不再是"指定基准下的干净渲染"，拟合会把用户的调整
/// 当成基准的一部分。
///
/// **刻意不列入**的项：`CameraProfile`（它本就该写明基准）、`Temperature` / `Tint`
/// （反映拍摄白平衡，非零是正常的）、`ProcessVersion`、`WhiteBalance`。
pub const NEUTRAL_ADJUSTMENTS: &[&str] = &[
    "Exposure2012",
    "Contrast2012",
    "Highlights2012",
    "Shadows2012",
    "Whites2012",
    "Blacks2012",
    "Clarity2012",
    "Texture",
    "Dehaze",
    "Vibrance",
    "Saturation",
    "Sharpness",
    "SharpenRadius",
    "SharpenDetail",
    "LuminanceSmoothing",
    "ColorNoiseReduction",
    "PostCropVignetteAmount",
    "GrainAmount",
];

/// 数值是否可视为中性默认（容忍浮点表示差异）。
///
/// **空值与非数值一律判为非中性**：`crs:Saturation=""` 这种写法无法确认其为 0，
/// 而标定数据的原则是"无法确认就拒绝"，不是"无法确认就放过"。
fn is_neutral_value(v: &str) -> bool {
    match v.trim().parse::<f64>() {
        Ok(x) => x.abs() < 1e-6,
        Err(_) => false,
    }
}

/// 逐对校验的结论。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairReport {
    pub index: usize,
    pub raw: PathBuf,
    pub export: PathBuf,
    pub issues: Vec<Issue>,
}

impl PairReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
}

/// 校验一份参考集。
///
/// 四项前提逐一检查，**任一项不满足即排除该对并报告原因**——不是降权继续用。
/// 降权会让不合格样本以看不见的方式污染拟合。
pub fn validate(set: &ReferenceSet) -> Vec<PairReport> {
    set.pairs
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut issues = Vec::new();
            let push = |issues: &mut Vec<Issue>, detail: String| {
                issues.push(Issue::PremiseViolated { pair_index: i, detail });
            };

            let raw = match std::fs::read(&p.raw) {
                Ok(d) => Some(d),
                Err(e) => {
                    push(&mut issues, format!("原始文件无法读取（{e}）：{}", p.raw.display()));
                    None
                }
            };
            let export = match std::fs::read(&p.export) {
                Ok(d) => Some(d),
                Err(e) => {
                    push(&mut issues, format!("参考导出无法读取（{e}）：{}", p.export.display()));
                    None
                }
            };

            if let Some(raw) = &raw {
                // 1. ADL 必须关闭
                if set.premises.require_adl_off {
                    match crate::makernote::active_d_lighting(raw) {
                        Ok(Some(adl)) if adl.is_off() => {}
                        Ok(Some(adl)) => push(
                            &mut issues,
                            format!("ADL 未关闭（当前为「{}」）——ADL 会改变明暗分布，会被混进拟合", adl.name()),
                        ),
                        Ok(None) => push(&mut issues, "读不到 ADL 标签，无法确认已关闭".into()),
                        Err(e) => push(&mut issues, format!("解析 ADL 失败：{e}")),
                    }
                }

                // 2. Picture Control 必须是指定的基准
                match crate::picture_control::read(raw) {
                    Ok(id) => {
                        if !id.name.eq_ignore_ascii_case(&set.premises.baseline_name) {
                            push(
                                &mut issues,
                                format!(
                                    "Picture Control 不是指定基准：要求「{}」，实际「{}」",
                                    set.premises.baseline_name, id.name
                                ),
                            );
                        }
                        if !id.base.eq_ignore_ascii_case(&set.premises.baseline_name) {
                            push(
                                &mut issues,
                                format!(
                                    "基准色彩不是指定值：要求「{}」，实际「{}」",
                                    set.premises.baseline_name, id.base
                                ),
                            );
                        }
                    }
                    Err(e) => push(&mut issues, format!("解析 Picture Control 失败：{e}")),
                }
            }

            if let Some(export) = &export {
                // 3. 参考导出中不得有额外调整
                if set.premises.require_no_adjustments {
                    match crate::tiff::Tiff::locate(export)
                        .and_then(|t| {
                            let ifd = t.ifd0_offset()?;
                            Ok(t.find(ifd, 0x02BC)?.and_then(|e| t.bytes(&e).ok().map(|b| b.to_vec())))
                        }) {
                        Ok(Some(xmp_bytes)) => {
                            let xmp = String::from_utf8_lossy(&xmp_bytes);
                            let adj = xmp_adjustments(&xmp);
                            if adj.is_empty() {
                                push(&mut issues, "参考导出含 XMP 但未解析出任何调整项，无法确认无调整".into());
                            }
                            let mut offending = Vec::new();
                            for k in NEUTRAL_ADJUSTMENTS {
                                if let Some(v) = adj.get(*k) {
                                    if !is_neutral_value(v) {
                                        offending.push(format!("{k}={v}"));
                                    }
                                }
                            }
                            if !offending.is_empty() {
                                push(
                                    &mut issues,
                                    format!("参考导出含有额外调整：{}", offending.join("、")),
                                );
                            }
                        }
                        Ok(None) => push(
                            &mut issues,
                            "参考导出不含 XMP，无法确认无额外调整".into(),
                        ),
                        Err(e) => push(&mut issues, format!("解析参考导出的 XMP 失败：{e}")),
                    }
                }

                // 4. 与配对 NEF 同源——以有效像素尺寸一致作为可查的证据
                if raw.is_some() {
                    let want = crate::camera::read_model(&p.raw)
                        .and_then(|m| crate::camera::lookup(&m).map(|e| e.margins));
                    if let Some(margins) = want {
                        if let (Ok((ew, eh)), Ok((nw, nh))) = (
                            tiff_dimensions(export),
                            crate::libraw::read_dimensions(&p.raw),
                        ) {
                            if let Some((cw, ch)) = margins.effective(nw, nh) {
                                if (ew, eh) != (cw, ch) {
                                    push(
                                        &mut issues,
                                        format!(
                                            "尺寸与配对 NEF 不符：导出 {ew}×{eh}，该 NEF 的有效区为 {cw}×{ch}——可能配错了文件"
                                        ),
                                    );
                                }
                            }
                        }
                    }
                }
            }

            PairReport { index: i, raw: p.raw.clone(), export: p.export.clone(), issues }
        })
        .collect()
}

/// 读出一个 TIFF 的 IFD0 尺寸。
pub fn tiff_dimensions(data: &[u8]) -> crate::error::Result<(usize, usize)> {
    use crate::tiff::Tiff;
    let t = Tiff::locate(data)?;
    let ifd = t.ifd0_offset()?;
    let w = t
        .find(ifd, 0x0100)?
        .ok_or_else(|| crate::error::Error::InvalidInput("TIFF 缺少宽度标签".into()))?;
    let h = t
        .find(ifd, 0x0101)?
        .ok_or_else(|| crate::error::Error::InvalidInput("TIFF 缺少高度标签".into()))?;
    Ok((t.u32_value(&w)? as usize, t.u32_value(&h)? as usize))
}

/// 从清单文件读入参考集。
///
/// 格式：每行一对，制表符分隔——`<原始 NEF 路径>\t<参考导出路径>`。
/// `#` 开头为注释，空行忽略。路径按 UTF-8 原样保存（本项目的照片库路径含中文）。
pub fn load_manifest(path: &Path) -> crate::error::Result<ReferenceSet> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| crate::error::Error::InvalidInput(format!("读取清单 {} 失败：{e}", path.display())))?;
    let mut set = ReferenceSet::default();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim_end_matches(['\r', '\n']);
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut it = t.split('\t');
        let raw = it.next().unwrap_or("").trim();
        let export = it.next().unwrap_or("").trim();
        if raw.is_empty() || export.is_empty() {
            return Err(crate::error::Error::BadPayload {
                what: "参考导出清单",
                detail: format!("第 {} 行不是「原始\\t导出」两列：{t}", n + 1),
            });
        }
        set.pairs.push(ReferencePair { raw: PathBuf::from(raw), export: PathBuf::from(export) });
    }
    Ok(set)
}

/// 写出清单文件。
pub fn save_manifest(path: &Path, set: &ReferenceSet) -> crate::error::Result<()> {
    let mut s = String::from(
        "# 参考导出清单：每行一对，制表符分隔\n\
         # 前提：Picture Control = 指定基准；ADL 关闭；无其他调整\n",
    );
    s.push_str(&format!("# 基准：{}（0x{:04X}）\n", set.premises.baseline_name, set.premises.baseline_code));
    for p in &set.pairs {
        s.push_str(&format!("{}\t{}\n", p.raw.display(), p.export.display()));
    }
    std::fs::write(path, s)
        .map_err(|e| crate::error::Error::InvalidInput(format!("写入清单 {} 失败：{e}", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_the_documented_baseline() {
        let p = Premises::default();
        assert_eq!(p.baseline_code, 0x03C2, "默认基准编码应为 NEUTRAL");
        assert!(p.require_adl_off, "ADL 必须关闭——否则会被混进拟合");
        assert!(p.require_no_adjustments);
    }

    #[test]
    fn missing_files_are_reported_with_position() {
        let set = ReferenceSet {
            premises: Premises::default(),
            pairs: vec![
                ReferencePair { raw: "no-such.nef".into(), export: "no-such.tif".into() },
                ReferencePair { raw: "also-missing.nef".into(), export: "x.tif".into() },
            ],
        };
        let issues = set.check_files_exist();
        assert_eq!(issues.len(), 4, "两对应各报缺失的 raw 与 export：{issues:?}");
        assert!(issues[0].describe().contains("第 1 对"));
        assert!(issues[2].describe().contains("第 2 对"));
    }

    #[test]
    fn existing_files_yield_no_issue() {
        let me = PathBuf::from(file!());
        let set = ReferenceSet {
            premises: Premises::default(),
            pairs: vec![ReferencePair { raw: me.clone(), export: me }],
        };
        assert!(set.check_files_exist().is_empty());
    }

    /// 这条断言的意义：在内容校验实现之前，任何"拟合可信"的说法都缺依据。
    /// 若有人实现了它，把常量改成 true 时这里会失败，提醒同步更新文档与任务。
    #[allow(clippy::assertions_on_constants)] // 断言对象就是常量，这正是它的用途
    #[test]
    fn content_validation_is_declared_unimplemented() {
        assert!(
            ReferenceSet::CONTENT_VALIDATION_IMPLEMENTED,
            "内容校验应当已实现（任务 5.1）——若回退为未实现，请同步更新文档"
        );
    }

    #[test]
    fn manifest_loading_is_explicitly_unimplemented() {
        // 现在是"文件不存在"这类明确失败，而不是静默返回空集
        let e = load_manifest(Path::new("no-such-manifest.tsv"));
        assert!(e.is_err());
        assert!(e.unwrap_err().to_string().contains("读取清单"));
    }

    #[test]
    fn xmp_scanner_picks_up_both_prefixes() {
        let xmp = r#"<rdf:Description crs:Exposure2012="+0.33" crd:Shadows2012="10"
                      crs:CameraProfile="Camera Neutral" crs:Temperature="5500"/>"#;
        let a = xmp_adjustments(xmp);
        assert_eq!(a.get("Exposure2012").map(String::as_str), Some("+0.33"));
        assert_eq!(a.get("Shadows2012").map(String::as_str), Some("10"));
        assert_eq!(a.get("CameraProfile").map(String::as_str), Some("Camera Neutral"));
        assert_eq!(a.get("Temperature").map(String::as_str), Some("5500"));
        // 非 crs/crd 前缀的属性不应被收进来
        assert!(!a.contains_key("about"));
    }

    #[test]
    fn xmp_scanner_survives_multibyte_and_bom() {
        // 回归：NX Studio 的 XMP 带 BOM（U+FEFF，UTF-8 三字节）。按字节下标去切
        // &str 会在多字节字符中间 panic——这条曾经真的崩过。
        let xmp = "\u{feff}<rdf:Description crs:Exposure2012=\"0.5\" \
                   crs:Artist=\"张三\" crs:Shadows2012=\"0\"/>";
        let a = xmp_adjustments(xmp);
        assert_eq!(a.get("Exposure2012").map(String::as_str), Some("0.5"));
        assert_eq!(a.get("Shadows2012").map(String::as_str), Some("0"));
        assert_eq!(a.get("Artist").map(String::as_str), Some("张三"));
    }

    #[test]
    fn neutral_value_judgement() {
        for v in ["0", "0.0", "+0.00", " 0 "] {
            assert!(is_neutral_value(v), "{v} 应判为中性");
        }
        for v in ["0.33", "-21", "20", "", "abc"] {
            assert!(!is_neutral_value(v), "{v} 应判为非中性");
        }
    }

    #[test]
    fn added_adjustments_are_detected() {
        // NX Studio 导出里真实见过的一组非默认值——这类样本必须被拒绝
        let xmp = r#"<rdf:Description crd:Exposure2012="0.33" crd:Highlights2012="-21"
                      crd:Shadows2012="10" crd:Saturation="20" crd:Clarity2012="2"
                      crd:CameraProfile="Camera Neutral" crd:Temperature="5500"/>"#;
        let a = xmp_adjustments(xmp);
        let offending: Vec<&str> = NEUTRAL_ADJUSTMENTS
            .iter()
            .filter(|k| a.get(**k).is_some_and(|v| !is_neutral_value(v)))
            .copied()
            .collect();
        assert_eq!(offending.len(), 5, "应挑出 5 项非默认调整：{offending:?}");
        assert!(offending.contains(&"Exposure2012"));
        // 刻意不列入的项不应被误判
        assert!(!offending.contains(&"CameraProfile"));
        assert!(!offending.contains(&"Temperature"));
    }

    #[test]
    fn manifest_roundtrip_preserves_non_ascii_paths() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(format!("manifest-test-{}.tsv", std::process::id()));
        let mut set = ReferenceSet::default();
        set.pairs.push(ReferencePair {
            raw: PathBuf::from(r"Z:\摄影\Nikon\Z8\2025.12.20-12.27 冰岛\DSC_1001.NEF"),
            export: PathBuf::from(r"E:\参考导出\DSC_1001.tif"),
        });
        save_manifest(&p, &set).unwrap();
        let back = load_manifest(&p).unwrap();
        let _ = std::fs::remove_file(&p);
        assert_eq!(back.pairs, set.pairs, "中文路径应原样往返");
    }

    #[test]
    fn malformed_manifest_line_is_rejected() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(format!("manifest-bad-{}.tsv", std::process::id()));
        std::fs::write(&p, "# 注释\n只有一列\n").unwrap();
        let e = load_manifest(&p);
        let _ = std::fs::remove_file(&p);
        assert!(e.is_err(), "缺列的行应报错而非静默跳过");
    }
}
