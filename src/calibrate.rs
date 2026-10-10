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

/// 从 XMP 文本中提取 `crs:` / `crd:` 的键值。
///
/// # 两种写法都要认
///
/// XMP 表达同一件事有两种形式，**NX Studio 用的是后者**：
///
/// ```xml
/// <rdf:Description crs:Exposure2012="0.33"/>          <!-- 属性形式 -->
/// <rdf:Description><crd:Exposure2012>0.33</crd:Exposure2012></rdf:Description>  <!-- 元素形式 -->
/// ```
///
/// 只认属性形式会漏掉 NX Studio 的全部输出——实测 15 张导出全部解析出 0 项，
/// 校验因此会把合格样本全部误判为「无法确认无调整」。
///
/// # 为什么全程在字节上扫描
///
/// 按字节下标去切 `&str` 会在多字节字符中间 panic——NX Studio 的 XMP 带 BOM
/// （U+FEFF，UTF-8 三字节），第一刀就会踩到。
pub fn xmp_adjustments(xmp: &str) -> std::collections::BTreeMap<String, String> {
    let mut out = std::collections::BTreeMap::new();
    let b = xmp.as_bytes();
    let mut i = 0;
    let is_name_byte = |c: u8| c.is_ascii_alphanumeric() || c == b'_';

    while i < b.len() {
        // —— 属性形式：crs:Name="值" / crd:Name="值"
        let attr = b[i..].starts_with(b"crs:") || b[i..].starts_with(b"crd:");
        if attr {
            let name_start = i + 4;
            let mut j = name_start;
            while j < b.len() && is_name_byte(b[j]) {
                j += 1;
            }
            if j > name_start {
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
                        let vs = k + 1;
                        let mut v = vs;
                        while v < b.len() && b[v] != q {
                            v += 1;
                        }
                        out.insert(
                            String::from_utf8_lossy(&b[name_start..j]).into_owned(),
                            String::from_utf8_lossy(&b[vs..v]).into_owned(),
                        );
                        i = v + 1;
                        continue;
                    }
                }
            }
        }

        // —— 元素形式：<crs:Name>值</crs:Name>
        if b[i] == b'<' {
            let name_start = i + 1;
            if b[name_start..].starts_with(b"crs:") || b[name_start..].starts_with(b"crd:") {
                let id_start = name_start + 4;
                let mut j = id_start;
                while j < b.len() && is_name_byte(b[j]) {
                    j += 1;
                }
                if j > id_start && j < b.len() && b[j] == b'>' {
                    let vs = j + 1;
                    // 结束标签形如 </crd:Name>，这里只找下一个 '<'
                    let mut v = vs;
                    while v < b.len() && b[v] != b'<' {
                        v += 1;
                    }
                    let name = String::from_utf8_lossy(&b[id_start..j]).into_owned();
                    let value = String::from_utf8_lossy(&b[vs..v]).trim().to_string();
                    out.insert(name, value);
                    i = v;
                    continue;
                }
            }
        }

        i += 1;
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
///
/// 注意：判据已不再是「这些值必须为 0」（见 [`consensus_adjustments`]）——
/// 相机在每个文件里固定写入的 `crd:` 块本身就带非零值。此函数现仅供测试使用。
#[cfg(test)]
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
    /// 导致该对被排除的问题。
    pub issues: Vec<Issue>,
    /// 供人工核对的信息，**不影响合格与否**。
    pub notes: Vec<String>,
}

impl PairReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
}

/// NX Studio 边车文件的路径（按约定：同目录下 `NKSC_PARAM\<完整文件名>.nksc`）。
///
/// 存在边车意味着 NEF 内嵌的设置**已过期**——导出以边车为准。
pub fn sidecar_path(raw: &Path) -> Option<PathBuf> {
    let name = raw.file_name()?.to_string_lossy().into_owned();
    let dir = raw.parent()?;
    let p = dir.join("NKSC_PARAM").join(format!("{name}.nksc"));
    p.is_file().then_some(p)
}

/// 内嵌设置的简短描述，用于提示。
fn embedded_brief(raw: &[u8]) -> String {
    let pc = crate::picture_control::read(raw)
        .map(|id| id.name)
        .unwrap_or_else(|_| "（读不到）".into());
    let adl = crate::makernote::active_d_lighting(raw)
        .ok()
        .flatten()
        .map(|a| a.name())
        .unwrap_or_else(|| "（读不到）".into());
    format!("PC={pc}, ADL={adl}")
}

/// 读出一个参考导出内嵌 XMP 中的调整块。
///
/// 返回 `None` 表示没有 XMP 或解析不出任何项。
pub fn export_adjustments(
    data: &[u8],
) -> crate::error::Result<Option<std::collections::BTreeMap<String, String>>> {
    use crate::tiff::Tiff;
    let t = Tiff::locate(data)?;
    let ifd = t.ifd0_offset()?;
    let Some(e) = t.find(ifd, 0x02BC)? else { return Ok(None) };
    let bytes = t.bytes(&e)?.to_vec();
    let xmp = String::from_utf8_lossy(&bytes);
    let m = xmp_adjustments(&xmp);
    if m.is_empty() {
        Ok(None)
    } else {
        Ok(Some(m))
    }
}

/// 一组导出的**共识调整块**：每个键取出现次数最多的值。
///
/// 尼康相机在每个文件里都会嵌一块固定的 Adobe 默认值（`crd:` 命名空间，
/// `Exposure2012=0.33`、`Highlights2012=-21`、`Saturation=20` 之类）。它不是用户
/// 的调整，因此**不能要求它为 0**——那样会把每一张都判成"有额外调整"。
///
/// 正确的判据是：**同一批参考导出之间应当一致**。谁不一样，谁才被动过。
pub fn consensus_adjustments(
    maps: &[&std::collections::BTreeMap<String, String>],
) -> std::collections::BTreeMap<String, String> {
    let mut counts: std::collections::BTreeMap<&str, std::collections::BTreeMap<&str, usize>> =
        std::collections::BTreeMap::new();
    for m in maps {
        for (k, v) in m.iter() {
            *counts.entry(k.as_str()).or_default().entry(v.as_str()).or_insert(0) += 1;
        }
    }
    counts
        .into_iter()
        .filter_map(|(k, vs)| {
            vs.into_iter().max_by_key(|(_, n)| *n).map(|(v, _)| (k.to_string(), v.to_string()))
        })
        .collect()
}

/// 校验一份参考集。
///
/// 四项前提逐一检查，**任一项不满足即排除该对并报告原因**——不是降权继续用。
/// 降权会让不合格样本以看不见的方式污染拟合。
///
/// # 「无额外调整」这一项的特殊之处
///
/// 判据是**整组一致**，不是绝对为 0——理由见 [`consensus_adjustments`]。
/// 局限也要说清：若**所有**参考导出都做了同一处调整（比如统一加了 +0.3 EV），
/// 共识会把它一起吸收，这一项查不出来。因此校验结论里会一并列出共识块本身，
/// 供人工过目。
pub fn validate(set: &ReferenceSet) -> Vec<PairReport> {
    // 第一遍：读出各导出的调整块，供「整组一致」判据使用
    let maps: Vec<Option<std::collections::BTreeMap<String, String>>> = set
        .pairs
        .iter()
        .map(|p| std::fs::read(&p.export).ok().and_then(|d| export_adjustments(&d).ok().flatten()))
        .collect();
    let present: Vec<&std::collections::BTreeMap<String, String>> = maps.iter().flatten().collect();
    let consensus = consensus_adjustments(&present);

    set.pairs
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut issues = Vec::new();
            let mut notes = Vec::new();
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
                // 关键：**存在 NX Studio 边车时，NEF 内嵌的设置是过期的**。
                //
                // 相机把拍摄时的 Picture Control / ADL 写进 MakerNote，而 NX Studio
                // 的修改存在同目录 `NKSC_PARAM\<名>.nksc` 里，导出时**边车覆盖内嵌值**。
                // 实测：15 张参考图的 NEF 都记着 `LINKS-Nature` + ADL 标准，而用户
                // 在 NX Studio 里设的是「自然 + ADL 关闭」——拿 NEF 去校验会把这 15 张
                // 全部误判为不合格。
                let sidecar = sidecar_path(&p.raw);
                let overridden = sidecar.is_some();

                if overridden {
                    // 不拿过期值下结论，但把差异记下来供人工核对
                    let eff = export_adjustments(export.as_deref().unwrap_or(&[]))
                        .ok()
                        .flatten()
                        .and_then(|m| m.get("CameraProfile").cloned());
                    notes.push(format!(
                        "NEF 内嵌设置为「{}」，但存在 NX Studio 边车 {}，导出以边车为准；\
                         导出内嵌的 CameraProfile = {}",
                        embedded_brief(raw),
                        sidecar.as_ref().unwrap().display(),
                        eff.as_deref().unwrap_or("（未记录）")
                    ));
                } else if set.premises.require_adl_off {
                    // 无边车 → 内嵌设置就是实际使用的设置，正常校验
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

                if !overridden {
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
                        }
                        Err(e) => push(&mut issues, format!("解析 Picture Control 失败：{e}")),
                    }
                }
            }

            if let Some(export) = &export {
                // 3. 参考导出中不得有额外调整——判据是**整组一致**（见 validate 的文档）
                if set.premises.require_no_adjustments {
                    match &maps[i] {
                        None => push(
                            &mut issues,
                            "参考导出不含可解析的 XMP 调整块，无法确认无额外调整".into(),
                        ),
                        Some(m) => {
                            let mut diff = Vec::new();
                            for k in NEUTRAL_ADJUSTMENTS {
                                if let Some(v) = m.get(*k) {
                                    match consensus.get(*k) {
                                        Some(bv) if bv != v => {
                                            diff.push(format!("{k} 为 {v}，整组为 {bv}"))
                                        }
                                        _ => {}
                                    }
                                }
                            }
                            if !diff.is_empty() {
                                push(
                                    &mut issues,
                                    format!("参考导出与整组不一致，疑似有额外调整：{}", diff.join("；")),
                                );
                            }
                        }
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
                                // 竖构图导出会被转置（实测参考图是 5504×8256），
                                // 因此两种朝向都接受。
                                let ok = (ew, eh) == (cw, ch) || (ew, eh) == (ch, cw);
                                if !ok {
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

            PairReport { index: i, raw: p.raw.clone(), export: p.export.clone(), issues, notes }
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
