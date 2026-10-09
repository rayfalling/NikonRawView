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

    /// 内容校验是否已实现。在它变真之前，拟合不应被当作可信。
    pub const CONTENT_VALIDATION_IMPLEMENTED: bool = false;

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }
}

/// 从清单文件读入参考集。
///
/// 尚未实现——格式与任务 5.1 一起定。
pub fn load_manifest(_path: &Path) -> crate::error::Result<ReferenceSet> {
    Err(crate::error::Error::InvalidInput(
        "参考导出清单的读取尚未实现（见任务 5.1）".into(),
    ))
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
            !ReferenceSet::CONTENT_VALIDATION_IMPLEMENTED,
            "内容校验已实现——请同步更新本模块文档、Premises 与任务 5.1 的状态"
        );
    }

    #[test]
    fn manifest_loading_is_explicitly_unimplemented() {
        let e = load_manifest(Path::new("x"));
        assert!(e.is_err(), "未实现时必须是明确失败，而不是返回空集");
    }
}
