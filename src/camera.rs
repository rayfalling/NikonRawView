//! 机型知识表。
//!
//! # 为什么需要这张表
//!
//! **相机有效像素区的边距不在文件里。** 实测确认：LibRaw 对 Nikon Z8 的 NEF 输出
//! 传感器全幅 `8280×5520`（含光学黑/掩蔽边框），而 Z8 的有效像素区是 `8256×5504`——
//! 左右各差 12 像素、上下各差 8 像素。这个差值来自尼康的机型规格，不是元数据。
//!
//! 因此边距必须作为**显式的机型知识**维护：把它隐含在代码判断里，换个机型就会
//! 静默裁错尺寸、画错边框。机型不在表中时应当明确报告，而不是猜一个边距。

/// 某机型相对传感器全幅的内缩边距（单位：像素）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Margins {
    pub left: usize,
    pub top: usize,
    pub right: usize,
    pub bottom: usize,
}

impl Margins {
    pub const ZERO: Self = Self { left: 0, top: 0, right: 0, bottom: 0 };

    /// 全幅尺寸减去边距后的有效尺寸，越界时返回 `None`。
    pub fn effective(&self, full_w: usize, full_h: usize) -> Option<(usize, usize)> {
        let w = full_w.checked_sub(self.left + self.right)?;
        let h = full_h.checked_sub(self.top + self.bottom)?;
        (w > 0 && h > 0).then_some((w, h))
    }

    /// 转过 90° 之后的边距。
    ///
    /// 边距表是针对**未旋转的传感器画幅**定义的。解码层会按 EXIF 方向旋转输出，
    /// 竖拍照片因此需要把左右与上下对调——否则会切错边。
    ///
    /// 实测印证：Z8 横构图边距为 `左12 上8 右12 下8`，全幅 `8280×5520` 裁后得
    /// `8256×5504`；竖拍输出为 `5520×8280`，用对调后的 `左8 上12 右8 下12` 裁，
    /// 得 `5504×8256`——与尼康工坊导出的参考图**完全一致**。
    pub fn rotated(self) -> Self {
        Self {
            left: self.top,
            top: self.right,
            right: self.bottom,
            bottom: self.left,
        }
    }
}

/// 一条相机 → 工作空间的线性色彩矩阵，连同它的适用条件与来源。
///
/// # 为什么矩阵要存在这里，而不是向解码层索要
///
/// 解码层有一个「相机 → 输出空间」的访问接口（`libraw_get_rgb_cam`），但它返回的
/// 矩阵**与输出色彩空间无关**——实测把输出分别设为 Camera / sRGB / ProPhoto，该接口
/// 给出的矩阵**逐位相同**，而像素两两平均绝对差达 254.7–550.5（16 位刻度）。也就是
/// 说那不是实际参与转换的矩阵。
///
/// 因此本项目的做法是**测**出来：用同一文件在「相机空间」与「工作空间」两种输出下
/// 的配对做最小二乘（见 [`crate::color::derive_matrix`]）。实测残差落在 1e-5 量级，
/// 且**换一张照片复验仍然成立**（6,348,068 样本，平均绝对差 0.000007）——这同时证明
/// 矩阵按机型固定，所以每个机型只需测一次。
#[derive(Debug, Clone, Copy)]
pub struct ColorMatrixEntry {
    /// 行主序 3×3。
    pub matrix: [[f32; 3]; 3],
    /// 该矩阵的目的空间。
    pub output_space: &'static str,
    /// 推导来源——含所用参考照片与残差，便于复核与再现。
    pub source: &'static str,
}

impl ColorMatrixEntry {
    /// 各行之和。相机矩阵把中性映射为中性时，三者都应接近 1.0。
    ///
    /// 注意这只是**必要条件而非充分条件**：实测正确的矩阵与解码层给的那个错矩阵
    /// 行和都是 1.0，仅凭这一条分辨不出真假。真正的判据是跨图复验的残差。
    pub fn row_sums(&self) -> [f32; 3] {
        let mut out = [0f32; 3];
        for (i, row) in self.matrix.iter().enumerate() {
            out[i] = row.iter().sum();
        }
        out
    }
}

/// 一条机型记录。
#[derive(Debug, Clone, Copy)]
pub struct Model {
    /// 用于匹配的机型标识（小写、去空白）。LibRaw 的 `model` 形如 `NIKON Z 8`。
    pub match_key: &'static str,
    /// 可读名称。
    pub display: &'static str,
    pub margins: Margins,
    /// 该条边距的来源依据——便于日后复核，也避免"数字从哪来的"变成谜。
    pub source: &'static str,
    /// 相机 → 工作空间的色彩矩阵；未标定时为 `None`。
    pub color_matrix: Option<ColorMatrixEntry>,
}

/// 已知机型表。
///
/// 新增机型时**必须**填 `source`：说明这个边距是从哪里得到的（官方规格、
/// 与尼康工坊导出对比、或其他可复核的途径）。
pub const MODELS: &[Model] = &[Model {
    match_key: "nikon z 8",
    display: "Nikon Z 8",
    margins: Margins { left: 12, top: 8, right: 12, bottom: 8 },
    source: "实测：LibRaw 输出传感器全幅 8280×5520，尼康工坊导出为 8256×5504，\
             差值左右各 12、上下各 8；与独立实现 rawler 报的裁切区 @(12,8) 吻合",
    color_matrix: Some(ColorMatrixEntry {
        matrix: [
            [0.688_90, 0.328_61, -0.017_58],
            [0.017_33, 1.261_10, -0.278_49],
            [-0.014_57, -0.107_67, 1.122_18],
        ],
        output_space: "ProPhoto 线性",
        source: "配对最小二乘推导。参考照片 simple/DSC_4143.NEF：6,529,215 个样本，\
                 rms 0.000005；跨图复验 DSC_3307.NEF：6,348,068 个样本，\
                 平均绝对差 0.000007、最大 0.000025。\
                 样本已排除任一侧饱和点（那里解码层做过裁切）",
    }),
}];

/// 把机型字符串规范化以便匹配：转小写、把连续空白折成单个空格、去首尾空白。
fn normalize(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ").to_lowercase()
}

/// 按机型字符串查表。
///
/// 匹配规则：规范化后**完全相等**或互为包含。用包含是为了容忍 LibRaw 与 EXIF 在
/// 厂商标注上的细微差异（如 `NIKON Z 8` 与 `Nikon Z 8`）。
pub fn lookup(model: &str) -> Option<&'static Model> {
    let key = normalize(model);
    if key.is_empty() {
        return None;
    }
    MODELS
        .iter()
        .find(|m| {
            let k = normalize(m.match_key);
            key == k || key.contains(&k) || k.contains(&key)
        })
}

/// 裁剪结果。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CropDecision {
    /// 已按机型边距裁切。
    Cropped { margins: Margins, source: &'static str },
    /// 调用方选择不裁切。
    Disabled,
    /// 机型不在表中，无法确定边距。
    UnknownModel { full_w: usize, full_h: usize },
    /// 边距大于全幅尺寸，数据不自洽。
    Inconsistent { full_w: usize, full_h: usize, margins: Margins },
}

impl CropDecision {
    /// 供结果记录使用的文字说明。
    pub fn describe(&self) -> String {
        match self {
            CropDecision::Cropped { margins, source } => format!(
                "已裁切：左{} 上{} 右{} 下{}（机型知识：{source}）",
                margins.left, margins.top, margins.right, margins.bottom
            ),
            CropDecision::Disabled => "未裁切：按调用方设置保留传感器全幅".to_string(),
            CropDecision::UnknownModel { full_w, full_h } => format!(
                "未裁切：机型不在边距表中，无法确定有效像素区（当前 {full_w}×{full_h}）"
            ),
            CropDecision::Inconsistent { full_w, full_h, margins } => format!(
                "未裁切：边距（左{} 上{} 右{} 下{}）与全幅 {full_w}×{full_h} 不自洽",
                margins.left, margins.top, margins.right, margins.bottom
            ),
        }
    }
}

/// 决定如何裁切。
///
/// `enabled` 为 `false` 时一律保留全幅——这是"裁切可开关"这条规范的入口。
/// 无论走哪条分支，返回的决定都带有可读说明，使调用方能把"为什么是这个尺寸"
/// 一并报告出去，而不是只给一个数字。
pub fn decide(model: &str, full_w: usize, full_h: usize, enabled: bool) -> CropDecision {
    if !enabled {
        return CropDecision::Disabled;
    }
    let Some(m) = lookup(model) else {
        return CropDecision::UnknownModel { full_w, full_h };
    };
    if m.margins.effective(full_w, full_h).is_none() {
        return CropDecision::Inconsistent { full_w, full_h, margins: m.margins };
    }
    CropDecision::Cropped { margins: m.margins, source: m.source }
}

/// 按决定裁切一份行主序的三通道 16 位像素。
///
/// 返回裁切后的像素与尺寸；未裁切时原样返回。
pub fn apply(
    pixels: &[u16],
    width: usize,
    height: usize,
    decision: &CropDecision,
) -> (Vec<u16>, usize, usize) {
    let CropDecision::Cropped { margins, .. } = decision else {
        return (pixels.to_vec(), width, height);
    };
    let Some((nw, nh)) = margins.effective(width, height) else {
        return (pixels.to_vec(), width, height);
    };
    debug_assert_eq!(pixels.len(), width * height * 3);

    let mut out = Vec::with_capacity(nw * nh * 3);
    for y in 0..nh {
        let sy = y + margins.top;
        let start = (sy * width + margins.left) * 3;
        out.extend_from_slice(&pixels[start..start + nw * 3]);
    }
    (out, nw, nh)
}

/// 从文件读出相机型号。
///
/// 走 `src/tiff.rs` 读 IFD0 的 `Model`(0x0110) 与 `Make`(0x010F)，**不经过解码层**——
/// 解码层的 `libraw_get_iparams` 返回的是一个大结构体指针，复刻它的布局会在
/// 版本升级时静默读到错误偏移。型号是标准 TIFF 标签，自己读更稳。
///
/// 只读文件前部：型号位于 IFD0，不需要整份文件。
pub fn read_model(path: &std::path::Path) -> Option<String> {
    use crate::tiff::Tiff;
    use std::io::Read;

    const PREFIX: u64 = 256 * 1024;
    let mut f = std::fs::File::open(path).ok()?;
    let mut buf = Vec::new();
    f.by_ref().take(PREFIX).read_to_end(&mut buf).ok()?;

    let t = Tiff::locate(&buf).ok()?;
    let ifd0 = t.ifd0_offset().ok()?;
    let read = |tag: u16| -> Option<String> {
        let e = t.find(ifd0, tag).ok()??;
        let b = t.bytes(&e).ok()?;
        let end = b.iter().position(|&c| c == 0).unwrap_or(b.len());
        let s = String::from_utf8_lossy(&b[..end]).trim().to_string();
        (!s.is_empty()).then_some(s)
    };
    let model = read(0x0110)?;
    // 有些机型的 Model 已含厂商名，避免出现 "NIKON NIKON Z 8"
    match read(0x010F) {
        Some(make) if !model.to_lowercase().contains(&make.to_lowercase()) => {
            Some(format!("{make} {model}"))
        }
        _ => Some(model),
    }
}

/// 按机型取相机 → 工作空间的色彩矩阵。
///
/// 未标定的机型返回 `None`——**不得用别的机型的矩阵顶替**，也不得退回解码层那个
/// 实测有误的矩阵。
pub fn color_matrix(model: &str) -> Option<&'static ColorMatrixEntry> {
    lookup(model).and_then(|m| m.color_matrix.as_ref())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn z8_margins_are_lookupable() {
        for name in ["NIKON Z 8", "Nikon Z 8", "nikon  z  8"] {
            let m = lookup(name).unwrap_or_else(|| panic!("应能匹配 {name:?}"));
            assert_eq!(m.display, "Nikon Z 8");
            assert_eq!(m.margins.left, 12);
            assert_eq!(m.margins.top, 8);
        }
    }

    #[test]
    fn unknown_model_is_not_guessed() {
        assert!(lookup("Canon EOS R5").is_none());
        assert!(lookup("").is_none());
        let d = decide("Canon EOS R5", 8280, 5520, true);
        assert!(matches!(d, CropDecision::UnknownModel { .. }));
        assert!(d.describe().contains("机型不在边距表中"));
    }

    #[test]
    fn disabled_keeps_full_frame() {
        let d = decide("NIKON Z 8", 8280, 5520, false);
        assert_eq!(d, CropDecision::Disabled);
        let px = vec![7u16; 4 * 3 * 3];
        let (out, w, h) = apply(&px, 4, 3, &d);
        assert_eq!((w, h), (4, 3));
        assert_eq!(out, px);
    }

    #[test]
    fn z8_crop_gives_effective_area() {
        let d = decide("NIKON Z 8", 8280, 5520, true);
        let CropDecision::Cropped { margins, .. } = d else {
            panic!("应为已裁切，得到 {d:?}");
        };
        assert_eq!(margins.effective(8280, 5520), Some((8256, 5504)));

        // 造一张可辨识的图：像素值 = 行号*1000 + 列号
        let (w, h) = (8280usize, 5520usize);
        let mut px = vec![0u16; w * h * 3];
        for y in 0..h {
            for x in 0..w {
                let i = (y * w + x) * 3;
                let v = (y as u16).wrapping_mul(1000).wrapping_add(x as u16);
                px[i] = v;
                px[i + 1] = v;
                px[i + 2] = v;
            }
        }
        let (out, nw, nh) = apply(&px, w, h, &d);
        assert_eq!((nw, nh), (8256, 5504));
        assert_eq!(out.len(), nw * nh * 3);
        // 首像素应来自 (left, top) = (12, 8)
        let expect = (8u16).wrapping_mul(1000).wrapping_add(12);
        assert_eq!(out[0], expect, "裁切起点应为左上角边距处");
        // 末像素应来自 (right-1, bottom-1)
        let last_y = h - 8 - 1;
        let last_x = w - 12 - 1;
        let expect_last = (last_y as u16).wrapping_mul(1000).wrapping_add(last_x as u16);
        assert_eq!(out[out.len() - 3], expect_last, "裁切终点应为右下角边距前一处");
    }

    #[test]
    fn inconsistent_margins_are_reported() {
        // 用一张比边距还小的图触发不自洽分支
        let d = decide("NIKON Z 8", 10, 10, true);
        assert!(matches!(d, CropDecision::Inconsistent { .. }), "得到 {d:?}");
    }

    #[test]
    fn every_model_entry_has_a_source() {
        for m in MODELS {
            assert!(!m.source.trim().is_empty(), "{} 缺少边距来源依据", m.display);
            if let Some(cm) = &m.color_matrix {
                assert!(!cm.source.trim().is_empty(), "{} 缺少矩阵来源依据", m.display);
                assert!(!cm.output_space.trim().is_empty(), "{} 未标注矩阵的目的空间", m.display);
            }
        }
    }

    #[test]
    fn z8_color_matrix_maps_neutral_to_neutral() {
        let cm = color_matrix("NIKON Z 8").expect("Z8 应有色彩矩阵");
        for (i, s) in cm.row_sums().iter().enumerate() {
            assert!(
                (s - 1.0).abs() < 1e-3,
                "第 {i} 行行和为 {s}，中性应被映射为中性"
            );
        }
    }

    #[test]
    fn color_matrix_is_not_the_decoder_one() {
        // 回归：解码层 libraw_get_rgb_cam 给出的矩阵与输出空间无关，不是实际使用的那个。
        // 若有人图省事换成它，这里会失败。
        let cm = color_matrix("NIKON Z 8").unwrap();
        let decoder_first_row = [1.393_10f32, -0.215_69, -0.177_41];
        let d: f32 = (0..3).map(|j| (cm.matrix[0][j] - decoder_first_row[j]).abs()).sum();
        assert!(
            d > 0.5,
            "矩阵第 0 行与解码层给出的过于接近（差 {d}）——那一个实测不是实际使用的矩阵"
        );
    }

    #[test]
    fn uncalibrated_model_returns_no_matrix() {
        assert!(color_matrix("Canon EOS R5").is_none());
        assert!(color_matrix("").is_none());
    }
}
