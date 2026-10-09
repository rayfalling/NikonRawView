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
        }
    }
}
