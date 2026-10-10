//! 暗场判据：**我们扣黑电平之后，纯黑画面是否还是纯黑？**
//!
//! # 为什么这条能定案
//!
//! 用户提供了一张**盖上镜头盖拍摄的纯黑画面**（`DSC_0561`），而且它有配套的尼康
//! 导出 TIF。于是我们有了一个直接可比的对象：
//!
//! ```text
//! 我们的解码（相机空间 → ProPhoto 线性）   vs   尼康导出的同一张
//! ```
//!
//! 两边都是纯黑输入，所以任何**非零的正值都是 pedestal**。
//!
//! # 背景（诊断 L 的源码结论）
//!
//! LibRaw 从 MakerNote tag `0x003d` 读到黑电平 `1008`（我们用自己的解析器复现了
//! 同一个值），**逐通道扣、只扣一次、整数、扣完钳到 0**，没有精度损失。
//!
//! 但源码看不出 **1008 本身准不准**。方向性是可判的：
//!
//! - `1008` 偏低 → 扣完仍有正残差 → **暗部被抬高**
//! - `1008` 偏高 → 扣完被钳到 0 → 暗部被压死（不是抬高）
//!
//! 而我们已知的症状是"近黑被抬到约 0.017，参考是 0.0045"——**正是抬高那个方向**。
//! 本测试就把这件事量出来。

mod common;

use common::samples_dir;

/// 统计结果：均值与全零占比——判读只需要这两个。
struct Stats {
    mean: f64,
    zero_frac: f64,
}

fn stats(px: &[[f32; 3]], label: &str) -> Option<Stats> {
    if px.is_empty() {
        return None;
    }
    // 用亮度（Rec.709）表征"有多黑"
    let mut lum: Vec<f64> = px
        .iter()
        .map(|p| 0.2126 * p[0] as f64 + 0.7152 * p[1] as f64 + 0.0722 * p[2] as f64)
        .collect();
    let n = lum.len();
    let mean = lum.iter().sum::<f64>() / n as f64;
    let zeros = lum.iter().filter(|v| **v <= 1e-7).count();
    lum.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // 最近秩分位数，与本项目其它报告一致
    let q = |f: f64| lum[((f * n as f64).ceil() as usize).clamp(1, n) - 1];
    eprintln!(
        "  {label:<28} n={n:<9} 均值={mean:.6} 中位={:.6} P99={:.6} 最大={:.6} 全零占比={:.2}%",
        q(0.5),
        q(0.99),
        lum[n - 1],
        100.0 * zeros as f64 / n as f64
    );
    Some(Stats { mean, zero_frac: zeros as f64 / n as f64 })
}

#[test]
fn dark_frame_reveals_whether_we_leave_a_pedestal() {
    let dir = samples_dir();
    let stem = "DSC_0561";
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        eprintln!("跳过：需要 {stem}.NEF 与 {stem}.TIF（暗场及其导出）");
        return;
    }

    // 黑电平（用我们自己的解析器，与 LibRaw 读的是同一个标签）
    match nikonrawview::makernote::black_level(&std::fs::read(&nef).unwrap()) {
        Ok(Some(bl)) => eprintln!("LibRaw 用的黑电平：{}\n", bl.describe()),
        _ => eprintln!("黑电平：读不到\n"),
    }
    if let Ok(Some(k)) = nikonrawview::makernote::color_temperature(&std::fs::read(&nef).unwrap()) {
        eprintln!("色温：{k} K\n");
    }

    eprintln!("=== 暗场对比（纯黑输入，任何正值都是 pedestal）===");

    let fx = common::fitted_fixture_for(stem, false);
    let Some(fx) = fx else {
        eprintln!("跳过：脚手架取不到样本（尺寸或配对不符）");
        return;
    };

    let ours: Vec<[f32; 3]> = fx.eval.iter().map(|s| s.ours_linear).collect();
    let theirs: Vec<[f32; 3]> = fx.eval.iter().map(|s| s.theirs_linear).collect();

    let a = stats(&ours, "我们的解码（相机空间）");
    let b = stats(&theirs, "尼康导出（参考）");

    let (Some(a), Some(b)) = (a, b) else {
        eprintln!("跳过：样本为空");
        return;
    };

    eprintln!("\n=== 判读 ===");
    let ratio = a.mean / b.mean.max(1e-12);
    eprintln!("  我们的均值 / 参考均值 = {ratio:.2}×");
    if b.mean <= 1e-6 && a.mean <= 1e-6 {
        eprintln!("  两边都是纯零：**黑电平正确，没有 pedestal**。");
    } else if a.mean > 3.0 * b.mean.max(1e-9) {
        eprintln!("  **我们的暗场明显比参考亮** —— 扣黑之后仍有正残差，");
        eprintln!("  即 1008 低于真实黑电平。这是暗部被抬高的直接证据。");
        eprintln!(
            "  残差量级：我们的均值 {:.6}（参考 {:.6}），全零占比 {:.2}%（参考 {:.2}%）",
            a.mean,
            b.mean,
            100.0 * a.zero_frac,
            100.0 * b.zero_frac
        );
    } else if a.mean < b.mean / 3.0 {
        eprintln!("  我们比参考更黑 —— 可能是我们扣多了，或参考端抬了黑。");
    } else {
        eprintln!("  两者量级相当：黑电平这条线**不能解释暗部偏亮**。");
    }

    eprintln!("\n附注：");
    eprintln!("  · 这张暗场的曝光参数与其他参考图未必相同，所以量级只能与**同图的**");
    eprintln!("    尼康导出相比，不能跨图比较。");
    eprintln!("  · 若两者都接近零但画面仍整体偏亮，说明问题不在暗场，而在色调映射。");
}
