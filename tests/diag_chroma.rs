//! 色度专项诊断：单张参考图的色度是否被正确还原，以及它本身够不够格当色度验证样本。
//!
//! # 为什么需要单独一条
//!
//! 团队诊断的结论是：暗部修好之后，**残余误差已转为彩度主导**（全体 ΔC* = +1.14，
//! 明度项只占 6%）。而要验证色度，参考图本身必须有足够的色度变化——原先那 15 张是
//! 42 秒内连拍的同一场景，占用前 5 的 LUT 格子全部贴着中性轴，色度那两维几乎没被
//! 用到。**用一组近乎中性的样本去验证色度，等于没验。**
//!
//! 所以这条诊断同时回答两件事：
//!
//! 1. **这张图的色度够不够格**——参考画面的 C*（Lab 彩度）分布有多宽，有多少比例
//!    落在"明显有色"的区间。若绝大多数像素 C* 很低，那么它测不出色度问题。
//! 2. **我们的色度还原得准不准**——ΔC* 与 ΔH' 的中位数、P95、以及它们在总误差量里
//!    的占比（用 `deltae::Decomposition` 的三口径）。
//!
//! 用法：设 `NRV_DIAG_STEM` 环境变量指定图名（默认 `DSC_0141`）。

mod common;

use common::{fitted_fixture_for, quantiles};
use nikonrawview::deltae;

/// 把 ProPhoto 线性的参考值转成 Lab，供色度统计。
fn lab(v: [f32; 3]) -> [f32; 3] {
    deltae::lab_from_prophoto_linear(v)
}

/// Lab 彩度。
fn chroma(l: [f32; 3]) -> f64 {
    ((l[1] as f64).powi(2) + (l[2] as f64).powi(2)).sqrt()
}

fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[v.len() / 2]
}

#[test]
fn chroma_accuracy_and_sample_adequacy() {
    let stem = std::env::var("NRV_DIAG_STEM").unwrap_or_else(|_| "DSC_0141".into());
    let Some(fx) = fitted_fixture_for(&stem, true) else {
        eprintln!("跳过：缺少 {stem} 的样本对");
        return;
    };

    eprintln!("=== {stem}：色度专项诊断 ===");
    eprintln!("拟合样本 {}，留出评估样本 {}", fx.fit_count, fx.eval.len());
    eprintln!("尺寸 {}×{}", fx.size.0, fx.size.1);

    // ---------- 一、这张图够不够格当色度验证样本 ----------
    let ref_c: Vec<f64> = fx.eval.iter().map(|s| chroma(lab(s.theirs_linear))).collect();
    let ours_c: Vec<f64> = fx.eval.iter().map(|s| chroma(lab(s.ours_linear))).collect();
    let ref_med = median(&mut ref_c.clone());
    let ours_med = median(&mut ours_c.clone());
    let colored = ref_c.iter().filter(|c| **c > 10.0).count();
    let rich = ref_c.iter().filter(|c| **c > 25.0).count();
    let n = ref_c.len().max(1);

    eprintln!("\n【参考画面的色度内容】");
    eprintln!("  参考 C* 中位 {ref_med:.2}   （我们的原始解码 C* 中位 {ours_med:.2}）");
    eprintln!("  参考 C* > 10（明显有色）: {colored} 个，占 {:.2}%", 100.0 * colored as f64 / n as f64);
    eprintln!("  参考 C* > 25（高彩度）  : {rich} 个，占 {:.2}%", 100.0 * rich as f64 / n as f64);

    let adequate = (colored as f64 / n as f64) > 0.10;
    eprintln!(
        "  判定：{}",
        if adequate {
            "**够格**——有色像素超过一成，足以测出色度问题"
        } else {
            "**不够格**——有色像素不足一成，用它验色度等于没验"
        }
    );

    // ---------- 二、色度还原得准不准 ----------
    let mut all: Vec<f64> = fx.eval.iter().map(|s| s.delta_e()).collect();
    let (m, p, x) = quantiles(&mut all);
    eprintln!("\n【ΔE00】中位 {m:.3}  P95 {p:.3}  最大 {x:.3}  （目标 ≤1.0 / ≤3.0）");

    // 分量拆解：用 decomp 的三口径
    let lab_pairs: Vec<([f32; 3], [f32; 3])> = fx
        .eval
        .iter()
        .map(|s| (lab(s.ours_transformed), lab(s.theirs_linear)))
        .collect();
    let sum = deltae::summarize(&lab_pairs);
    eprintln!("\n【误差量归因（{stem}）】");
    for (name, share) in deltae::COMPONENTS.iter().zip(sum.total_share.iter()) {
        eprintln!("  {name:<10} 合计占比 {:>6.1}%", 100.0 * share);
    }
    eprintln!("  —— 主因票数：");
    for (name, frac) in deltae::COMPONENTS.iter().zip(sum.primary_fraction.iter()) {
        eprintln!("    {name:<10} {:>6.1}%", 100.0 * frac);
    }

    // ΔC* 与 ΔH' 的实际分布（不只占比）
    let mut dc: Vec<f64> = Vec::with_capacity(fx.eval.len());
    let mut dl: Vec<f64> = Vec::with_capacity(fx.eval.len());
    for s in &fx.eval {
        let d = deltae::decompose(lab(s.ours_transformed), lab(s.theirs_linear));
        dc.push(d.dc.abs());
        dl.push(d.dl.abs());
    }
    let dc_med = median(&mut dc.clone());
    let dl_med = median(&mut dl.clone());
    eprintln!("\n【分量的绝对量级（中位）】");
    eprintln!("  |ΔL'/SL| = {dl_med:.3}    |ΔC'/SC| = {dc_med:.3}");

    // 色度只在高彩度区才有意义——单独看有色像素的色度误差
    let colored_dc: Vec<f64> = fx
        .eval
        .iter()
        .filter(|s| chroma(lab(s.theirs_linear)) > 10.0)
        .map(|s| {
            deltae::decompose(lab(s.ours_transformed), lab(s.theirs_linear)).dc.abs()
        })
        .collect();
    if !colored_dc.is_empty() {
        let cm = median(&mut colored_dc.clone());
        let mut sorted = colored_dc.clone();
        let (_, cp95, _) = quantiles(&mut sorted);
        eprintln!("  仅看有色像素（C*>10，{} 个）：|ΔC'/SC| 中位 {cm:.3}，P95 {cp95:.3}", colored_dc.len());
    }

    eprintln!("\n判读：色度若要符合验证要求，需要同时满足两点——");
    eprintln!("  1. 样本本身有足够的有色像素（否则测不出色度）；");
    eprintln!("  2. 有色像素上的 |ΔC'/SC| 与 ΔH' 处于 JND 以下（ΔE00 中位 ≤1.0）。");
}
