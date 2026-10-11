//! 选项 A：**逐通道曲线是否真的能在跨图协议下兑现。**
//!
//! # 为什么先做这个，不改格式
//!
//! 诊断 N 给出：跨图迁移惩罚（transfer penalty）的 **96.6% 由「输出端逐通道色调/尺度差」
//! 解释**，而共享曲线只有 18.6%。诊断 D（chroma-scale）独立指向同一件事。
//!
//! 但那两级都是**神谕**——用留出图自己的参考值拟合，量的是 penalty 的**形状**，
//! **不是可达增益**。而且本项目刚有过一次「中间量很好看、端到端没兑现」的先例：
//! `w` 在自拟合协议下吃掉 95% 的色度改善量，端到端却几乎为零。
//!
//! **所以先拿端到端数字，再决定是否改 `BaseTransform` 的格式。**
//!
//! # 协议
//!
//! 用**留一图**（协议 B，唯一有效的那种）：
//!
//! - 标定集 = 其余所有图的样本池化
//! - 留出图 = 固定不动地套上去
//! - 两个变体：**(i) 1 条共享曲线 + 17³ LUT**、**(ii) 3 条逐通道曲线 + 17³ LUT**
//! - 两个变体都接上 `w(R/G, B/G)`（它已被跨图协议验证有效）
//!
//! 逐通道变体的施加顺序：`ours --逐通道曲线--> x --17³ LUT--> 结果`。
//! 这样两边都有 LUT，差别只在曲线是一条还是三条。

mod common;

use common::{fitted_fixture_for, quantiles};
use nikonrawview::deltae;
use nikonrawview::mat3::{self, Mat3};

/// 把"矩阵之前"的逐通道增益搬到"矩阵之后"的等价形式（相似变换）。
fn to_post_matrix(gain: [f64; 3], m0: Mat3) -> Option<Mat3> {
    let inv = mat3::inverse(m0)?;
    let d: Mat3 = [
        [gain[0] as f32, 0.0, 0.0],
        [0.0, gain[1] as f32, 0.0],
        [0.0, 0.0, gain[2] as f32],
    ];
    Some(mat3::mul(m0, mat3::mul(d, inv)))
}

/// 恒等曲线——用来让 LUT 单独工作（逐通道变体里曲线已经单独施加过了）。
fn identity_curve(bins: usize) -> Vec<f32> {
    (0..=bins).map(|i| i as f32 / bins as f32).collect()
}

fn delta_med(pairs: &[([f32; 3], [f32; 3])]) -> f64 {
    let mut d: Vec<f64> = pairs
        .iter()
        .map(|(a, b)| deltae::delta_e_prophoto(*a, *b))
        .collect();
    quantiles(&mut d).0
}

#[test]
fn per_channel_curves_versus_shared_curve_on_held_out_images() {
    let m0 = match nikonrawview::camera::color_matrix("Nikon Z 8") {
        Some(e) => e.matrix,
        None => {
            eprintln!("跳过：机型表里没有 Z 8 的矩阵");
            return;
        }
    };

    let stems = [
        "DSC_0001", "DSC_0010", "DSC_0141", "DSC_0314", "DSC_0507", "DSC_0540", "DSC_0567",
        "DSC_0570", "DSC_8562",
    ];

    struct Info {
        stem: String,
        samples: Vec<nikonrawview::fit::Sample>,
    }
    let mut all: Vec<Info> = Vec::new();
    for stem in stems {
        let Some(fx) = fitted_fixture_for(stem, false) else {
            eprintln!("  跳过 {stem}：样本不可用");
            continue;
        };
        let nef = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("simple")
            .join(format!("{stem}.NEF"));
        let Ok(wb) = nikonrawview::libraw::read_wb(&nef) else { continue };
        let gain = nikonrawview::color::wb_adaptation(
            (wb[0] / wb[1]) as f64,
            (wb[2] / wb[1]) as f64,
        );
        let Some(sim) = to_post_matrix(gain, m0) else { continue };
        // 抽稀，并把 w 的相似变换预先施加到 ours 上（两边变体都一样）
        let samples: Vec<nikonrawview::fit::Sample> = fx
            .eval
            .iter()
            .step_by(7)
            .map(|s| nikonrawview::fit::Sample {
                ours: mat3::mul_vec(sim, s.ours_linear),
                theirs: s.theirs_linear,
            })
            .collect();
        if samples.len() < 2000 {
            continue;
        }
        all.push(Info { stem: stem.to_string(), samples });
    }
    if all.len() < 4 {
        eprintln!("跳过：可用图不足（{} 张）", all.len());
        return;
    }

    let edge = 17;
    let bins = 128;

    eprintln!("=== 留一图：共享曲线 vs 逐通道曲线（两边都接 w + 17³ LUT）===");
    eprintln!("目标：中位 ≤ 1.0，P95 ≤ 3.0\n");
    eprintln!(
        "  {:<11} {:>10} {:>12} {:>10}",
        "留出图", "共享曲线", "逐通道曲线", "改善"
    );

    let mut shared_all = Vec::new();
    let mut per_all = Vec::new();

    for i in 0..all.len() {
        let mut pooled: Vec<nikonrawview::fit::Sample> = Vec::new();
        for (j, info) in all.iter().enumerate() {
            if j != i {
                pooled.extend_from_slice(&info.samples);
            }
        }

        // (i) 共享曲线 + LUT
        let shared_curve = nikonrawview::fit::fit_curve(&pooled, bins);
        let shared_lut = nikonrawview::fit::fit_lut(&pooled, &shared_curve, edge);
        let t_shared = nikonrawview::transform::BaseTransform::from_fit(
            0,
            "shared",
            vec![],
            shared_curve,
            edge,
            shared_lut,
        );

        // (ii) 3 条逐通道曲线 + LUT（先把逐通道曲线施加到样本上，再拟 LUT）
        let per_curves = nikonrawview::fit::fit_curves_per_channel(&pooled, bins);
        let pre: Vec<nikonrawview::fit::Sample> = pooled
            .iter()
            .map(|s| nikonrawview::fit::Sample {
                ours: nikonrawview::fit::forward_curves(&per_curves, s.ours),
                theirs: s.theirs,
            })
            .collect();
        let idc = identity_curve(bins);
        let per_lut = nikonrawview::fit::fit_lut(&pre, &idc, edge);

        let held = &all[i];
        let sp: Vec<([f32; 3], [f32; 3])> =
            held.samples.iter().map(|s| (t_shared.apply(s.ours), s.theirs)).collect();
        let pp: Vec<([f32; 3], [f32; 3])> = held
            .samples
            .iter()
            .map(|s| {
                // 逐通道曲线 → 恒等曲线（不做事）→ 查 LUT
                let x = nikonrawview::fit::forward_curves(&per_curves, s.ours);
                (nikonrawview::fit::apply_lut(&per_lut, edge, x), s.theirs)
            })
            .collect();

        let ms = delta_med(&sp);
        let mp = delta_med(&pp);
        eprintln!(
            "  {:<11} {:>10.3} {:>12.3} {:>9.1}%",
            held.stem,
            ms,
            mp,
            100.0 * (ms - mp) / ms.max(1e-9)
        );
        shared_all.push(ms);
        per_all.push(mp);
    }

    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    let mut ss = shared_all.clone();
    let mut pps = per_all.clone();
    ss.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    pps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    eprintln!("\n=== 汇总（{} 张留出图）===", shared_all.len());
    eprintln!("  {:<18} {:>10} {:>10}", "", "共享曲线", "逐通道曲线");
    eprintln!("  {:<18} {:>10.3} {:>10.3}", "ΔE00 中位均值", mean(&shared_all), mean(&per_all));
    eprintln!("  {:<18} {:>10.3} {:>10.3}", "ΔE00 中位的中位", ss[ss.len() / 2], pps[pps.len() / 2]);
    let improved = per_all.iter().zip(&shared_all).filter(|(a, b)| **a < **b).count();
    eprintln!("  {:<18} {:>10}", "改善张数", format!("{improved}/{}", shared_all.len()));

    let gain = (mean(&shared_all) - mean(&per_all)) / mean(&shared_all).max(1e-9);
    eprintln!("\n判读：");
    eprintln!("  跨图场景下逐通道曲线改善 **{:.1}%**。", gain * 100.0);
    if gain > 0.20 {
        eprintln!("  ⇒ **确实兑现**，值得改 BaseTransform 的格式（曲线从 1 条变 3 条）。");
    } else if gain > 0.05 {
        eprintln!("  ⇒ 有改善但不显著。**在改格式之前应当先问为什么兑现不了**——");
        eprintln!("    诊断 N 的 L3 只是个上界（神谕），这里才是可达值。");
    } else {
        eprintln!("  ⇒ **几乎没有兑现**。那 'L3 吃掉 96.6%' 就不是可达增益，");
        eprintln!("    不要据此改格式；应回到'缺的是别的环节'这条线。");
    }
    eprintln!("\n限制：9 张图、两类场景，留一图仍共享场景类别，**跨场景这一维依然未验证**。");
}
