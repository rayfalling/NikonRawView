//! 端到端复测：**接上 `w(R/G, B/G)` 之后，ΔE00 到底降到多少。**
//!
//! # 这是第一次让整条链真正合起来
//!
//! 之前所有数字都是"色度尺"上的中间量（尺度不变的 `Δlog2` 比值差）。而最终判据是
//! **端到端 ΔE00**（目标：中位 ≤1.0、P95 ≤3.0）。
//!
//! # 为什么这里用相似变换
//!
//! `w` 的实测位置是**矩阵之前**：`p = M0·diag(w)·c`。而脚手架的 `DiagSample::ours_linear`
//! 是**矩阵之后**的值（LibRaw 已经施加过矩阵，且已被证明与 `M0` 一致到 5e-6）。
//!
//! 于是把它搬到可用的那一级：
//!
//! ```text
//! ours_linear = M0·c                    （脚手架给的就是这个）
//! c           = M0⁻¹·ours_linear
//! p           = M0·diag(w)·c = (M0·diag(w)·M0⁻¹)·ours_linear
//! ```
//!
//! 中间那个 `M0·diag(w)·M0⁻¹` 是相似变换，**一般不是对角的**——所以不能简单地把
//! `w` 乘到 `ours_linear` 上。这一步搞错会让结果完全不对，故在此写明。

mod common;

use common::{fitted_fixture_for, quantiles};
use nikonrawview::deltae;
use nikonrawview::mat3::{self, Mat3};

/// 把"矩阵之前"的逐通道增益搬到"矩阵之后"的等价形式。
fn to_post_matrix(gain: [f64; 3], m0: Mat3) -> Option<Mat3> {
    let inv = mat3::inverse(m0)?;
    let d: Mat3 = [
        [gain[0] as f32, 0.0, 0.0],
        [0.0, gain[1] as f32, 0.0],
        [0.0, 0.0, gain[2] as f32],
    ];
    Some(mat3::mul(m0, mat3::mul(d, inv)))
}

/// # 跨图测试：**一套固定变换套到白平衡不同的图上**
///
/// 上面那个测试是**每张图自己拟合一套曲线+LUT**——那种用法下逐图 LUT 早就把白平衡
/// 差异吸收掉了，`w` 自然冗余（实测中位 3.529 → 3.539，几乎没动）。
///
/// 但本工程的目标是另一种用法：**从参考图标定出**一套**基准变换，套用到任意照片**。
/// 在那个场景下，一张白平衡与标定集不同的照片会被系统性偏色——**`w` 正是修这个的**。
///
/// 所以这里做留一图：用**其余所有图**的样本池化拟合一整套（曲线+LUT），固定不动，
/// 再套到留出的那张上，比较接 `w` 与不接 `w`。
///
/// 这才是 `w` 该被评判的地方。
#[test]
fn cross_image_transfer_with_and_without_wb_adaptation() {
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

    // 先把所有样本取出来（每张图解码一次），并按固定步长抽稀以免拟合过慢
    struct Info {
        stem: String,
        rg: f64,
        bg: f64,
        sim: Mat3,
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
        let rg = (wb[0] / wb[1]) as f64;
        let bg = (wb[2] / wb[1]) as f64;
        let gain = nikonrawview::color::wb_adaptation(rg, bg);
        let Some(sim) = to_post_matrix(gain, m0) else { continue };
        // 抽稀：每 7 个取 1 个，保留足够的分布形状
        let samples: Vec<nikonrawview::fit::Sample> = fx
            .eval
            .iter()
            .step_by(7)
            .map(|s| nikonrawview::fit::Sample { ours: s.ours_linear, theirs: s.theirs_linear })
            .collect();
        if samples.len() < 2000 {
            continue;
        }
        all.push(Info { stem: stem.to_string(), rg, bg, sim, samples });
    }
    if all.len() < 4 {
        eprintln!("跳过：可用图不足（{} 张）", all.len());
        return;
    }

    eprintln!("=== 跨图测试：一套固定变换套到白平衡不同的留出图上 ===");
    eprintln!("标定集 = 除留出图之外的全部；每张图轮流留出\n");
    eprintln!(
        "  {:<11} {:>7} {:>7} {:>9} {:>9} {:>10}",
        "留出图", "R/G", "B/G", "不接 w", "接 w", "改善"
    );

    let mut plain_all = Vec::new();
    let mut adap_all = Vec::new();

    for i in 0..all.len() {
        // 池化其余图的样本作为标定集
        let mut pooled: Vec<nikonrawview::fit::Sample> = Vec::new();
        for (j, info) in all.iter().enumerate() {
            if j != i {
                pooled.extend_from_slice(&info.samples);
            }
        }
        let curve = nikonrawview::fit::fit_curve(&pooled, 128);
        let edge = 17;
        let lut = nikonrawview::fit::fit_lut(&pooled, &curve, edge);
        let t = nikonrawview::transform::BaseTransform::from_fit(
            0,
            "cross",
            vec![all[i].stem.clone()],
            curve,
            edge,
            lut,
        );

        let held = &all[i];
        let mut d0: Vec<f64> = Vec::with_capacity(held.samples.len());
        let mut d1: Vec<f64> = Vec::with_capacity(held.samples.len());
        for s in &held.samples {
            d0.push(deltae::delta_e_prophoto(t.apply(s.ours), s.theirs));
            let adapted = mat3::mul_vec(held.sim, s.ours);
            d1.push(deltae::delta_e_prophoto(t.apply(adapted), s.theirs));
        }
        let (m0_, _, _) = quantiles(&mut d0);
        let (m1_, _, _) = quantiles(&mut d1);
        eprintln!(
            "  {:<11} {:>7.4} {:>7.4} {:>9.3} {:>9.3} {:>9.1}%",
            held.stem, held.rg, held.bg, m0_, m1_, 100.0 * (m0_ - m1_) / m0_.max(1e-9)
        );
        plain_all.push(m0_);
        adap_all.push(m1_);
    }

    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    let mut sp = plain_all.clone();
    let mut sa = adap_all.clone();
    sp.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    sa.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    eprintln!("\n=== 汇总（{} 张留出图）===", plain_all.len());
    eprintln!("  {:<16} {:>10} {:>10}", "", "不接 w", "接 w");
    eprintln!("  {:<16} {:>10.3} {:>10.3}", "ΔE00 中位均值", mean(&plain_all), mean(&adap_all));
    eprintln!("  {:<16} {:>10.3} {:>10.3}", "ΔE00 中位的中位", sp[sp.len() / 2], sa[sa.len() / 2]);
    let improved = adap_all.iter().zip(&plain_all).filter(|(a, b)| **a < **b).count();
    eprintln!("  {:<16} {:>10}", "改善张数", format!("{improved}/{}", plain_all.len()));

    eprintln!("\n判读：");
    let gain = (mean(&plain_all) - mean(&adap_all)) / mean(&plain_all).max(1e-9);
    if gain > 0.10 {
        eprintln!("  跨图场景下 `w` 有实质贡献（均值改善 {:.1}%）——**这就是它的用武之地**。", gain * 100.0);
    } else if gain > 0.02 {
        eprintln!("  跨图场景下 `w` 贡献有限（{:.1}%）。", gain * 100.0);
    } else {
        eprintln!("  **跨图场景下 `w` 也没有贡献（{:.1}%）** —— 那它就不可迁移，", gain * 100.0);
        eprintln!("  留一状态 CV 的好数字可能只反映了'同一组数据互相预测'。");
    }
    eprintln!("\n注意：这里的标定集与留出图**共享场景类别**（同一批参考图），");
    eprintln!("所以仍然测不出'换场景还成不成立'——那需要用户再提供别的场景的样张。");
}

/// 用给定的前端配对拟合并评估完整后端（曲线 128 箱 + 17³ LUT），返回 ΔE00 分位数。
fn fit_eval(pairs: &[([f32; 3], [f32; 3])]) -> (f64, f64, f64) {
    let samples: Vec<nikonrawview::fit::Sample> = pairs
        .iter()
        .map(|(c, t)| nikonrawview::fit::Sample { ours: *c, theirs: *t })
        .collect();
    let curve = nikonrawview::fit::fit_curve(&samples, 128);
    let edge = 17;
    let lut = nikonrawview::fit::fit_lut(&samples, &curve, edge);
    let t = nikonrawview::transform::BaseTransform::from_fit(0, "e2e", vec![], curve, edge, lut);
    let mut des: Vec<f64> = samples
        .iter()
        .map(|s| deltae::delta_e_prophoto(t.apply(s.ours), s.theirs))
        .collect();
    quantiles(&mut des)
}

#[test]
fn end_to_end_delta_e_with_and_without_wb_adaptation() {
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

    eprintln!("=== 端到端 ΔE00：接上 w(R/G, B/G) 前后 ===");
    eprintln!("目标：中位 ≤ 1.0，P95 ≤ 3.0\n");
    eprintln!(
        "  {:<11} {:>8} {:>8} {:>9}   {:>8} {:>8} {:>9}   {:>9}",
        "图", "R/G", "B/G", "w_R w_B", "基线中位", "基线P95", "接入中位", "改善"
    );

    let mut base = Vec::new();
    let mut adap = Vec::new();
    let mut rows = 0;

    for stem in stems {
        let Some(fx) = fitted_fixture_for(stem, false) else {
            eprintln!("  {stem:<11} （样本不可用，跳过）");
            continue;
        };
        if fx.eval.is_empty() {
            continue;
        }
        // cam_mul 的两个比值 —— 与渲染时用的是同一个来源
        let Some(wb) = nikonrawview::libraw::read_wb(
            &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("simple").join(format!("{stem}.NEF")),
        )
        .ok() else {
            continue;
        };
        let rg = (wb[0] / wb[1]) as f64;
        let bg = (wb[2] / wb[1]) as f64;
        let gain = nikonrawview::color::wb_adaptation(rg, bg);
        let Some(sim) = to_post_matrix(gain, m0) else {
            eprintln!("  {stem:<11} （矩阵不可逆，跳过）");
            continue;
        };

        let plain: Vec<([f32; 3], [f32; 3])> =
            fx.eval.iter().map(|s| (s.ours_linear, s.theirs_linear)).collect();
        let adapted: Vec<([f32; 3], [f32; 3])> = fx
            .eval
            .iter()
            .map(|s| (mat3::mul_vec(sim, s.ours_linear), s.theirs_linear))
            .collect();

        let (m0_, p0_, _) = fit_eval(&plain);
        let (m1_, p1_, _) = fit_eval(&adapted);

        eprintln!(
            "  {:<11} {:>8.4} {:>8.4} {:>4.3} {:>4.3}   {:>8.3} {:>8.3} {:>9.3}   {:>8.1}%",
            stem, rg, bg, gain[0], gain[2], m0_, p0_, m1_, 100.0 * (m0_ - m1_) / m0_.max(1e-9)
        );
        base.push((m0_, p0_));
        adap.push((m1_, p1_));
        rows += 1;
    }

    if rows == 0 {
        eprintln!("\n跳过：没有可用样本");
        return;
    }

    let med = |v: &[(f64, f64)]| {
        let mut x: Vec<f64> = v.iter().map(|p| p.0).collect();
        x.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        x[x.len() / 2]
    };
    let worst = |v: &[(f64, f64)]| {
        v.iter().map(|p| p.0).fold(f64::MIN, f64::max)
    };

    eprintln!("\n=== 汇总（{rows} 张）===");
    eprintln!("  {:<16} {:>10} {:>10}", "", "基线", "接入 w");
    eprintln!("  {:<16} {:>10.3} {:>10.3}", "ΔE00 中位的中位", med(&base), med(&adap));
    eprintln!("  {:<16} {:>10.3} {:>10.3}", "最差图的中位", worst(&base), worst(&adap));

    eprintln!("\n判读：");
    let improved = adap.iter().zip(&base).filter(|(a, b)| a.0 < b.0).count();
    eprintln!("  {rows} 张里 {improved} 张改善。");
    if med(&adap) <= 1.0 {
        eprintln!("  **中位已达标（≤1.0）**。");
    } else {
        eprintln!(
            "  **中位仍未达标**（{:.3} > 1.0）——w 只解决白平衡适配那一层，",
            med(&adap)
        );
        eprintln!("  剩下的要靠色调映射或其它环节。");
    }
}
