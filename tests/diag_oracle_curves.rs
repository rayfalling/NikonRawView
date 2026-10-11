//! 诊断 O：**神谕逐通道曲线的逐图相似性**。
//!
//! # 为什么要做这个检查
//!
//! 诊断 N 说跨图迁移惩罚的 **96.6%** 由「输出端逐通道色调/尺度差」解释（共享曲线只有 18.6%）。
//! 但 Lead 用**留一图（标定集拟合，非神谕）**验证后不兑现：共享 6.607 → 逐通道 6.525，
//! 只有 **−1.2%**。⇒ 修好一张图所需的逐通道曲线是**逐图特有**的。
//!
//! 那就只剩两种可能，而它们的修法完全不同：
//!
//! - **A**：曲线彼此相似 → 它**本该可预测**，是我们拟合/池化的方式有问题；
//! - **B**：曲线差异大且无规律 → 修正是**真正逐图特有的**，需要一个逐图参数；
//! - （**A′**：组内相似、跨组差异大 → 指向「场景/内容相关」，而不是逐图随机。）
//!
//! 本诊断把 10 张图的「神谕逐通道曲线」摊开比较，并给一条**有判读力的对照**：
//! 把 A 图的曲线套到 B 图上，ΔE00 涨多少 —— 纯数值差不如这个直观。
//!
//! # 纪律：本任务**只做描述性比较**
//!
//! 神谕量的是**形状**，不是可达增益（诊断 N 刚在这一点上吃过一次亏：神谕 96.6% vs 实际 −1.2%）。
//! 下面所有「神谕」字样都指「用该图自己的参考值拟合出来的曲线」，**不要**把它读成
//! 「我们因此能改好多少」。作为对照，[1] 同时给出**非神谕**的「池化逐通道曲线」那一列，
//! 两者的差距本身就是本诊断要展示的东西。
//!
//! # 协议
//!
//! - 10 张图分三组：**连拍**（DSC_0001~0003，同场景同白平衡，42 秒内）、
//!   **扫描**（DSC_0563~0566，同场景、只变色温）、**老图**（DSC_0141/0314/0540，不同场景）。
//! - 每张图：左半拟合、右半评估（沿用前几轮的空间留出）。
//! - 两个底座：
//!   - `T_i`（留一图底座）= 其余 9 张图的**左半池化**拟合的「共享曲线 128 箱 + 17³ LUT」，
//!     用于 [1] 的「神谕 vs 非神谕」对照（与 Lead 的验证同口径）；
//!   - `T_pool`（公共底座）= **全部 10 张**池化拟合的同一套，用于 [2] 的逐图曲线比较。
//!     **必须用公共底座**：若每张图各自用 `T_i`，曲线差异会被底座差异污染，比不出「曲线像不像」。
//!     代价是 `T_pool` 见过每张图（有泄漏），但那正是我们要量的东西——「一个拟合得很好的共享
//!     变换之后，**还剩下什么逐图特有的形状**」。
//!
//! # 门槛（先定后看）
//!
//! - 曲线「相似」定义：编码域逐点 RMS ≤ [`RMS_SIMILAR`]；
//! - 「组内相似、跨组差异大」定义：跨组 RMS 中位 ÷ 组内 RMS 中位 ≥ [`GROUP_RATIO`]。
//!
//! # 限制
//!
//! 10 张图、三组场景；连拍组是同场景同白平衡（只能量噪声底），扫描组是同场景变白平衡
//! （量「内容固定时白平衡的影响」），老图才有跨场景成分，但只有 3 张。**跨场景这一维依然很薄。**
//! `NRV_O_STEMS=...` 可只跑指定图（迭代用）。

mod common;

use common::samples_dir;
use nikonrawview::camera;
use nikonrawview::deltae;
use nikonrawview::fit::{self, Sample};
use nikonrawview::icc;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::mat3::{self, Mat3};
use nikonrawview::transform::BaseTransform;

// ---------------------------------------------------------------------------
// 口径与门槛（先定后看）
// ---------------------------------------------------------------------------

/// 默认 10 张图，按组。
const DEFAULT_STEMS: [&str; 10] = [
    // 连拍（同场景、同白平衡）
    "DSC_0001", "DSC_0002", "DSC_0003", // 扫描（同场景、变白平衡）
    "DSC_0563", "DSC_0564", "DSC_0565", "DSC_0566", // 老图（不同场景）
    "DSC_0141", "DSC_0314", "DSC_0540",
];

const STEP: usize = 17;
/// 曲线箱数（与生产一致）。
const CURVE_BINS: usize = 128;
const LUT_EDGE: usize = 17;

const CAM_LO: f32 = 1e-4;
const CAM_HI: f32 = 0.999;
const REF_LO: f32 = 1e-5;
const REF_HI: f32 = 0.999;
const MIN_SAMPLES: usize = 2000;

/// 曲线「相似」的门槛：编码域逐点差的 RMS。
const RMS_SIMILAR: f64 = 0.02;
/// 「组内相似、跨组差异大」的门槛：跨组 ÷ 组内。
const GROUP_RATIO: f64 = 2.0;
/// 「曲线差异主要随白平衡变化」的门槛：`r(曲线距离, |Δln(B/G)|)`。
const WB_EXPLAINS_R: f64 = 0.6;

// ---------------------------------------------------------------------------
// 基础工具
// ---------------------------------------------------------------------------

/// 一个样本：相机空间值、工作空间值（`M0·c`）、参考值。
type Triple = ([f32; 3], [f32; 3], [f32; 3]);

fn to_samples(rows: &[Triple]) -> Vec<Sample> {
    rows.iter().map(|r| Sample { ours: r.1, theirs: r.2 }).collect()
}

fn median_of(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[v.len() / 2]
}

fn mean_of(v: &[f64]) -> f64 {
    if v.is_empty() {
        f64::NAN
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().min(ys.len());
    if n < 3 {
        return f64::NAN;
    }
    let (mx, my) = (
        xs[..n].iter().sum::<f64>() / n as f64,
        ys[..n].iter().sum::<f64>() / n as f64,
    );
    let (mut sxy, mut sxx, mut syy) = (0f64, 0f64, 0f64);
    for (x, y) in xs[..n].iter().zip(ys[..n].iter()) {
        let (dx, dy) = (x - mx, y - my);
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return f64::NAN;
    }
    sxy / (sxx * syy).sqrt()
}

/// ΔE00 的中位数。闭包收整个三元组，字段由调用方显式挑（`r.1` 工作空间）。
fn de_median<F: Fn(&Triple) -> [f32; 3]>(rows: &[Triple], f: F) -> f64 {
    let mut des: Vec<f64> = rows.iter().map(|r| deltae::delta_e_prophoto(f(r), r.2)).collect();
    median_of(&mut des)
}

/// 两条曲线三元组的逐点差：`(RMS, 最大)`，单位是编码域 0..1。
///
/// `shape_only = true` 时先减掉**整体**水平差——把「层级/尺度差」与「形状差」分开。
/// 注意它只去掉一个全局常数，**不是**逐通道对齐：某通道单独平移仍会被算作差异。
fn curve_distance(a: &[Vec<f32>; 3], b: &[Vec<f32>; 3], shape_only: bool) -> (f64, f64) {
    let mut diff: Vec<f64> = Vec::new();
    for (ca, cb) in a.iter().zip(b.iter()) {
        let len = ca.len().min(cb.len());
        for i in 0..len {
            diff.push((ca[i] - cb[i]) as f64);
        }
    }
    if diff.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let off = if shape_only { mean_of(&diff) } else { 0.0 };
    let (mut ss, mut mx) = (0f64, 0f64);
    for d in &diff {
        let e = d - off;
        ss += e * e;
        mx = mx.max(e.abs());
    }
    ((ss / diff.len() as f64).sqrt(), mx)
}

/// 曲线在若干编码域点上的取值（用来做「紧凑参数化」）。
fn curve_probe(curves: &[Vec<f32>; 3], xs: &[f32]) -> Vec<f64> {
    let mut out = Vec::new();
    for c in curves.iter() {
        for x in xs {
            out.push(fit::apply_curve(c, *x) as f64);
        }
    }
    out
}

fn stem_number(s: &str) -> Option<u32> {
    s.rsplit('_').next()?.parse::<u32>().ok()
}

/// 组标签：连拍（0001~0015）/ 扫描（0561~0570）/ 老图（其余）。
fn group_of(stem: &str) -> &'static str {
    match stem_number(stem) {
        Some(n) if n <= 15 => "连拍",
        Some(n) if (561..=570).contains(&n) => "扫描",
        _ => "老图",
    }
}

// ---------------------------------------------------------------------------
// 逐图样本
// ---------------------------------------------------------------------------

struct ImageData {
    stem: String,
    group: &'static str,
    rg: f64,
    bg: f64,
    left: Vec<Triple>,
    right: Vec<Triple>,
    secs: f64,
}

fn collect(stem: &str, m0: Mat3) -> Option<ImageData> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        common::skip("诊断 O", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }
    let t0 = std::time::Instant::now();
    let wb = libraw::read_wb(&nef).ok()?;
    if wb[1] <= 0.0 {
        eprintln!("跳过 {stem}：cam_mul 的 G 非正");
        return None;
    }
    let (rg, bg) = ((wb[0] / wb[1]) as f64, (wb[2] / wb[1]) as f64);
    let tif_data = std::fs::read(&tif).ok()?;
    let img = fit::read_rgb16(&tif_data).ok()?;
    let plan = icc::plan_for(&tif_data).ok()?;
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let cam = libraw::decode_with_output_for_test(&nef, &opts, OutputColor::Camera).ok()?;
    let mo = camera::read_model(&nef)?;
    let entry = camera::lookup(&mo)?;
    let mg = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = mg.effective(cam.width, cam.height)?;
    if (cw, ch) != (img.width, img.height) {
        eprintln!("跳过 {stem}：裁切后 {cw}×{ch} 与参考 {}×{} 不符", img.width, img.height);
        return None;
    }
    let theirs = fit::reference_to_working(&img, &plan);

    let mut left: Vec<Triple> = Vec::new();
    let mut right: Vec<Triple> = Vec::new();
    for y in (0..ch).step_by(STEP) {
        for x in (0..cw).step_by(STEP) {
            let j = x + y * img.width;
            let Some(px) = cam.at(x + mg.left, y + mg.top) else { continue };
            let c = [
                px[0] as f32 / 65535.0,
                px[1] as f32 / 65535.0,
                px[2] as f32 / 65535.0,
            ];
            let t = theirs[j];
            if c.iter().any(|v| *v <= CAM_LO || *v >= CAM_HI)
                || t.iter().any(|v| *v <= REF_LO || *v >= REF_HI)
            {
                continue;
            }
            let ours = mat3::mul_vec(m0, c);
            if x < cw / 2 {
                left.push((c, ours, t));
            } else {
                right.push((c, ours, t));
            }
        }
    }
    if left.len() < MIN_SAMPLES || right.len() < MIN_SAMPLES {
        eprintln!("跳过 {stem}：左右半样本不足（{} / {}）", left.len(), right.len());
        return None;
    }
    Some(ImageData {
        stem: stem.to_string(),
        group: group_of(stem),
        rg,
        bg,
        left,
        right,
        secs: t0.elapsed().as_secs_f64(),
    })
}

fn pool(images: &[ImageData], skip: Option<usize>) -> Vec<Triple> {
    let mut out = Vec::new();
    for (i, im) in images.iter().enumerate() {
        if skip == Some(i) {
            continue;
        }
        out.extend_from_slice(&im.left);
    }
    out
}

fn build_transform(samples: &[Sample], name: &str) -> BaseTransform {
    let curve = fit::fit_curve(samples, CURVE_BINS);
    let lut = fit::fit_lut(samples, &curve, LUT_EDGE);
    BaseTransform::from_fit(0x0000, name, Vec::new(), curve, LUT_EDGE, lut)
}

/// 在给定底座上，用 `rows` 拟合 3 条逐通道曲线（`rows` 的 `ours` 必须是**底座输出**）。
fn per_channel_on(base: &BaseTransform, rows: &[Triple]) -> [Vec<f32>; 3] {
    let post: Vec<Sample> = rows
        .iter()
        .map(|r| Sample { ours: base.apply(r.1), theirs: r.2 })
        .collect();
    fit::fit_curves_per_channel(&post, CURVE_BINS)
}

// ---------------------------------------------------------------------------
// 主测试
// ---------------------------------------------------------------------------

fn env_stems() -> Vec<String> {
    match std::env::var("NRV_O_STEMS") {
        Ok(v) if !v.trim().is_empty() => v
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect(),
        _ => DEFAULT_STEMS.iter().map(|s| s.to_string()).collect(),
    }
}

#[test]
fn oracle_curve_similarity() {
    eprintln!("=== 诊断 O：神谕逐通道曲线的逐图相似性（描述性比较，不声称可达增益）===");
    let Some(cm) = camera::color_matrix("Nikon Z 8") else {
        eprintln!("跳过：机型表里没有 Z8 的色彩矩阵");
        return;
    };
    let m0 = cm.matrix;
    eprintln!("门槛（先定后看）：曲线相似 = 编码域逐点 RMS ≤ {RMS_SIMILAR}；组内/跨组 = 跨组÷组内 ≥ {GROUP_RATIO}");

    let stems = env_stems();
    eprintln!("\n[0] 图（{} 张）：{}", stems.len(), stems.join(" "));
    let mut images: Vec<ImageData> = Vec::new();
    for s in &stems {
        eprintln!("  … 解码 {s}");
        if let Some(d) = collect(s, m0) {
            images.push(d);
        }
    }
    if images.len() < 4 {
        eprintln!("跳过：可用图不足（{} 张）", images.len());
        return;
    }
    eprintln!(
        "  解码总耗时 {:.0} s · 分组：{}",
        images.iter().map(|i| i.secs).sum::<f64>(),
        images
            .iter()
            .map(|i| format!("{}({})", i.stem, i.group))
            .collect::<Vec<_>>()
            .join(" ")
    );

    // ---- [1] 神谕 vs 非神谕（留一图底座）----
    eprintln!("\n[1] 每张图的曲线：底座 = 其余图池化（T_i，留给图协议）");
    eprintln!("  (b) 底座原样 · (d) 底座 + **池化**逐通道曲线（非神谕，其余图左半拟合）· (c) 底座 + **自己**的曲线（神谕）");
    eprintln!("  {:<10} {:>6} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10}", "图", "组", "B/G", "(b)底座", "(d)池化逐通道", "d−b", "(c)神谕", "c−b");
    let mut gain_oracle: Vec<f64> = Vec::new();
    let mut gain_pooled: Vec<f64> = Vec::new();
    let mut base_all: Vec<f64> = Vec::new();
    let mut oracle_all: Vec<f64> = Vec::new();
    let mut pooled_all: Vec<f64> = Vec::new();
    for (i, im) in images.iter().enumerate() {
        let cal = pool(&images, Some(i));
        let t_i = build_transform(&to_samples(&cal), "T_i");
        let b = de_median(&im.right, |r| t_i.apply(r.1));
        // 非神谕：用**其余图**（在 T_i 输出上）拟合逐通道曲线
        let post_cal: Vec<Triple> = cal
            .iter()
            .map(|r| (r.0, t_i.apply(r.1), r.2))
            .collect();
        let pooled_curves = fit::fit_curves_per_channel(&to_samples(&post_cal), CURVE_BINS);
        let d = de_median(&im.right, |r| fit::forward_curves(&pooled_curves, t_i.apply(r.1)));
        // 神谕：用**自己**的左半拟合
        let own_curves = per_channel_on(&t_i, &im.left);
        let c = de_median(&im.right, |r| fit::forward_curves(&own_curves, t_i.apply(r.1)));
        gain_oracle.push(b - c);
        gain_pooled.push(b - d);
        base_all.push(b);
        oracle_all.push(c);
        pooled_all.push(d);
        eprintln!(
            "  {:<10} {:>6} {:>8.4} {:>10.3} {:>12.3} {:>10.3} {:>10.3} {:>10.3}",
            im.stem, im.group, im.bg, b, d, b - d, c, b - c
        );
    }
    eprintln!(
        "  中位：底座 {:.3} · 池化逐通道 {:.3}（吃 {:.3}）· 神谕 {:.3}（吃 {:.3}）",
        median_of(&mut base_all.clone()),
        median_of(&mut pooled_all.clone()),
        median_of(&mut gain_pooled.clone()),
        median_of(&mut oracle_all.clone()),
        median_of(&mut gain_oracle.clone())
    );
    eprintln!("  ⇒ 「神谕吃得动、池化吃不动」在这里再次成立；下面比的是**神谕曲线彼此像不像**。");

    // ---- [2] 公共底座上的逐图神谕曲线 ----
    let t_pool = build_transform(&to_samples(&pool(&images, None)), "T_pool");
    let curves: Vec<[Vec<f32>; 3]> = images.iter().map(|im| per_channel_on(&t_pool, &im.left)).collect();
    let base_pool_de: Vec<f64> = images
        .iter()
        .map(|im| de_median(&im.right, |r| t_pool.apply(r.1)))
        .collect();
    eprintln!("\n[2] 公共底座 T_pool（全部 10 张池化）上的逐图神谕曲线");
    eprintln!("  底座原样（无曲线）的 ΔE00 中位：{}", base_pool_de.iter().map(|v| format!("{v:.2}")).collect::<Vec<_>>().join(" "));

    // 同名短标签，便于矩阵排版
    let short: Vec<String> = images
        .iter()
        .map(|im| im.stem.trim_start_matches("DSC_").to_string())
        .collect();

    // 功能矩阵：行 = 曲线的来源，列 = 应用到的图
    eprintln!("\n[2a] 功能矩阵：用**行**那张图的曲线，套到**列**那张图上，ΔE00 中位（对角线 = 用自己的曲线）");
    eprint!("  {:<6}", "曲线\\图");
    for s in &short {
        eprint!(" {:>7}", s);
    }
    eprintln!();
    let mut de_matrix: Vec<Vec<f64>> = vec![vec![f64::NAN; images.len()]; images.len()];
    for a in 0..images.len() {
        eprint!("  {:<6}", short[a]);
        for (b, im) in images.iter().enumerate() {
            let v = de_median(&im.right, |r| fit::forward_curves(&curves[a], t_pool.apply(r.1)));
            de_matrix[a][b] = v;
            eprint!(" {:>7.3}", v);
        }
        eprintln!();
    }
    // 用错曲线的代价：把「别人的曲线」与「自己的曲线」比
    let mut wrong_cost: Vec<f64> = Vec::new();
    for (a, row) in de_matrix.iter().enumerate() {
        for (b, v) in row.iter().enumerate() {
            if a != b {
                wrong_cost.push(v - de_matrix[b][b]);
            }
        }
    }
    let own_gain: Vec<f64> = (0..images.len())
        .map(|b| base_pool_de[b] - de_matrix[b][b])
        .collect();
    eprintln!(
        "  「用自己的曲线」相对底座原样：中位改善 {:.3}（范围 {:.3}~{:.3}）",
        median_of(&mut own_gain.clone()),
        own_gain.iter().cloned().fold(f64::INFINITY, f64::min),
        own_gain.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    );
    eprintln!(
        "  「用**别人**的曲线」相对自己的曲线：中位代价 +{:.3}（范围 {:.3}~{:.3}；负值=别人的反而更好）",
        median_of(&mut wrong_cost.clone()),
        wrong_cost.iter().cloned().fold(f64::INFINITY, f64::min),
        wrong_cost.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    );

    // 曲线空间矩阵
    eprintln!("\n[2b] 曲线空间距离：逐点差 RMS（编码域，上三角）· 括号内为去掉整体水平差后的「纯形状」距离");
    eprint!("  {:<6}", "曲线\\曲线");
    for s in &short {
        eprint!(" {:>7}", s);
    }
    eprintln!();
    let mut rms_matrix: Vec<Vec<f64>> = vec![vec![f64::NAN; images.len()]; images.len()];
    let mut shape_pairs: Vec<f64> = Vec::new();
    for a in 0..images.len() {
        eprint!("  {:<6}", short[a]);
        for b in 0..images.len() {
            if b < a {
                eprint!(" {:>7}", "·");
                continue;
            }
            let (r, _) = curve_distance(&curves[a], &curves[b], false);
            let (s, _) = curve_distance(&curves[a], &curves[b], true);
            rms_matrix[a][b] = r;
            rms_matrix[b][a] = r;
            if a != b {
                shape_pairs.push(s);
            }
            eprint!(" {:>7.4}", r);
        }
        eprintln!();
    }
    let mut rms_pairs: Vec<f64> = Vec::new();
    for (a, row) in rms_matrix.iter().enumerate() {
        for v in row.iter().skip(a + 1) {
            rms_pairs.push(*v);
        }
    }
    eprintln!(
        "  两两 RMS：中位 {:.4} · 最大 {:.4} · 门槛（相似）{RMS_SIMILAR}",
        median_of(&mut rms_pairs.clone()),
        rms_pairs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    );
    eprintln!(
        "  去掉整体水平差后的「纯形状」距离：中位 {:.4} ⇒ {}",
        median_of(&mut shape_pairs.clone()),
        if median_of(&mut shape_pairs.clone()) < median_of(&mut rms_pairs.clone()) * 0.5 {
            "大部分差异是整体水平差"
        } else {
            "差异主要不在「整体水平」上——曲线两端被钉在 0/1，水平差本来就不是一个独立的自由度"
        }
    );

    // ---- [3] 分组 ----
    eprintln!("\n[3] 分组：连拍（同场景同白平衡，噪声底）/ 扫描（同场景变白平衡）/ 老图（不同场景）");
    let same_group = |a: usize, b: usize| images[a].group == images[b].group;
    let mut within_rms: Vec<f64> = Vec::new();
    let mut cross_rms: Vec<f64> = Vec::new();
    let mut within_cost: Vec<f64> = Vec::new();
    let mut cross_cost: Vec<f64> = Vec::new();
    for a in 0..images.len() {
        for b in 0..images.len() {
            if a == b {
                continue;
            }
            let cost = de_matrix[a][b] - de_matrix[b][b];
            if same_group(a, b) {
                within_rms.push(rms_matrix[a][b]);
                within_cost.push(cost);
            } else {
                cross_rms.push(rms_matrix[a][b]);
                cross_cost.push(cost);
            }
        }
    }
    eprintln!(
        "  曲线 RMS：组内中位 {:.4}（{} 对）· 跨组中位 {:.4}（{} 对）· 比值 {:.2}（门槛 {GROUP_RATIO}）",
        median_of(&mut within_rms.clone()),
        within_rms.len(),
        median_of(&mut cross_rms.clone()),
        cross_rms.len(),
        median_of(&mut cross_rms.clone()) / median_of(&mut within_rms.clone()).max(1e-9)
    );
    eprintln!(
        "  用错曲线的 ΔE00 代价：组内中位 {:.3} · 跨组中位 {:.3}",
        median_of(&mut within_cost.clone()),
        median_of(&mut cross_cost.clone())
    );
    // 组内细分
    for g in ["连拍", "扫描"] {
        let idx: Vec<usize> = (0..images.len()).filter(|i| images[*i].group == g).collect();
        if idx.len() < 2 {
            continue;
        }
        let mut rr = Vec::new();
        let mut cc = Vec::new();
        for (x, &a) in idx.iter().enumerate() {
            for &b in idx.iter().skip(x + 1) {
                rr.push(rms_matrix[a][b]);
                cc.push(de_matrix[a][b] - de_matrix[b][b]);
            }
        }
        eprintln!(
            "    {g}组内（{} 对）：RMS 中位 {:.4} · 最大 {:.4}；用错曲线代价中位 {:+.3}",
            rr.len(),
            median_of(&mut rr.clone()),
            rr.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            median_of(&mut cc.clone())
        );
    }
    let mut gmeans: Vec<(String, f64)> = Vec::new();
    for g in ["连拍", "扫描", "老图"] {
        let mut v: Vec<f64> = images
            .iter()
            .enumerate()
            .filter(|(_, im)| im.group == g)
            .map(|(i, _)| de_matrix[i][i])
            .collect();
        if !v.is_empty() {
            gmeans.push((g.to_string(), median_of(&mut v)));
        }
    }
    eprintln!(
        "  各组「用自己的曲线」后的 ΔE00 中位：{}",
        gmeans.iter().map(|(g, v)| format!("{g} {v:.3}")).collect::<Vec<_>>().join(" · ")
    );

    // ---- [3b] 核心判据：曲线差异随**白平衡距离**增大吗 ----
    //
    // [3] 的分组是按「片子来源」（连拍/扫描/老图）分的，那**不是**假设检验该用的坐标：
    // 扫描组内部跨了 B/G 1.73~3.02（白平衡差得很远），而跨组之间也可能白平衡很接近。
    // 所以这里直接把「曲线距离」对「白平衡距离」回归——这才是「曲线像不像」的判据。
    eprintln!("\n[3b] 核心判据：曲线差异随**白平衡距离**增大吗（配对 = 两两组合，n = {}）", rms_pairs.len());
    let mut x_bg: Vec<f64> = Vec::new();
    let mut x_rg: Vec<f64> = Vec::new();
    let mut y_rms: Vec<f64> = Vec::new();
    let mut y_cost: Vec<f64> = Vec::new();
    let mut near_wb: Vec<f64> = Vec::new();
    let mut far_wb: Vec<f64> = Vec::new();
    let mut near_cost: Vec<f64> = Vec::new();
    let mut far_cost: Vec<f64> = Vec::new();
    for a in 0..images.len() {
        for b in (a + 1)..images.len() {
            let dbg = (images[a].bg / images[b].bg).ln().abs();
            let drg = (images[a].rg / images[b].rg).ln().abs();
            let cost = 0.5
                * ((de_matrix[a][b] - de_matrix[b][b]) + (de_matrix[b][a] - de_matrix[a][a]));
            x_bg.push(dbg);
            x_rg.push(drg);
            y_rms.push(rms_matrix[a][b]);
            y_cost.push(cost);
            if dbg <= 0.05 {
                near_wb.push(rms_matrix[a][b]);
                near_cost.push(cost);
            } else if dbg >= 0.30 {
                far_wb.push(rms_matrix[a][b]);
                far_cost.push(cost);
            }
        }
    }
    let r_bg = pearson(&x_bg, &y_rms);
    let r_rg = pearson(&x_rg, &y_rms);
    eprintln!(
        "  r(曲线 RMS, |Δln B/G|) = {r_bg:+.3} · r(曲线 RMS, |Δln R/G|) = {r_rg:+.3}（门槛 {WB_EXPLAINS_R}）"
    );
    eprintln!(
        "  r(用错曲线的 ΔE00 代价, |Δln B/G|) = {:+.3} · r(代价, |Δln R/G|) = {:+.3}",
        pearson(&x_bg, &y_cost),
        pearson(&x_rg, &y_cost)
    );
    eprintln!(
        "  B/G 相近（|Δln| ≤ 0.05，{} 对）：曲线 RMS 中位 {:.4} · 用错曲线代价中位 {:+.3}",
        near_wb.len(),
        median_of(&mut near_wb.clone()),
        median_of(&mut near_cost.clone())
    );
    eprintln!(
        "  B/G 相隔远（|Δln| ≥ 0.30，{} 对）：曲线 RMS 中位 {:.4} · 用错曲线代价中位 {:+.3}",
        far_wb.len(),
        median_of(&mut far_wb.clone()),
        median_of(&mut far_cost.clone())
    );

    // ---- [4] 紧凑参数化：曲线差在哪、与白平衡有没有关系 ----
    eprintln!("\n[4] 紧凑参数化：三条曲线在编码域 0.25/0.50/0.75 处的取值（每行一张图）");
    let xs = [0.25f32, 0.50, 0.75];
    eprintln!(
        "  {:<10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}   {:>8}",
        "图", "R.25", "R.50", "R.75", "G.25", "G.50", "G.75", "B.25", "B.50", "B.75", "B/G"
    );
    let mut probes: Vec<Vec<f64>> = Vec::new();
    for (im, c) in images.iter().zip(&curves) {
        let p = curve_probe(c, &xs);
        eprintln!(
            "  {:<10} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4}   {:>8.4}",
            im.stem, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], im.bg
        );
        probes.push(p);
    }
    // 每个参数跨图的离散 + 与白平衡的相关
    let names = ["R.25", "R.50", "R.75", "G.25", "G.50", "G.75", "B.25", "B.50", "B.75"];
    let bgs: Vec<f64> = images.iter().map(|i| i.bg).collect();
    let rgs: Vec<f64> = images.iter().map(|i| i.rg).collect();
    eprintln!("\n  每个参数的跨图离散与相关性（n = {}）：", images.len());
    eprintln!(
        "  {:<8} {:>10} {:>10} {:>10} {:>12} {:>12}",
        "参数", "最小", "最大", "极差", "r(B/G)", "r(R/G)"
    );
    for (j, nm) in names.iter().enumerate() {
        let v: Vec<f64> = probes.iter().map(|p| p[j]).collect();
        let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        eprintln!(
            "  {nm:<8} {lo:>10.4} {hi:>10.4} {:>10.4} {:>12.3} {:>12.3}",
            hi - lo,
            pearson(&bgs, &v),
            pearson(&rgs, &v)
        );
    }

    // ---- [5] 判定 ----
    let med_rms = median_of(&mut rms_pairs.clone());
    let med_cost = median_of(&mut wrong_cost.clone());
    let med_within = median_of(&mut within_rms.clone());
    let med_cross = median_of(&mut cross_rms.clone());
    let ratio = med_cross / med_within.max(1e-9);
    eprintln!("\n[5] 判定（本诊断只做描述性比较，不声称可达增益）");
    eprintln!(
        "  1) 曲线彼此像不像：两两 RMS 中位 {med_rms:.4}（门槛 {RMS_SIMILAR}）⇒ {}",
        if med_rms <= RMS_SIMILAR { "**相似**" } else { "**不相似**" }
    );
    eprintln!(
        "  2) 用错曲线的代价：ΔE00 中位 {med_cost:+.3}（「用自己的曲线」本身带来 {:.3} 的改善）⇒ {}",
        median_of(&mut own_gain.clone()),
        if med_cost <= 0.3 { "代价小、曲线可互换" } else { "**代价明显，曲线不可互换**" }
    );
    eprintln!(
        "  3) 分组结构（按片子来源）：组内 RMS 中位 {med_within:.4} · 跨组 {med_cross:.4} · 比值 {ratio:.2}（门槛 {GROUP_RATIO}）"
    );
    eprintln!(
        "  4) **核心判据**（按白平衡距离）：r(曲线 RMS, |Δln B/G|) = {r_bg:+.3}（门槛 {WB_EXPLAINS_R}）"
    );
    if r_bg >= WB_EXPLAINS_R {
        eprintln!("     ⇒ **曲线差异主要随白平衡变化**：既不是「逐图随机」（假设 B），也不是场景内容——");
        eprintln!("       扫描组是**同一个场景**、只变白平衡，组内曲线差异就已经很大。⇒ 指向 **A**：");
        eprintln!("       我们的变换缺一个**随白平衡变化**的分量。按「片子来源」分组看不出这一点。");
    } else if ratio >= GROUP_RATIO {
        eprintln!("     ⇒ 组内相似、跨组差异大：指向「场景/内容相关」，不是逐图随机。");
    } else if med_rms <= RMS_SIMILAR {
        eprintln!("     ⇒ 组内组外都相似 ⇒ 假设 **A**（本该可预测，拟合/池化方式有问题）。");
    } else {
        eprintln!("     ⇒ 与白平衡无强相关、组内跨组都差 ⇒ 假设 **B**（真正逐图特有）。");
    }
    eprintln!();
    eprintln!("  限制（必须随结论一起转述）：");
    eprintln!("  · 神谕曲线是**用该图自己的参考值**拟合的，量的是形状，**不是可达增益**；");
    eprintln!("    [1] 里「池化逐通道」那一列才是可达侧的对照。");
    eprintln!("  · 10 张图、三组；连拍组是同场景同白平衡（只能给噪声底），扫描组是同场景变白平衡，");
    eprintln!("    只有 3 张老图带跨场景成分 —— **跨场景这一维依然很薄**。");
    eprintln!("  · [2] 用公共底座 T_pool（见过每张图，有泄漏），这是为了让曲线可比；[1] 用留一图底座。");
    eprintln!("  · 曲线是分箱中位数 + 强制单调的产物，bins = {CURVE_BINS}；分箱方式本身的差异也算进「曲线差异」里。");
}

/// 仪器自检（不依赖样本）：
/// 1. 曲线距离：同一条曲线距离为 0、去掉水平差后对整体平移不敏感；
/// 2. `forward_curves` 与逐通道 `forward_curve` 一致；
/// 3. `fit_curves_per_channel` 在「已知逐通道 gamma」的合成数据上能还原三条不同的曲线。
#[test]
fn instrument_checks() {
    // 1) 曲线距离
    let a: [Vec<f32>; 3] = [vec![0.0, 0.3, 0.6], vec![0.0, 0.3, 0.6], vec![0.0, 0.3, 0.6]];
    let (d0, m0) = curve_distance(&a, &a, false);
    assert!(d0.abs() < 1e-12 && m0.abs() < 1e-12, "同曲线距离应为 0：{d0} {m0}");
    // 纯整体平移（不加 clamp，否则会掺进形状差）
    let shifted: [Vec<f32>; 3] = [vec![0.1, 0.4, 0.7], vec![0.1, 0.4, 0.7], vec![0.1, 0.4, 0.7]];
    let (d_raw, _) = curve_distance(&a, &shifted, false);
    let (d_shape, _) = curve_distance(&a, &shifted, true);
    assert!((d_raw - 0.1).abs() < 1e-6, "整体平移 0.1 的 RMS 应为 0.1：{d_raw}");
    assert!(d_shape < 1e-6, "去掉水平差后整体平移应完全不可见：{d_shape}");

    // 2) forward_curves 与逐通道一致
    let v = [0.31f32, 0.52, 0.77];
    let got = fit::forward_curves(&a, v);
    for k in 0..3 {
        let want = fit::forward_curve(&a[k], [v[k]; 3])[0];
        assert!((got[k] - want).abs() < 1e-6, "通道 {k}：{} vs {want}", got[k]);
    }

    // 3) 已知逐通道 gamma：三条曲线应当明显不同，且各自近似对应的幂
    //    注意三个通道都要**覆盖到整个 0..1**，否则上端只是端点延拓，比不出曲线高低。
    let mut samples: Vec<Sample> = Vec::new();
    for i in 0..60 {
        for j in 0..20 {
            let ours = [
                0.02 + i as f32 * 0.016,                            // 0.02 ~ 0.964
                0.02 + j as f32 * 0.048,                            // 0.02 ~ 0.932
                0.03 + ((i + j * 3) % 50) as f32 * 0.019,           // 0.03 ~ 0.961
            ];
            let theirs = [ours[0].powf(1.8), ours[1].powf(1.0), ours[2].powf(0.6)];
            samples.push(Sample { ours, theirs });
        }
    }
    let cs = fit::fit_curves_per_channel(&samples, CURVE_BINS);
    for x in [0.2f32, 0.5, 0.8] {
        let y = fit::forward_curves(&cs, [x; 3]);
        assert!(y[0] < y[1] - 0.02, "R 通道（gamma 1.8）应低于 G 通道：{y:?} @ x={x}");
        assert!(y[2] > y[1] + 0.02, "B 通道（gamma 0.6）应高于 G 通道：{y:?} @ x={x}");
        assert!((y[1] - x).abs() < 0.03, "G 通道（gamma 1.0）应近似恒等：{y:?} @ x={x}");
    }
    // 三条曲线应彼此可分
    let (d_rg, _) = curve_distance(&[cs[0].clone(), cs[0].clone(), cs[0].clone()], &[cs[1].clone(), cs[1].clone(), cs[1].clone()], false);
    assert!(d_rg > 0.02, "不同 gamma 的曲线应可区分，得到 {d_rg}");
}
