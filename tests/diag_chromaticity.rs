//! 诊断 D：**尺度不变**的色度比对——我们究竟「偏淡」，还是「只是偏暗」。
//!
//! # 要回答的问题
//!
//! 三张参考图的实测（`tasks.md` 5.6a）：残余误差已从明度转向彩度，而且**我们解码的
//! C\* 中位数只有参考的一半左右**（DSC_0141 9.92 vs 18.57、DSC_8562 5.46 vs 8.69）。
//! 但「偏淡」与「偏暗」是两种完全不同的病、修法也不同，而那张表**分不开它们**。
//!
//! # 为什么 C\* 分不开
//!
//! Lab 的 a\*、b\* 是 XYZ 经立方根（暗部是线性段）非线性之后的结果，**同一色度在更亮的
//! 像素上读出更大的 C\***。拿「亮度还没对齐」的两组像素去比 C\*，等于把「谁更亮」混进
//! 了「谁更艳」里——而我们已知中间调偏暗，这个混淆**恰好指向「偏淡」这个错误结论**。
//!
//! 本测试用同一批像素实测这个弹性：把像素整体乘 k 使它的 L\* 对齐到参考的 L\*，再看 C\*
//! 变了多少。弹性 ≈ 1 表示该亮度段里 C\* 与亮度几乎成正比（Lab 线性段，深阴影就是），
//! ≈ 1/3 表示已进入立方根段。
//!
//! # 判据：线性域里尺度不变的通道比值
//!
//! `ln(R/G)`、`ln(B/G)` 在**任意整体增益**下都不变：两张图各自乘任何正常数，这两个数
//! 一个都不动。于是逐像素相减得到的差**不可能**是「谁更亮」造成的，取中位数后也不用
//! 担心错位带来的对称噪声。再按**参考亮度分箱**看这个差是否在各段持续存在：
//!
//! - 各亮度箱内比值差都显著非零 ⇒ **确实是色度错**，把亮度对齐也修不掉；
//! - 各箱内都趋于零、而 C\* 差很大 ⇒ **只是亮度差造成的假象**，该继续修明度而非色度。
//!
//! 还有一层必须交代清楚：**逐通道各自一条单调曲线**（= 现在的 1D 曲线能做的事，也是「每个
//! 通道各自偏暗/偏亮」能造出的全部效果）同样会改变通道比值——整体增益让掉之后剩下的差，
//! 未必是「跨通道」的色度错。所以再补一步控制：用 3 条独立单调曲线把 ours 逐通道映射到
//! theirs，看比值差是否塌到 1.0000。塌了 ⇒ 这个差与逐通道模型相容、与亮度差同源；没塌 ⇒
//! 才需要跨通道的修法（色彩矩阵 / 3D LUT / 饱和度）。
//!
//! 判定一律用 `DiagSample::ours_linear`（**未经基准变换的原始解码**）：要判的是解码
//! 环节，不是变换之后。变换后的口径另出一段作对照，看现在的曲线 + LUT 把这个差修掉了
//! 多少。
//!
//! 样本缺席时跳过而非失败。只跑本文件：
//! `cargo test --offline --test diag_chromaticity -- --nocapture`。

mod common;

use common::{fitted_fixture_for, quantiles, DiagSample, FittedFixture};
use nikonrawview::deltae;
use nikonrawview::fit::encode_srgb;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// 参数
// ---------------------------------------------------------------------------

/// 参与诊断的图。`DSC_0001` 是竖构图、画面接近中性的对照；另两张是场景完全不同的横构图。
const STEMS: [&str; 3] = ["DSC_0001", "DSC_0141", "DSC_8562"];

/// 亮度分箱数，键 = **参考**亮度的 sRGB 编码值（与 `fit::fit_curve` 的定义域一致）。
const BINS: usize = 8;

/// 比值下限：16 位量化下 1 个计数是 1/65535，只有几个计数的值其比值由量化噪声主导。
const MIN_SIGNAL: f32 = 4.0 / 65535.0;

/// 比值上限：解码与参考都在 1.0 处剪切，剪切过的通道其比值已经失真。
const MAX_SIGNAL: f32 = 0.999;

/// 判定「比值差可察觉」的乘性阈值。
///
/// 取 1%：实测在 ProPhoto→Lab 尺度上，R/G 差 1% 折算到中性锚点已有 ΔE00 约 1.0 的量级
/// （见每箱末列 `dE_ratio`），与达标线 ≤1.0 同级。阈值定得比它松就测不出问题。
const RATIO_EPS: f64 = 0.01;

/// 一个箱至少要有这么多个通过信号门槛的样本才参与判定。
const MIN_BIN: usize = 200;

/// 在 (ln R/G, ln B/G) 平面上，参考彩度方向短于这个长度时谈「沿彩度方向」没有意义。
const MIN_CHROMA_MAG: f64 = 0.05;

/// 弹性只在彩度足够大的像素上算：接近中性的像素 a\*、b\* 都是小量，比值不稳。
const MIN_C_FOR_ELASTICITY: f64 = 0.5;

/// 「逐通道单调映射」控制的取样箱数（sRGB 编码域等距，与 `fit::fit_curve` 同域同量级）。
const MONO_BINS: usize = 128;

// ---------------------------------------------------------------------------
// 基本量
// ---------------------------------------------------------------------------

/// ProPhoto 线性 → CIELAB(D50)。
fn lab(v: [f32; 3]) -> [f32; 3] {
    deltae::lab_from_prophoto_linear(v)
}

/// CIE L\*。
fn lstar(v: [f32; 3]) -> f64 {
    lab(v)[0] as f64
}

/// Lab 彩度 C\*。
fn chroma(v: [f32; 3]) -> f64 {
    let l = lab(v);
    ((l[1] as f64).powi(2) + (l[2] as f64).powi(2)).sqrt()
}

/// Rec.709 亮度（与脚手架 `DiagSample::lum_ours` 同一口径）。
fn lum(v: [f32; 3]) -> f64 {
    (0.2126 * v[0] + 0.7152 * v[1] + 0.0722 * v[2]) as f64
}

/// 求缩放系数 k，使 `k·v` 的 L\* 等于 `target`（L\* 对 k 单调，二分即可）。
///
/// 这是本次诊断的核心工具：C\* 不是尺度不变量，只有把两边的**明度对齐**之后再比 C\*，
/// 读出的彩度差才不是亮度差的影子。
fn scale_for_lstar(v: [f32; 3], target: f64) -> Option<f64> {
    if !target.is_finite() || v.iter().any(|&x| !x.is_finite()) || v.iter().all(|&x| x <= 0.0) {
        return None;
    }
    let at = |k: f64| lstar(v.map(|x| (x as f64 * k) as f32));
    let (mut lo, mut hi) = (1e-4_f64, 1e4_f64);
    if !(at(lo) < target && target < at(hi)) {
        return None;
    }
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if at(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// 把像素整体缩放到明度对齐后，再读 C\*。
fn chroma_scaled(v: [f32; 3], k: f64) -> f64 {
    chroma(v.map(|x| (x as f64 * k) as f32))
}

/// 一个像素的**尺度不变**色度差。
///
/// 全部由通道比值的对数差构成：把两张图各自整体乘任意正常数，这里一个数都不变，因此它
/// 读出的差异不可能来自「谁更亮」。
#[derive(Debug, Clone, Copy)]
struct Shift {
    /// `ln(R/G)`：我们 − 参考。
    d_rg: f64,
    /// `ln(B/G)`：我们 − 参考。
    d_bg: f64,
    /// 差向量长度：尺度不变的色度偏差总幅度。
    d_mag: f64,
    /// 沿参考彩度方向的分量：正 = 我们更饱和，负 = 我们更淡。
    radial: f64,
    /// 垂直于参考彩度方向的分量：色相偏移。
    tangential: f64,
}

/// 计算尺度不变的通道比值差；未过信号门槛时返回 `None`（该像素不参与比值统计）。
///
/// 门槛两边都要过：任一侧有通道低于 `MIN_SIGNAL`（量化噪声主导）或高于 `MAX_SIGNAL`
/// （已剪切）时，通道比值就不再代表色度。
fn shift_of(o: [f32; 3], t: [f32; 3]) -> Option<Shift> {
    let usable =
        |v: [f32; 3]| v.iter().all(|&x| x.is_finite() && (MIN_SIGNAL..=MAX_SIGNAL).contains(&x));
    if !usable(o) || !usable(t) {
        return None;
    }
    let axes = |v: [f32; 3]| {
        (
            (v[0] as f64 / v[1] as f64).ln(),
            (v[2] as f64 / v[1] as f64).ln(),
        )
    };
    let (u_t, v_t) = axes(t);
    let (u_o, v_o) = axes(o);
    let d_rg = u_o - u_t;
    let d_bg = v_o - v_t;
    let mag = (u_t * u_t + v_t * v_t).sqrt();
    let (radial, tangential) = if mag >= MIN_CHROMA_MAG {
        (
            (d_rg * u_t + d_bg * v_t) / mag,
            (u_t * d_bg - v_t * d_rg) / mag,
        )
    } else {
        (f64::NAN, f64::NAN)
    };
    Some(Shift {
        d_rg,
        d_bg,
        d_mag: (d_rg * d_rg + d_bg * d_bg).sqrt(),
        radial,
        tangential,
    })
}

/// 把比值差折算成同一把尺子上的读数。
///
/// 在中性锚点上把 R、B 各乘测得的因子，与未扰动的灰求 ΔE00：这样「R/G 差 1%」就有了
/// 与达标线可比的量级（ΔE00 本身含扰动引起的微量明度变化，故不拆项）。
fn anchor_delta_e(y: f64, d_rg: f64, d_bg: f64) -> f64 {
    let g = y.clamp(MIN_SIGNAL as f64, 1.0) as f32;
    let base = [g, g, g];
    let shifted = [
        (g as f64 * d_rg.exp()) as f32,
        g,
        (g as f64 * d_bg.exp()) as f32,
    ];
    deltae::delta_e_prophoto(shifted, base)
}

// ---------------------------------------------------------------------------
// 分布统计
// ---------------------------------------------------------------------------

/// 一个分布：样本数、中位数、P5/P95，以及中位数自身的标准误。
struct Dist {
    n: usize,
    med: f64,
    p05: f64,
    p95: f64,
    /// `1.253·σ/√n`，σ 由四分位距推得（`IQR/1.349`）。用来判断「千分之几」的系统性偏移
    /// 是不是真的——上万样本时它不是。
    se: f64,
}

/// 最近秩分位数，与 `common::quantiles` 同一口径（`ceil(q·n) − 1`）。
///
/// 这里另写一份是因为还需要 P5/P25/P75；分位数的定义必须与仓库其余部分一致——5.5 曾因
/// `round` 与 `ceil` 的差别把中位数报到有偏的一侧。
fn pct(sorted: &[f64], q: f64) -> f64 {
    let rank = (q * sorted.len() as f64).ceil().max(1.0) as usize;
    sorted[rank.saturating_sub(1).min(sorted.len() - 1)]
}

/// 汇总一组数；非有限值（缺失项）直接丢弃，空集返回 `None`。
fn dist(values: &[f64]) -> Option<Dist> {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let sigma = (pct(&v, 0.75) - pct(&v, 0.25)) / 1.349;
    Some(Dist {
        n: v.len(),
        med: pct(&v, 0.5),
        p05: pct(&v, 0.05),
        p95: pct(&v, 0.95),
        se: 1.253 * sigma / (v.len() as f64).sqrt(),
    })
}

/// 取一个分布的中位数；缺失返回 NaN（打印用的中性值）。
fn med_of(d: &Option<Dist>) -> f64 {
    d.as_ref().map_or(f64::NAN, |x| x.med)
}

/// 定宽数值列；缺失打 `—`。
fn col(opt: Option<f64>, prec: usize) -> String {
    match opt {
        Some(v) if v.is_finite() => format!("{v:>w$.p$}", w = prec + 6, p = prec),
        _ => format!("{:>w$}", "—", w = prec + 6),
    }
}

// ---------------------------------------------------------------------------
// 逐样本观测
// ---------------------------------------------------------------------------

/// 一个像素上收集到的观测（`ours` 由调用方选定的口径决定）。
struct Rec {
    /// 分箱键 = 参考亮度的 sRGB 编码值。
    key: f64,
    /// 我们的亮度（变换前 / 变换后，取决于口径）。
    y_o: f64,
    /// 参考亮度。
    y_t: f64,
    /// 我们的 C\*。
    c_o: f64,
    /// 参考的 C\*。
    c_t: f64,
    /// 我们把明度对齐到参考之后的 C\*。
    c_o_aligned: f64,
    /// C\* 对亮度的弹性 `Δln C*/Δln k`（k 为把 L\* 对齐所需的整体缩放）。
    elasticity: f64,
    /// 尺度不变的通道比值差；未过信号门槛时为 `None`。
    shift: Option<Shift>,
}

/// 收集一组样本的观测。`ours` 取 `|s| s.ours_linear`（判定口径）或 `|s| s.ours_transformed`。
fn recs_of(eval: &[DiagSample], ours: fn(&DiagSample) -> [f32; 3]) -> Vec<Rec> {
    eval.iter()
        .map(|s| {
            let o = ours(s);
            let t = s.theirs_linear;
            let c_o = chroma(o);
            let c_t = chroma(t);
            let k_up = scale_for_lstar(o, lstar(t));
            let c_o_aligned = k_up.map_or(f64::NAN, |k| chroma_scaled(o, k));
            let elasticity = match k_up {
                Some(k)
                    if c_o >= MIN_C_FOR_ELASTICITY
                        && c_o_aligned.is_finite()
                        && (k - 1.0).abs() > 1e-6 =>
                {
                    (c_o_aligned / c_o).ln() / k.ln()
                }
                _ => f64::NAN,
            };
            Rec {
                key: encode_srgb(lum(t).clamp(0.0, 1.0) as f32) as f64,
                y_o: lum(o),
                y_t: lum(t),
                c_o,
                c_t,
                c_o_aligned,
                elasticity,
                shift: shift_of(o, t),
            }
        })
        .collect()
}

/// 信号门槛的排除原因统计：`(太暗、已剪切、非有限)`。
fn guard_stats(eval: &[DiagSample], ours: fn(&DiagSample) -> [f32; 3]) -> (usize, usize, usize) {
    // 排除原因取「最靠近源头」的一条：非有限 → 太暗 → 已剪切（同一像素可能同时命中，
    // 但它是被排除的，分类只影响报告里的解释，不影响统计）。
    let (mut dark, mut clipped, mut bad) = (0usize, 0usize, 0usize);
    for s in eval {
        let o = ours(s);
        let t = s.theirs_linear;
        let mut reason = 0u8;
        for &x in o.iter().chain(t.iter()) {
            if !x.is_finite() {
                reason = 1;
                break;
            }
            if x < MIN_SIGNAL {
                reason = 2;
                break;
            }
            if x > MAX_SIGNAL {
                reason = 3;
            }
        }
        match reason {
            1 => bad += 1,
            2 => dark += 1,
            3 => clipped += 1,
            _ => {}
        }
    }
    (dark, clipped, bad)
}

/// 一个口径下的整体数字。
struct Overall {
    pass: usize,
    y_o: Option<Dist>,
    y_t: Option<Dist>,
    c_o: Option<Dist>,
    c_t: Option<Dist>,
    /// 我们把明度对齐到参考之后的 C\*。
    c_o_aligned: Option<Dist>,
    /// 逐像素 `C*_ours / C*_ref` 的分布。
    c_ratio: Option<Dist>,
    /// 逐像素 `C*_ours(明度对齐) / C*_ref` 的分布——**扣除亮度差之后**的彩度比。
    c_aligned_ratio: Option<Dist>,
    d_rg: Option<Dist>,
    d_bg: Option<Dist>,
    d_mag: Option<Dist>,
    radial: Option<Dist>,
    tangential: Option<Dist>,
}

fn dist_of(recs: &[Rec], f: fn(&Rec) -> f64) -> Option<Dist> {
    dist(&recs.iter().map(f).collect::<Vec<_>>())
}

fn overall(recs: &[Rec]) -> Overall {
    Overall {
        pass: recs.iter().filter(|r| r.shift.is_some()).count(),
        y_o: dist_of(recs, |r| r.y_o),
        y_t: dist_of(recs, |r| r.y_t),
        c_o: dist_of(recs, |r| r.c_o),
        c_t: dist_of(recs, |r| r.c_t),
        c_o_aligned: dist_of(recs, |r| r.c_o_aligned),
        c_ratio: dist_of(recs, |r| r.c_o / r.c_t),
        c_aligned_ratio: dist_of(recs, |r| r.c_o_aligned / r.c_t),
        d_rg: dist_of(recs, |r| r.shift.map_or(f64::NAN, |s| s.d_rg)),
        d_bg: dist_of(recs, |r| r.shift.map_or(f64::NAN, |s| s.d_bg)),
        d_mag: dist_of(recs, |r| r.shift.map_or(f64::NAN, |s| s.d_mag)),
        radial: dist_of(recs, |r| r.shift.map_or(f64::NAN, |s| s.radial)),
        tangential: dist_of(recs, |r| r.shift.map_or(f64::NAN, |s| s.tangential)),
    }
}

// ---------------------------------------------------------------------------
// 分箱
// ---------------------------------------------------------------------------

/// 一个亮度箱里的样本。
#[derive(Default)]
struct Bin {
    count: usize,
    pass: usize,
    /// 亮度比 `Y_ours / Y_ref`。
    k: Vec<f64>,
    y_t: Vec<f64>,
    c_t: Vec<f64>,
    c_o: Vec<f64>,
    c_ratio: Vec<f64>,
    c_aligned_ratio: Vec<f64>,
    d_rg: Vec<f64>,
    d_bg: Vec<f64>,
}

/// 该箱的比值差是否超出可察觉阈；样本不足或没有通过门槛的样本时无法判定。
fn bin_over(b: &Bin) -> Option<bool> {
    if b.pass < MIN_BIN {
        return None;
    }
    let rg = dist(&b.d_rg)?.med.exp();
    let bg = dist(&b.d_bg)?.med.exp();
    Some((rg - 1.0).abs() > RATIO_EPS || (bg - 1.0).abs() > RATIO_EPS)
}

/// 按参考亮度分箱。
fn bin_all(recs: &[Rec]) -> Vec<Bin> {
    let mut bins: Vec<Bin> = (0..BINS).map(|_| Bin::default()).collect();
    for r in recs {
        let i = ((r.key * BINS as f64) as usize).min(BINS - 1);
        let b = &mut bins[i];
        b.count += 1;
        b.k.push(r.y_o / r.y_t);
        b.y_t.push(r.y_t);
        b.c_t.push(r.c_t);
        b.c_o.push(r.c_o);
        b.c_ratio.push(r.c_o / r.c_t);
        b.c_aligned_ratio.push(r.c_o_aligned / r.c_t);
        if let Some(s) = r.shift {
            b.pass += 1;
            b.d_rg.push(s.d_rg);
            b.d_bg.push(s.d_bg);
        }
    }
    bins
}

// ---------------------------------------------------------------------------
// 逐通道单调映射控制
// ---------------------------------------------------------------------------

/// 三条独立的「ours_c → theirs_c」单调曲线（sRGB 编码域等距取样，编码域内线性插值）。
///
/// # 这一步回答什么
///
/// 尺度不变的比值差只能说「两张图差不只是整体增益」。但**逐通道各自一条单调曲线**这一族
/// （= 现在的 1D 曲线能做的事，也是「每个通道各自偏暗/偏亮」能造出的全部效果）同样会改变
/// 通道比值。把这一族让掉之后还剩多少，才分得开两件事：
///
/// - 让掉之后就塌到 1.000 ⇒ 解码与参考的通道关系**在逐通道单调意义下已被解释**，
///   读数里的「偏淡」与「偏暗」是同一个逐通道差异的两种说法，不需要另立色彩矩阵的修法；
/// - 让掉之后仍停在同一量级 ⇒ 存在**跨通道**的色度错，1D 曲线无论怎么拟合都够不着。
///
/// 注意这是**上界**：三条自由曲线能吸收的东西很多，所以「塌了」只说明差异与逐通道模型
/// 相容，「没塌」才是强结论。
struct MonoMap {
    /// 每通道 `MONO_BINS + 1` 个节点，节点值 = 该编码位置处参考线性值的**中位数**。
    nodes: [Vec<f32>; 3],
}

impl MonoMap {
    /// 用一半样本拟合（每通道独立）。
    fn fit(half: &[DiagSample], ours: fn(&DiagSample) -> [f32; 3]) -> MonoMap {
        let mut nodes: [Vec<f32>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for k in 0..3 {
            let mut buckets: Vec<Vec<f32>> = vec![Vec::new(); MONO_BINS];
            for s in half {
                let o = ours(s);
                let at = ((encode_srgb(o[k]) * MONO_BINS as f32) as usize).min(MONO_BINS - 1);
                buckets[at].push(s.theirs_linear[k]);
            }
            // 箱内中位数；空箱两侧线性插值、两端外推——与 `fit::fit_curve` 同一套做法。
            let mut med: Vec<Option<f32>> = buckets
                .iter_mut()
                .map(|v| {
                    if v.is_empty() {
                        None
                    } else {
                        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                        Some(v[v.len() / 2])
                    }
                })
                .collect();
            let known: Vec<usize> = (0..MONO_BINS).filter(|i| med[*i].is_some()).collect();
            for i in 0..MONO_BINS {
                if med[i].is_some() {
                    continue;
                }
                let lo = known.iter().rev().find(|k2| **k2 < i).copied();
                let hi = known.iter().find(|k2| **k2 > i).copied();
                med[i] = match (lo, hi) {
                    (Some(l), Some(h)) => {
                        let f = (i - l) as f32 / (h - l) as f32;
                        Some(med[l].unwrap_or(0.0) * (1.0 - f) + med[h].unwrap_or(0.0) * f)
                    }
                    (Some(l), None) => med[l],
                    (None, Some(h)) => med[h],
                    (None, None) => Some(0.0),
                };
            }
            // 采样成 MONO_BINS+1 个点，并强制单调不减（非单调只可能来自噪声）。
            let mut n: Vec<f32> = (0..=MONO_BINS)
                .map(|i| med[i.min(MONO_BINS - 1)].unwrap_or(0.0))
                .collect();
            for i in 1..n.len() {
                if n[i] < n[i - 1] {
                    n[i] = n[i - 1];
                }
            }
            nodes[k] = n;
        }
        MonoMap { nodes }
    }

    /// 施加该映射（编码域内线性插值）。
    fn apply(&self, v: [f32; 3]) -> [f32; 3] {
        if v.iter().any(|x| !x.is_finite()) {
            return v;
        }
        let mut out = [0.0_f32; 3];
        for (k, node) in self.nodes.iter().enumerate() {
            let x = encode_srgb(v[k]) * MONO_BINS as f32;
            let i = (x as usize).min(MONO_BINS - 1);
            let f = x - i as f32;
            out[k] = node[i] + (node[i + 1] - node[i]) * f;
        }
        out
    }
}

/// 逐通道单调映射的控制结果。
struct MonoControl {
    /// 参与比较的样本数（评估半区）。
    n: usize,
    /// 映射前 `(R/G, B/G)` 的乘性因子。
    before: (f64, f64),
    /// 映射后同样的因子。
    after: (f64, f64),
}

/// 分半做控制：偶数号样本拟合映射，奇数号样本评估。
///
/// 同一批像素既当学生又当考官会高估模型的解释力，所以宁可只用一半样本评估。
fn mono_control(eval: &[DiagSample], ours: fn(&DiagSample) -> [f32; 3]) -> Option<MonoControl> {
    let fit: Vec<DiagSample> = eval.iter().step_by(2).copied().collect();
    let test: Vec<DiagSample> = eval.iter().skip(1).step_by(2).copied().collect();
    if fit.len() < 1000 || test.is_empty() {
        return None;
    }
    let map = MonoMap::fit(&fit, ours);
    let (mut before_rg, mut before_bg) = (Vec::new(), Vec::new());
    let (mut after_rg, mut after_bg) = (Vec::new(), Vec::new());
    for s in &test {
        let o = ours(s);
        let t = s.theirs_linear;
        // 只有映射前后都过门槛的样本才进对照，否则两组样本不同、没法比。
        if let (Some(b), Some(a)) = (shift_of(o, t), shift_of(map.apply(o), t)) {
            before_rg.push(b.d_rg);
            before_bg.push(b.d_bg);
            after_rg.push(a.d_rg);
            after_bg.push(a.d_bg);
        }
    }
    let f = |d: Option<Dist>| d.map_or(f64::NAN, |x| x.med.exp());
    Some(MonoControl {
        n: before_rg.len(),
        before: (f(dist(&before_rg)), f(dist(&before_bg))),
        after: (f(dist(&after_rg)), f(dist(&after_bg))),
    })
}

// ---------------------------------------------------------------------------
// 打印
// ---------------------------------------------------------------------------

/// C\* 面板：先把「C\* 不是尺度不变量」这件事用同一批像素量出来。
fn print_c_panel(recs: &[Rec], ov: &Overall) {
    let c_t = med_of(&ov.c_t);
    let c_o = med_of(&ov.c_o);
    let c_o_aligned = med_of(&ov.c_o_aligned);
    // 两个口径都给：报告里的「我们 9.92 / 参考 18.57 = 0.53」是**中位之比**；而逐像素比值的
    // 中位数是另一个量（接近中性的像素上比值不稳，会把两者拉开）。混用这两个数会得出
    // 互相矛盾的结论，所以并排打印。
    let med_ratio = c_o / c_t;
    let px_ratio = med_of(&ov.c_ratio);
    let med_ratio_aligned = c_o_aligned / c_t;
    let px_ratio_aligned = med_of(&ov.c_aligned_ratio);
    let elasticity = med_of(&dist_of(recs, |r| r.elasticity));

    println!();
    println!("【一、C* 的口径陷阱：它随亮度变化，不是尺度不变量】");
    println!(
        "  C* 中位：参考 {c_t:.3} · 我们解码 {c_o:.3} → 我们/参考 = {med_ratio:.4}（逐像素比值的中位 {px_ratio:.4}）"
    );
    println!(
        "  明度对齐后：我们 {c_o_aligned:.3} → 我们/参考 = {med_ratio_aligned:.4}（逐像素比值的中位 {px_ratio_aligned:.4}）"
    );
    println!(
        "  亮度中位：我们 {:.5} · 参考 {:.5}（相对差 {:+.1}%）",
        med_of(&ov.y_o),
        med_of(&ov.y_t),
        (med_of(&ov.y_o) / med_of(&ov.y_t) - 1.0) * 100.0
    );
    println!("  C* 对亮度的弹性（中位）= {elasticity:.3}   （1 = 与亮度成正比，1/3 = 立方根区；实测的像素级弹性）");
    if med_ratio.is_finite() && med_ratio_aligned.is_finite() && med_ratio > 0.0 && med_ratio_aligned > 0.0
    {
        let raw = med_ratio.ln().abs();
        if raw < 0.05 {
            println!("  ⇒ 未对齐亮度的差距本身已接近 1（两个方向的误差相互抵消），「解释了多少」的百分比无意义；");
            println!("    要紧的是下一行：明度对齐后反而变成 ×{med_ratio_aligned:.4}。");
        } else {
            let remaining = med_ratio_aligned.ln().abs() / raw;
            println!(
                "  ⇒ 单靠亮度差就能解释 C* 差距的 {:.1}%，余下 {:.1}%（明度对齐后仍在的那部分）才是真的彩度差",
                (1.0 - remaining) * 100.0,
                remaining * 100.0
            );
        }
    }
    println!("  注：弹性由「整体乘 k 使 L* 对齐参考」实测，同一像素两次读数之比，不含任何模型假设。");
}

/// 尺度不变面板：本诊断的判据。
fn print_ratio_panel(stem: &str, n: usize, ov: &Overall, guard: (usize, usize, usize)) {
    let (dark, clipped, bad) = guard;
    println!();
    println!("【二、尺度不变口径（判定的依据）：通道比值的对数差】");
    println!(
        "  通过信号门槛 {}/{} 个（{:.2}%）；排除：太暗 {} 个、已剪切 {} 个、非有限 {} 个",
        ov.pass,
        n,
        ov.pass as f64 / n as f64 * 100.0,
        dark,
        clipped,
        bad
    );
    println!("  口径：ln(R/G)_我们 − ln(R/G)_参考（B/G 同理）。两图各自乘任意正常数，这个数都不变。");
    for (name, d) in [("ln(R/G) 差", &ov.d_rg), ("ln(B/G) 差", &ov.d_bg)] {
        if let Some(d) = d {
            println!(
                "  {name}：{:+.5}（SE {:.5}）· P5 {:+.4} · P95 {:+.4} · 乘性因子 ×{:.4}",
                d.med,
                d.se,
                d.p05,
                d.p95,
                d.med.exp()
            );
        }
    }
    if let (Some(m), Some(rg), Some(bg)) = (&ov.d_mag, &ov.d_rg, &ov.d_bg) {
        println!(
            "  差向量长度 |(Δln R/G, Δln B/G)| 中位 {:.4} · P95 {:.4}；在参考亮度中位处折算 ΔE00 ≈ {:.3}",
            m.med,
            m.p95,
            anchor_delta_e(med_of(&ov.y_t), rg.med, bg.med)
        );
        println!("    （折算口径：把中性灰按测得的 R、B 因子扰动，再用同一把 ΔE00 尺子读数）");
    }
    if let (Some(r), Some(t)) = (&ov.radial, &ov.tangential) {
        println!(
            "  饱和度/色相分解（{} 个有色像素，彩度方向 ≥{MIN_CHROMA_MAG}）：径向 {:+.5}（正 = 我们更饱和）· 切向 {:+.5}",
            r.n, r.med, t.med
        );
        println!("    径向是「整体加浓/减淡」，切向是「色相偏移」；两者都不随整体亮度变化。");
    }
    println!("  判定用 {stem} 的原始解码（ours_linear），未经基准变换。");
}

/// 逐通道单调映射控制面板。
fn print_mono_panel(raw: &Option<MonoControl>, tf: &Option<MonoControl>) {
    println!();
    println!("【三、逐通道单调映射控制：把「每通道各自偏暗/偏亮」这一族让掉之后还剩多少】");
    println!("  做法：sRGB 编码域 128 箱上，用偶数号样本拟合 3 条独立单调曲线 ours_c → theirs_c");
    println!("        （箱内取参考值中位、空箱插值、强制单调不减），在奇数号样本上复测通道比值差。");
    println!("        三条自由单调曲线 = 「曝光 / 每通道各自的曲线 / 白平衡」全都允许不同的**上界**。");
    if let Some(m) = raw {
        println!(
            "  原始解码口径（n={}）：映射前 R/G ×{:.4} · B/G ×{:.4} → 映射后 R/G ×{:.4} · B/G ×{:.4}",
            m.n, m.before.0, m.before.1, m.after.0, m.after.1
        );
    }
    if let Some(m) = tf {
        println!(
            "  变换后口径  （n={}）：映射前 R/G ×{:.4} · B/G ×{:.4} → 映射后 R/G ×{:.4} · B/G ×{:.4}",
            m.n, m.before.0, m.before.1, m.after.0, m.after.1
        );
    }
    println!("  读法：映射后塌到 1.0000 附近 ⇒ 这个差**与逐通道模型相容**，它与亮度差同源，1D 曲线");
    println!("        就是对症的修法；仍停在同一量级 ⇒ 存在跨通道（矩阵/饱和度）的真实色度错。");
}

/// 分箱表：比值差是否在各亮度段持续存在。
fn print_bin_table(bins: &[Bin], n: usize) {
    println!();
    println!("【四、按参考亮度分箱：差是「处处都在」还是「只在某一段」】");
    println!("bin  enc_lo  enc_hi    count  pass%  share%   k_o/ref   Yt_med  C*_ref  C*_our  C*rat  C*match     RGx      BGx  dE_ratio");
    for (i, b) in bins.iter().enumerate() {
        let lo = i as f64 / BINS as f64;
        let hi = (i + 1) as f64 / BINS as f64;
        let rg = dist(&b.d_rg).map(|d| d.med.exp());
        let bg = dist(&b.d_bg).map(|d| d.med.exp());
        let y_med = med_of(&dist(&b.y_t));
        let de = match (rg, bg) {
            (Some(x), Some(y)) => Some(anchor_delta_e(y_med, x.ln(), y.ln())),
            _ => None,
        };
        let flag = match bin_over(b) {
            Some(true) => " *",
            Some(false) => "",
            None => " -",
        };
        println!(
            "{i:>3}  {lo:>6.4}  {hi:>6.4}  {:>7}  {:>5.1}%  {:>5.2}%  {}  {}  {}  {}  {}  {}  {}  {}  {}{flag}",
            b.count,
            b.pass as f64 / b.count.max(1) as f64 * 100.0,
            b.count as f64 / n as f64 * 100.0,
            col(dist(&b.k).map(|d| d.med), 3),
            col(Some(y_med), 5),
            col(dist(&b.c_t).map(|d| d.med), 3),
            col(dist(&b.c_o).map(|d| d.med), 3),
            col(dist(&b.c_ratio).map(|d| d.med), 3),
            col(dist(&b.c_aligned_ratio).map(|d| d.med), 3),
            col(rg, 4),
            col(bg, 4),
            col(de, 3),
        );
    }
    println!("  列义：count 样本数 · pass% 通过信号门槛的比例（深暗箱必然低）· share% 占本次评估集");
    println!("        k_o/ref 该箱亮度比中位（1 = 亮度已对齐）· Yt_med 参考亮度中位（线性域）");
    println!("        C*rat = C*_our / C*_ref 逐像素比值的中位（**未对齐亮度**）");
    println!("        C*match = 明度对齐后同样比值的中位（**扣除亮度差之后**）");
    println!("        RGx / BGx = exp(ln 比值差的中位)，1.0000 = 与参考一致；dE_ratio = 该箱比值差折算的 ΔE00");
    println!(
        "        * = 该箱比值差超 {:.0}% 可察觉阈；- = 通过门槛的样本不足 {MIN_BIN} 个，不参与判定",
        RATIO_EPS * 100.0
    );
}

/// 单图结论。
fn print_verdict(
    stem: &str,
    bins: &[Bin],
    ov: &Overall,
    mono: &Option<MonoControl>,
) -> (&'static str, usize, usize) {
    let flags: Vec<Option<bool>> = bins.iter().map(bin_over).collect();
    let valid = flags.iter().filter(|f| f.is_some()).count();
    let over = flags.iter().filter(|f| **f == Some(true)).count();

    let verdict = if valid == 0 {
        "样本不足，无法判定"
    } else if over == valid {
        "确实是色度错（整体增益解释不掉）"
    } else if over == 0 {
        "只是亮度差造成的假象"
    } else {
        "两者都有（部分亮度段超阈）"
    };

    println!();
    println!("── 结论[{stem}] ──");
    println!(
        "  1) C*（中位之比，与报告的 9.92/18.57 同口径）：未对齐亮度 ×{:.3}；**明度对齐后 ×{:.3}**；逐像素比值中位 ×{:.3} → ×{:.3}",
        med_of(&ov.c_o) / med_of(&ov.c_t),
        med_of(&ov.c_o_aligned) / med_of(&ov.c_t),
        med_of(&ov.c_ratio),
        med_of(&ov.c_aligned_ratio)
    );
    if let (Some(rg), Some(bg)) = (&ov.d_rg, &ov.d_bg) {
        println!(
            "  2) 尺度不变比值差：R/G ×{:.4}（SE {:.5}，n={}）· B/G ×{:.4}（SE {:.5}，n={}）",
            rg.med.exp(),
            rg.se,
            rg.n,
            bg.med.exp(),
            bg.se,
            bg.n
        );
    }
    println!("  3) 是否各亮度段都持续存在：可判定的 {valid} 个箱中 {over} 个超阈。");
    if over > 0 && over < valid {
        let miss: Vec<String> = flags
            .iter()
            .enumerate()
            .filter(|(_, f)| **f == Some(false))
            .map(|(i, _)| i.to_string())
            .collect();
        println!("     未超阈的箱：{}", miss.join(", "));
    }
    if let Some(m) = mono {
        let still = if (m.after.0 - 1.0).abs().max((m.after.1 - 1.0).abs()) > RATIO_EPS {
            "仍超阈"
        } else {
            "已落到阈内"
        };
        println!(
            "  4) 逐通道单调映射后：R/G ×{:.4} · B/G ×{:.4}（映射前 ×{:.4} / ×{:.4}）——{still}",
            m.after.0, m.after.1, m.before.0, m.before.1
        );
    }
    println!("  ⇒ 判定：**{verdict}**");
    if let Some(m) = mono {
        let collapsed =
            (m.after.0 - 1.0).abs() <= RATIO_EPS && (m.after.1 - 1.0).abs() <= RATIO_EPS;
        if collapsed {
            println!(
                "     附注：上面这个差**被逐通道单调映射吃掉了**（×{:.4} / ×{:.4}）——它与「每通道各自偏暗/偏亮」",
                m.after.0, m.after.1
            );
            println!("           是同一个来源，不需要跨通道的色彩矩阵修法；要跨通道手段的是变换后剩的那点残差（见【五】）。");
        } else {
            println!(
                "     附注：逐通道单调映射吸收不掉（映射后仍 ×{:.4} / ×{:.4}）⇒ 这里有跨通道的真实色度错。",
                m.after.0, m.after.1
            );
        }
    }
    (verdict, over, valid)
}

/// 跨图汇总。
struct ImageReport {
    stem: &'static str,
    size: (usize, usize),
    n: usize,
    verdict: &'static str,
    over: usize,
    valid: usize,
    c_ratio: f64,
    c_aligned: f64,
    rg: f64,
    bg: f64,
    /// 逐通道单调映射之后的 `(R/G, B/G)` 因子。
    mono: (f64, f64),
}

fn print_cross(reports: &[ImageReport]) {
    println!();
    println!("【六、跨图汇总】");
    println!("stem           size               n   C*rat  C*match     RGx      BGx   M-RGx    M-BGx  bins  verdict");
    for r in reports {
        println!(
            "{:<10}  {:<14}  {:>7}  {:>6.4}  {:>7.4}  {:>6.4}  {:>6.4}  {:>6.4}  {:>7.4}  {:>2}/{:<2}  {}",
            r.stem,
            format!("{}×{}", r.size.0, r.size.1),
            r.n,
            r.c_ratio,
            r.c_aligned,
            r.rg,
            r.bg,
            r.mono.0,
            r.mono.1,
            r.over,
            r.valid,
            r.verdict
        );
    }
    println!("C*rat = 未对齐亮度的彩度比；C*match = 明度对齐后的彩度比；bins = 超阈箱数/可判定箱数。");
    println!("M-RGx / M-BGx = 逐通道单调映射之后的通道比值因子（1.0000 = 该差已被逐通道模型解释掉）。");
}

// ---------------------------------------------------------------------------
// 主流程
// ---------------------------------------------------------------------------

/// 跑一张图：打印三个面板 + 结论 + 变换后对照，返回汇总行。
fn run_one(stem: &'static str, fx: &FittedFixture) -> ImageReport {
    let n = fx.eval.len();

    // 判定口径：原始解码。
    let recs = recs_of(&fx.eval, |s| s.ours_linear);
    let ov = overall(&recs);
    let bins = bin_all(&recs);
    assert_eq!(
        bins.iter().map(|b| b.count).sum::<usize>(),
        n,
        "分箱必须覆盖全部样本（漏样本会让所有占比失去意义）"
    );
    if let Some(d) = &ov.d_rg {
        assert!(d.med.is_finite(), "比值差中位数不应是 NaN");
    }

    println!();
    println!("════ 图 {stem}（{}×{}）════", fx.size.0, fx.size.1);
    println!(
        "样本：留出评估集 {n} 个 · 拟合集 {} 个（两者用互不相交的像素子集）",
        fx.fit_count
    );

    print_c_panel(&recs, &ov);
    print_ratio_panel(stem, n, &ov, guard_stats(&fx.eval, |s| s.ours_linear));
    // 变换后的口径：现在管线里跑的曲线 + LUT 把这个差修掉了多少。
    let recs_t = recs_of(&fx.eval, |s| s.ours_transformed);
    let ov_t = overall(&recs_t);
    let mono_raw = mono_control(&fx.eval, |s| s.ours_linear);
    let mono_tf = mono_control(&fx.eval, |s| s.ours_transformed);
    print_mono_panel(&mono_raw, &mono_tf);
    print_bin_table(&bins, n);
    let (verdict, over, valid) = print_verdict(stem, &bins, &ov, &mono_raw);

    let (dark_t, clip_t, _) = guard_stats(&fx.eval, |s| s.ours_transformed);
    let mut de: Vec<f64> = fx.eval.iter().map(|s| s.delta_e()).collect();
    let (de_med, de_p95, de_max) = quantiles(&mut de);
    let pairs: Vec<([f32; 3], [f32; 3])> = fx
        .eval
        .iter()
        .map(|s| (s.ours_transformed, s.theirs_linear))
        .collect();
    let sum = deltae::summarize_prophoto(&pairs);

    println!();
    println!("【五、变换后口径的对照（现在管线里的曲线 + LUT）】");
    println!(
        "  ΔE00：中位 {de_med:.3} · P95 {de_p95:.3} · 最大 {de_max:.3}（目标 ≤1.0 / ≤3.0）；明度合计占比 {:.1}%",
        sum.total_share[0] * 100.0
    );
    println!(
        "  变换后 C*：未对齐 ×{:.4} → 明度对齐后 ×{:.4}",
        med_of(&ov_t.c_ratio),
        med_of(&ov_t.c_aligned_ratio)
    );
    if let (Some(rg), Some(bg)) = (&ov_t.d_rg, &ov_t.d_bg) {
        println!(
            "  变换后尺度不变比值差：R/G ×{:.4}（SE {:.5}）· B/G ×{:.4}（SE {:.5}）；通过门槛 {}/{}（太暗 {}、已剪切 {}）",
            rg.med.exp(),
            rg.se,
            bg.med.exp(),
            bg.se,
            ov_t.pass,
            n,
            dark_t,
            clip_t
        );
    }
    // 变换前后比值差的变化，直接回答「基准变换修掉了多少」。
    if let (Some(a), Some(b)) = (&ov.d_rg, &ov_t.d_rg) {
        let before = a.med.exp();
        let after = b.med.exp();
        let kept = if before.abs() > 1e-9 {
            let ratio = (after - 1.0) / (before - 1.0);
            if ratio < 0.0 {
                "已越过 1.0（符号翻转，原残差被修掉且略过头）".to_string()
            } else {
                format!("{:.0}%", ratio * 100.0)
            }
        } else {
            "—".to_string()
        };
        println!("  → R/G 因子由变换前的 ×{before:.4} 变为 ×{after:.4}（残差保留 {kept}）");
    }

    ImageReport {
        stem,
        size: fx.size,
        n,
        verdict,
        over,
        valid,
        c_ratio: med_of(&ov.c_ratio),
        c_aligned: med_of(&ov.c_aligned_ratio),
        rg: med_of(&ov.d_rg).exp(),
        bg: med_of(&ov.d_bg).exp(),
        mono: mono_raw.map_or((f64::NAN, f64::NAN), |m| m.after),
    }
}

/// 尺度不变的色度比对：我们究竟「偏淡」，还是「只是偏暗」。
///
/// 只做诊断：断言的都是内部自洽性（分箱覆盖全部样本、统计量非 NaN），不断言达标与否
/// ——当前已知未达标，把结论写死成断言只会让用例变红而不会带来信息。
#[test]
fn scale_invariant_chroma_comparison() {
    println!();
    println!("======================================================================");
    println!("诊断 D：尺度不变色度比对——「确实是色度错」还是「只是亮度差」");
    println!("======================================================================");
    println!("判据：ln(R/G)、ln(B/G) 的对数差。两张图各自乘任意正常数，它们一个都不动，");
    println!("      所以这里读出的差异不可能是「谁更亮」造成的。判定一律用 ours_linear。");

    let mut reports = Vec::new();
    for stem in STEMS {
        let Some(fx) = fitted_fixture_for(stem, true) else {
            continue;
        };
        if fx.eval.is_empty() {
            println!("跳过 {stem}：评估集为空");
            continue;
        }
        reports.push(run_one(stem, &fx));
    }

    if reports.is_empty() {
        println!();
        println!("三张图都没有可用样本——本测试按仓库惯例跳过而非失败。");
        return;
    }
    print_cross(&reports);
}
