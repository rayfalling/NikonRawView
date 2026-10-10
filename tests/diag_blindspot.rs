//! 诊断 F：**色度盲区量化**——近乎中性的参考样本对「乘性饱和度偏差」有多不敏感。
//!
//! # 为什么需要这一条
//!
//! 团队诊断把疑似缺陷锁定为**乘性饱和度偏差**（我们原始解码的 C* 只有参考的
//! 0.53~0.63 倍）。而原先用于验收的 15 张竖构图样本是 42 秒内连拍的同一近中性场景，
//! 它们把这个偏差**掩盖**了：`DSC_0001` 的 ΔE00 中位数只有 1.698，而两张横构图彩色
//! 样本是 3.543 / 4.767。**用一组几乎没色度的样本去验收色度，等于没验。**
//!
//! 这条诊断把"掩盖"量化成一个倍数：
//!
//! 1. 量化 `DSC_0001` 的参考 C* 分布（中位 / P90 / P99 / C\*>10 与 C\*>25 的占比）；
//! 2. **注入实验**：人为把我们的原始解码值按固定倍率 `k` 缩向灰度
//!    （`v' = lum + k·(v − lum)`，`lum` 取 Rec.709 权重），对每个 `k`
//!    **重新拟合曲线与 3D LUT**，再看评估集的 ΔE00 中位数；
//! 3. 对横构图的 `DSC_0141` 做同一实验作对照；
//! 4. 两边"要让 ΔE00 中位数恶化 1.0 所需的饱和度偏差"之比，就是**盲区倍数**；
//! 5. 另给一条**不重拟合**的对照列，量出"重拟合到底补偿掉了多少"——任务书警告的那个坑，
//!    这里把它变成数字；
//! 6. 顺带把"我们解码 C* 只有参考的 0.53~0.63 倍"这条前提核一遍：C* 倍率会被亮度差
//!    **冒充**（Lab 的 a\*/b\* 随亮度立方根缩放），所以同时给出扣掉亮度差后的残余倍率，
//!    以及残余倍率随参考 C* 的走势。
//!
//! # 本机实测（DSC_0001 / DSC_0141）
//!
//! * **盲区倍数 0.9×**（主口径，重拟合）、0.6×（近基线灵敏度）、1.2×（不重拟合对照）
//!   —— 三者都在 1 附近，**"那 15 张样本对色度缺陷很钝"这个预期不成立**；
//! * `DSC_0001` 并不是"几乎没色度"：参考 C* 中位 **6.62**、P90 11.79、P99 32.66，
//!   C\*>10 占 **26.55%**；但 C\*>15 只占 1.53%，而 `DSC_0141` 是 55.51%；
//! * 扣掉亮度差后，`DSC_0001` 的残余色度倍率 **1.166**（色度不缺，反而略过），
//!   `DSC_0141` 为 **0.519**（真亏）；且 `DSC_0141` 在 C\* 10–15 这一档就已经是
//!   0.453，而同一档 `DSC_0001` 是 1.169 ⇒ 这个亏欠**不是**由像素彩度高低触发的，
//!   更像**逐图**性质。**"盲"来自那 15 张覆盖不到出问题的画面，而不是度量不灵敏。**
//! * 重拟合的补偿量：`DSC_0001` ≈ 0（−2%~+1%），`DSC_0141` 4%~18%（随偏差增大）。
//!   ⇒ "必须重拟合"这条要求是对的，但它的作用主要落在彩色图上。
//!
//! # 两处最容易做错的地方
//!
//! * **注入后必须重新拟合。** 基准变换（曲线 + LUT）本身会补偿掉一部分系统性色度
//!   偏差；沿用原变换不重拟合，测到的是"原变换直接套在偏移数据上"的误差，量级完全
//!   不同。本文件每个 `k` 都从头 `fit::fit_curve` + `fit::fit_lut`，不重拟合只作对照列。
//! * **样本只解码一次。** 一张 45 MP 的 NEF 解码约 20 秒、参考 TIF 有 260 MB，
//!   19 个 `k` × 2 张图若各自解码一遍要跑十几分钟。这里用
//!   `fitted_fixture_for(stem, false)`（`do_fit = false`，不自动拟合）取一次样本，
//!   所有 `k` 复用同一批 `ours_linear` / `theirs_linear`。
//!
//! # 口径
//!
//! * 只走基准变换的**曲线 + LUT** 两段（`BaseTransform::apply` 就是先 `forward_curve`
//!   再 `apply_lut`），配方那几步不参与——这与 `tests/common/mod.rs` 的脚手架、以及
//!   正式基线 1.698 / 3.543 的口径一致。
//! * 达标线 1.0 是**绝对**阈值，而两张图的基线本来就已在它之上，所以"中位数超过
//!   1.0"这个字面口径本身是退化的（偏差为 0 就已经"超过"）。本诊断因此给出三套口径，
//!   主口径是**相对基线再高 1.0**——只有它给出可比的偏差量。
//! * 本文件只断言**测量有效性**（注入确实发生且保亮度、拟合收敛、扫描把穿越点夹住、
//!   基线没跑偏），不断言"盲区一定存在"——那是本诊断要测的量，不是前提。
//!
//! 样本缺席时跳过而非失败（见 `common::skip`）。

mod common;

use common::{fitted_fixture_for, quantiles, DiagSample};
use nikonrawview::deltae;
use nikonrawview::fit::{self, Sample};
use nikonrawview::transform::BaseTransform;

// ---------------------------------------------------------------------------
// 口径常量
// ---------------------------------------------------------------------------

/// 拟合口径：128 箱曲线 + 17³ LUT——与 `tests/common/mod.rs` 的脚手架一致，
/// 好让 `k = 1.0` 那一行能与正式基线直接对照。
const CURVE_BINS: usize = 128;
const LUT_EDGE: usize = 17;

/// 正式达标线：ΔE00 中位数 ≤ 1.0。
const TARGET_MEDIAN: f64 = 1.0;

/// 参与对照的两张图：竖构图近中性（代表原来那 15 张）与横构图彩色。
const STEMS: [(&str, &str); 2] = [
    ("DSC_0001", "竖构图近中性 —— 原来那 15 张的代表"),
    ("DSC_0141", "横构图彩色 —— 对照"),
];

/// 已知的正式基线（全密度拟合口径的 ΔE00 中位数），用来粗查本诊断的口径没跑偏。
const PUBLISHED_MEDIAN: [(&str, f64); 2] = [("DSC_0001", 1.698), ("DSC_0141", 3.543)];

/// 注入扫描的色度保留倍率。
///
/// 任务书要求的 `{1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60}` **全部在内**；
/// 另在 0.9 附近加密（穿越点大概率落在这里，需要插值精度），并向下延伸到 0.00
/// （全灰，必定到达上界）以确保把两条穿越点都夹住。
const K_SWEEP: [f32; 19] = [
    1.00, 0.98, 0.95, 0.92, 0.90, 0.87, 0.85, 0.82, 0.80, 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.30,
    0.20, 0.10, 0.00,
];

// ---------------------------------------------------------------------------
// 注入
// ---------------------------------------------------------------------------

/// Rec.709 加权和（与 `DiagSample::lum_ours` 同一约定）。
///
/// 权重向量在 ProPhoto 线性值上并不是严格的物理亮度，这里沿用它只是为了给"色度"
/// 定一条**固定的投影轴**，与仓库里其它诊断保持一致。
fn lum709(v: [f32; 3]) -> f32 {
    0.2126 * v[0] + 0.7152 * v[1] + 0.0722 * v[2]
}

/// 保亮度的饱和度缩放：`v' = lum + k·(v − lum)`。
///
/// * `k = 1.0` 原状，`k = 0.0` 全灰；`k ∈ [0,1]` 时是 `lum` 与 `v` 的凸组合，
///   因此**不会**产生越界值或负值，无需裁切。
/// * 权重之和为 1，所以 `lum(v') = lum(v)` **精确成立**——注入只动色度、不动亮度。
fn sat_scale(v: [f32; 3], k: f32) -> [f32; 3] {
    let lum = lum709(v);
    [
        lum + k * (v[0] - lum),
        lum + k * (v[1] - lum),
        lum + k * (v[2] - lum),
    ]
}

/// 注入的平均强度：所有通道上 `|v' − v|` 的均值。用来确认"注入真的发生了"。
fn injection_strength(samples: &[Sample], k: f32) -> f64 {
    if samples.is_empty() {
        return 0.0;
    }
    let mut acc = 0.0f64;
    for s in samples {
        let v = sat_scale(s.ours, k);
        acc += v
            .iter()
            .zip(s.ours.iter())
            .map(|(a, b)| (a - b).abs() as f64)
            .sum::<f64>();
    }
    acc / (samples.len() * 3) as f64
}

// ---------------------------------------------------------------------------
// 色度统计
// ---------------------------------------------------------------------------

/// Lab 彩度 `C* = √(a² + b²)`（ProPhoto 线性 → CIELAB D50）。
fn c_star(v: [f32; 3]) -> f64 {
    let lab = deltae::lab_from_prophoto_linear(v);
    ((lab[1] as f64).powi(2) + (lab[2] as f64).powi(2)).sqrt()
}

/// nearest-rank 百分位（口径与 `src/deltae.rs` 一致：`ceil(q·n) − 1`），`values` 须已升序。
fn at_rank(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let rank = (q * sorted.len() as f64).ceil().max(1.0) as usize;
    sorted[rank.saturating_sub(1).min(sorted.len() - 1)]
}

/// 已升序切片中大于 `x` 的元素个数。
fn count_gt(sorted: &[f64], x: f64) -> usize {
    sorted.len() - sorted.partition_point(|v| *v <= x)
}

/// 中位数（升序口径同 [`at_rank`]）。
fn median_of(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    at_rank(values, 0.5)
}

/// 参考 C* 的分段——用来判断"色度亏欠是乘性常数，还是随彩度升高才出现"。
///
/// 这一区分直接决定结论往哪走：若亏欠只在高彩度段出现，那么近中性样本测不出来并不是
/// **度量钝**，而是**画面里根本没有那段内容**。
const C_BINS: [(f64, f64); 6] = [
    (10.0, 15.0),
    (15.0, 20.0),
    (20.0, 25.0),
    (25.0, 30.0),
    (30.0, 35.0),
    (35.0, 100.0),
];

/// 一个参考 C* 区间上的统计。
struct ChromaBin {
    lo: f64,
    hi: f64,
    n: usize,
    /// 扣亮度后的残余色度倍率中位（空段为 NaN）。
    detoned: f64,
}

/// 参考画面的色度内容——量化"这张图到底有多中性"。
struct ChromaStats {
    n: usize,
    median: f64,
    p90: f64,
    p99: f64,
    /// `C* > 10`（明显有色）的占比。
    gt10: f64,
    /// `C* > 25`（高彩度）的占比。
    gt25: f64,
    /// `C* > 10` 子集的样本数（倍率的可信度取决于它）。
    n_chromatic: usize,
    /// 该子集上「我们原始解码 C* ÷ 参考 C*」的中位比值。
    ours_ratio: f64,
    /// 同一子集上的 L* 中位（参考 / 我们原始解码）。
    ///
    /// **必须和 C* 倍率一起看**：Lab 的 a\*、b\* 大致随亮度的立方根缩放，亮度差会
    /// **冒充**色度差。这一对数就是用来判断"色度只有参考的百分之多少"里有多少其实是
    /// 亮度造成的（尺规不变的口径由 `diag_chromaticity` 负责）。
    l_theirs: f64,
    l_ours: f64,
    /// 扣掉亮度差之后**残余**的色度倍率——粗算，逐像素扣一阶项。
    ///
    /// `= (C*_ours / C*_theirs) ÷ (Y_ours / Y_theirs)^(1/3)`，`Y = ((L*+16)/116)³`。
    /// ≈1 表示"色度没缺，缺的只是亮度"；明显 <1 才是真正的饱和度损失。
    detoned_ratio: f64,
    n_detoned: usize,
    /// 残余色度倍率随参考 C* 的走势。
    bins: Vec<ChromaBin>,
}

fn chroma_stats(samples: &[DiagSample]) -> ChromaStats {
    // L* → 相对亮度（CIELAB 的定义式，只对 L* > 8 的三次根段有效，故下面设了门槛）。
    let y_of = |l: f64| ((l + 16.0) / 116.0).powi(3);
    let mut refs: Vec<f64> = Vec::with_capacity(samples.len());
    let mut ratios: Vec<f64> = Vec::new();
    let mut l_theirs: Vec<f64> = Vec::new();
    let mut l_ours: Vec<f64> = Vec::new();
    let mut detoned: Vec<f64> = Vec::new();
    let mut bin_vals: Vec<Vec<f64>> = vec![Vec::new(); C_BINS.len()];
    for s in samples {
        let lt = deltae::lab_from_prophoto_linear(s.theirs_linear);
        let t = ((lt[1] as f64).powi(2) + (lt[2] as f64).powi(2)).sqrt();
        refs.push(t);
        if t > 10.0 {
            let lo = deltae::lab_from_prophoto_linear(s.ours_linear);
            let o = ((lo[1] as f64).powi(2) + (lo[2] as f64).powi(2)).sqrt();
            ratios.push(o / t);
            l_theirs.push(lt[0] as f64);
            l_ours.push(lo[0] as f64);
            if lt[0] > 8.0 && lo[0] > 8.0 {
                let s_scale = (y_of(lo[0] as f64) / y_of(lt[0] as f64)).cbrt();
                if s_scale > 1e-6 {
                    let r = o / t / s_scale;
                    detoned.push(r);
                    if let Some(i) = C_BINS.iter().position(|(a, b)| t >= *a && t < *b) {
                        bin_vals[i].push(r);
                    }
                }
            }
        }
    }
    let n = refs.len();
    refs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let denom = n.max(1) as f64;
    let bins: Vec<ChromaBin> = C_BINS
        .iter()
        .zip(bin_vals.iter_mut())
        .map(|(&(lo, hi), v)| ChromaBin {
            lo,
            hi,
            n: v.len(),
            detoned: median_of(v),
        })
        .collect();
    ChromaStats {
        n,
        median: at_rank(&refs, 0.5),
        p90: at_rank(&refs, 0.9),
        p99: at_rank(&refs, 0.99),
        gt10: count_gt(&refs, 10.0) as f64 / denom,
        gt25: count_gt(&refs, 25.0) as f64 / denom,
        n_chromatic: ratios.len(),
        ours_ratio: median_of(&mut ratios),
        l_theirs: median_of(&mut l_theirs),
        l_ours: median_of(&mut l_ours),
        n_detoned: detoned.len(),
        detoned_ratio: median_of(&mut detoned),
        bins,
    }
}

// ---------------------------------------------------------------------------
// 扫描结果
// ---------------------------------------------------------------------------

/// 一个 `k` 的扫描结果。
struct Row {
    k: f32,
    /// **重新拟合**后的 ΔE00 中位数——本诊断的主口径。
    median: f64,
    /// 不重拟合（拿 `k = 1.0` 那份"出厂标定"直接套）的 ΔE00 中位数——对照。
    median_frozen: f64,
    p95: f64,
    max: f64,
    /// 注入后的**原始解码** C* 相对参考的倍率（`C* > 10` 子集上的中位比值）。
    chroma_ratio: f64,
}

impl Row {
    /// 饱和度偏差 `d = 1 − k`。
    fn deviation(&self) -> f64 {
        1.0 - self.k as f64
    }
}

/// 一张图的完整诊断结果。
struct StemReport {
    stem: &'static str,
    role: &'static str,
    size: (usize, usize),
    fit_n: usize,
    eval_n: usize,
    chroma: ChromaStats,
    /// `k = 0.80` 时的平均注入强度（`|Δ通道|` 均值）。
    strength_080: f64,
    rows: Vec<Row>,
}

/// 把 `samples` 按注入倍率 `k` 缩放色度后**重新拟合**一份基准变换。
///
/// 这是本诊断的核心动作：每个 `k` 都从头拟合，拟合才有机会去补偿注入的偏差。
fn fit_at(samples: &[Sample], k: f32, stem: &str) -> BaseTransform {
    let injected: Vec<Sample> = samples
        .iter()
        .map(|s| Sample {
            ours: sat_scale(s.ours, k),
            theirs: s.theirs,
        })
        .collect();
    let curve = fit::fit_curve(&injected, CURVE_BINS);
    let lut = fit::fit_lut(&injected, &curve, LUT_EDGE);
    BaseTransform::from_fit(
        0x0000,
        "色度盲区注入实验",
        vec![format!("simple/{stem}.NEF + simple/{stem}.TIF")],
        curve,
        LUT_EDGE,
        lut,
    )
}

/// 跑一张图的注入扫描。样本缺席返回 `None`（调用方跳过）。
fn run_stem(stem: &str, role: &'static str) -> Option<StemReport> {
    // do_fit = false：这里只要样本，拟合留给下面每个 k 自己做。
    let fx = fitted_fixture_for(stem, false)?;

    // 留出集再对半劈：偶数下标拟合、奇数下标评估。两半在空间上交错（同一行相邻列
    // 会分属两边），既不会拿同一批像素既训又考，也不会让两半落到画面不同区域。
    let mut fit_side: Vec<Sample> = Vec::with_capacity(fx.eval.len() / 2 + 1);
    let mut eval_side: Vec<Sample> = Vec::with_capacity(fx.eval.len() / 2 + 1);
    for (i, s) in fx.eval.iter().enumerate() {
        let smp = Sample {
            ours: s.ours_linear,
            theirs: s.theirs_linear,
        };
        if i % 2 == 0 {
            fit_side.push(smp);
        } else {
            eval_side.push(smp);
        }
    }
    assert!(
        !fit_side.is_empty() && !eval_side.is_empty(),
        "{stem}：样本不足，无法划分拟合/评估集"
    );

    let chroma = chroma_stats(&fx.eval);
    let ref_cs: Vec<f64> = eval_side.iter().map(|s| c_star(s.theirs)).collect();

    // 注入必须"只动色度、不动亮度"——两条断言钉住投影轴没写错。
    let probe = fit_side[0];
    let lum_shift = (lum709(sat_scale(probe.ours, 0.6)) - lum709(probe.ours)).abs();
    assert!(
        lum_shift < 1e-5,
        "{stem}：注入改变了 Rec.709 亮度（Δ={lum_shift:.3e}），投影轴写错了"
    );
    let strength_080 = injection_strength(&fit_side, 0.8);
    assert!(
        strength_080 > 1e-4,
        "{stem}：k=0.80 的注入没有改变任何样本（平均 |Δ通道| = {strength_080:.3e}）"
    );

    // "出厂标定"：按原状（k = 1.0）拟合出来的那一份。对照用——把它直接套到注入后的
    // 数据上，就得到"不重拟合"口径，用来量出重拟合到底补偿掉了多少。
    let frozen = fit_at(&fit_side, 1.0, stem);

    // ---- 逐个 k：注入 → **重新拟合** → 评估 ----
    let mut rows = Vec::with_capacity(K_SWEEP.len());
    for &k in K_SWEEP.iter() {
        let t = fit_at(&fit_side, k, stem);

        let mut des: Vec<f64> = Vec::with_capacity(eval_side.len());
        let mut des_frozen: Vec<f64> = Vec::with_capacity(eval_side.len());
        let mut ratios: Vec<f64> = Vec::with_capacity(chroma.n_chromatic);
        for (s, ref_c) in eval_side.iter().zip(ref_cs.iter()) {
            let injected_eval = sat_scale(s.ours, k);
            des.push(deltae::delta_e_prophoto(t.apply(injected_eval), s.theirs));
            des_frozen.push(deltae::delta_e_prophoto(
                frozen.apply(injected_eval),
                s.theirs,
            ));
            if *ref_c > 10.0 {
                ratios.push(c_star(injected_eval) / ref_c);
            }
        }
        let (median, p95, max) = quantiles(&mut des);
        rows.push(Row {
            k,
            median,
            median_frozen: median_of(&mut des_frozen),
            p95,
            max,
            chroma_ratio: median_of(&mut ratios),
        });
    }

    assert!(
        (rows[0].k - 1.0).abs() < 1e-6,
        "{stem}：扫描表必须以 k = 1.0（基线）开头"
    );
    // 一致性护栏：k = 1.0 那一行，"出厂标定"与当轮重新拟合本就是同一次拟合，
    // 两列必须逐位一致。若不等，说明对照列读错了变换（本诊断最容易悄悄错的地方）。
    assert!(
        (rows[0].median - rows[0].median_frozen).abs() < 1e-12,
        "{stem}：k=1.0 处重拟合与不重拟合两列不一致（{:.6} vs {:.6}），对照口径写错了",
        rows[0].median,
        rows[0].median_frozen
    );
    for r in &rows {
        assert!(
            r.median.is_finite() && r.p95.is_finite(),
            "{stem}：k={} 的 ΔE00 不是有限值",
            r.k
        );
    }

    // 粗查口径：k = 1.0 那一行应该贴近正式基线（全密度拟合）。
    if let Some(pub_median) = PUBLISHED_MEDIAN
        .iter()
        .find(|(s, _)| *s == stem)
        .map(|(_, m)| *m)
    {
        let diff = rows[0].median - pub_median;
        assert!(
            diff.abs() < 1.5,
            "{stem}：k=1.0 基线 {:.3} 与正式基线 {pub_median:.3} 差 {diff:+.3}，口径跑偏了",
            rows[0].median
        );
    }

    Some(StemReport {
        stem: STEMS
            .iter()
            .find(|(s, _)| *s == stem)
            .map_or("?", |(s, _)| *s),
        role,
        size: fx.size,
        fit_n: fit_side.len(),
        eval_n: eval_side.len(),
        chroma,
        strength_080,
        rows,
    })
}

// ---------------------------------------------------------------------------
// 穿越点
// ---------------------------------------------------------------------------

/// 「ΔE00 中位数达到某阈值」所需的饱和度偏差。
#[derive(Debug, Clone, Copy)]
enum Crossing {
    /// 基线（k = 1.0）已在阈值之上——该阈值口径退化，给不出"需要多大偏差"。
    AlreadyAbove,
    /// 首次达到阈值；`k_lo`（未达）与 `k_hi`（已达）是夹住它的两行，`k_lo > k_hi`。
    At { d: f64, k_lo: f32, k_hi: f32 },
    /// 扫到最大偏差仍未达到。
    Beyond { d_max: f64 },
}

/// 在扫描表上求穿越点（表须按 `k` 降序、即 `d` 升序排列，线性插值）。
///
/// `of` 选定读哪一列（`|r| r.median` 为主口径，`|r| r.median_frozen` 为不重拟合的对照）。
fn crossing(rows: &[Row], target: f64, of: fn(&Row) -> f64) -> Crossing {
    let Some(first) = rows.first() else {
        return Crossing::Beyond { d_max: 0.0 };
    };
    if of(first) >= target {
        return Crossing::AlreadyAbove;
    }
    for w in rows.windows(2) {
        let (lo, hi) = (&w[0], &w[1]);
        if of(hi) >= target {
            let (d_lo, d_hi) = (lo.deviation(), hi.deviation());
            let span = of(hi) - of(lo);
            let d = if span > 0.0 {
                d_lo + (target - of(lo)) / span * (d_hi - d_lo)
            } else {
                d_hi
            };
            return Crossing::At {
                d: d.clamp(d_lo, d_hi),
                k_lo: lo.k,
                k_hi: hi.k,
            };
        }
    }
    Crossing::Beyond {
        d_max: 1.0 - rows[rows.len() - 1].k as f64,
    }
}

/// 穿越点写成一行可读文本。
fn cell(c: Crossing) -> String {
    match c {
        Crossing::AlreadyAbove => "退化（基线已在其上）".to_string(),
        Crossing::At { d, .. } => format!("d={d:.3} → k={:.3}", 1.0 - d),
        Crossing::Beyond { d_max } => format!("扫描内未达到（上限 d={d_max:.2}）"),
    }
}

/// 近基线灵敏度：`d ∈ [0, 0.10]`（k 从 1.00 到 0.90）之间 ΔE00 中位数的增量 ÷ 0.10。
fn secant(r: &StemReport) -> Option<f64> {
    let lo = r.rows.iter().find(|x| (x.k - 0.9).abs() < 1e-6)?;
    let d = lo.deviation();
    if d <= 0.0 {
        return None;
    }
    Some((lo.median - r.rows[0].median) / d)
}

// ---------------------------------------------------------------------------
// 打印
// ---------------------------------------------------------------------------

fn print_stem(r: &StemReport) {
    let c = &r.chroma;
    eprintln!("\n── {} ──", r.stem);
    eprintln!("角色：{}", r.role);
    eprintln!(
        "尺寸 {}×{}   拟合样本 {}   评估样本 {}（留出集再对半劈：偶拟合 / 奇评估）",
        r.size.0, r.size.1, r.fit_n, r.eval_n
    );
    eprintln!(
        "注入强度：k=0.80 时平均 |Δ通道| = {:.5}（亮度守恒，仅色度被缩）",
        r.strength_080
    );

    eprintln!("\n【参考 C* 分布】（留出集 {} 个样本）", c.n);
    eprintln!("  中位 {:.2}   P90 {:.2}   P99 {:.2}", c.median, c.p90, c.p99);
    eprintln!(
        "  C* > 10（明显有色）{:>6.2}%      C* > 25（高彩度）{:>6.2}%",
        c.gt10 * 100.0,
        c.gt25 * 100.0
    );
    eprintln!(
        "  我们的原始解码 C* ÷ 参考 C* 中位：{:.3}（C*>10 子集 {} 个）",
        c.ours_ratio, c.n_chromatic
    );
    eprintln!(
        "  同一子集的 L* 中位：参考 {:.1}  vs  我们原始解码 {:.1}",
        c.l_theirs, c.l_ours
    );
    eprintln!(
        "  扣掉亮度差之后的**残余**色度倍率：{:.3}（{} 个样本；≈1 表示色度没缺、缺的只是亮度）",
        c.detoned_ratio, c.n_detoned
    );
    eprintln!("  ⚠ C* 倍率与 L* 差必须一起看：Lab 的 a*/b* 随亮度立方根缩放，亮度差会冒充色度差。");

    eprintln!("  残余色度倍率随参考 C* 的走势（>1 = 我们更艳，<1 = 我们更淡，100 表示 35 以上）：");
    for b in &c.bins {
        if b.n == 0 {
            eprintln!("    参考C* {:>4.0}–{:<4.0}  样本 0        ——（画面里没有这段内容）", b.lo, b.hi);
        } else {
            eprintln!(
                "    参考C* {:>4.0}–{:<4.0}  样本 {:>6}   倍率 {:.3}",
                b.lo, b.hi, b.n, b.detoned
            );
        }
    }
    let solid: Vec<&ChromaBin> = c.bins.iter().filter(|b| b.n >= 100).collect();
    if solid.len() >= 2 {
        let (first, last) = (solid[0], solid[solid.len() - 1]);
        eprintln!(
            "    走势：低彩度段 {:.3} → 高彩度段 {:.3}（{}）",
            first.detoned,
            last.detoned,
            if last.detoned < first.detoned * 0.9 {
                "随彩度升高明显变淡 —— 不是乘性常数偏差，而是高彩度段被压"
            } else {
                "基本持平 —— 看不出随彩度恶化"
            }
        );
    }

    eprintln!("\n【注入扫描】k = 色度保留倍率，d = 1 − k；主口径每个 k 都**重新拟合**曲线 + LUT");
    eprintln!(
        "       k      d   重拟合中位   不重拟合中位     P95        最大    注入后C*/参考C*"
    );
    for r in &r.rows {
        let tag = if (r.k - 1.0).abs() < 1e-6 { "  ←原状" } else { "" };
        eprintln!(
            "   {:>5.2}  {:>5.2}   {:>10.3}  {:>12.3}  {:>8.3}  {:>9.3}   {:>10.3}{tag}",
            r.k,
            r.deviation(),
            r.median,
            r.median_frozen,
            r.p95,
            r.max,
            r.chroma_ratio
        );
    }
    eprintln!(
        "   注：k=0.00 处中位数可能回落——输入彻底变灰后 LUT 沿中性轴给出一个平均色即可，\n\
         而 k=0.10 时残留的微小色度被放大 10 倍，反而落进互相冲突的格子。属极端外推，\n\
         不影响 d ≤ 0.6 区间的穿越点。"
    );

    let base = r.rows[0].median;
    if let Some(pub_median) = PUBLISHED_MEDIAN
        .iter()
        .find(|(s, _)| *s == r.stem)
        .map(|(_, m)| *m)
    {
        eprintln!(
            "\n基线对照：本诊断 k=1.0 得 {base:.3}，正式全密度拟合口径为 {pub_median:.3}（差 {:+.3}）",
            base - pub_median
        );
    }
}

fn print_conclusion(a: &StemReport, b: &StemReport) {
    let target_a = a.rows[0].median + TARGET_MEDIAN;
    let target_b = b.rows[0].median + TARGET_MEDIAN;
    let main_a = crossing(&a.rows, target_a, |r| r.median);
    let main_b = crossing(&b.rows, target_b, |r| r.median);

    eprintln!("\n════════ 穿越点：让 ΔE00 中位数达到阈值所需的饱和度偏差 d = 1 − k ════════");
    eprintln!("  阈值口径                {:<26}{}", a.stem, b.stem);
    for t in [TARGET_MEDIAN, 2.0, 3.0, 5.0] {
        eprintln!(
            "  绝对 中位 ≥ {t:<4}         {:<26}{}",
            cell(crossing(&a.rows, t, |r| r.median)),
            cell(crossing(&b.rows, t, |r| r.median))
        );
    }
    eprintln!(
        "  相对 基线 + {TARGET_MEDIAN:.1}（主）   {:<26}{}",
        cell(main_a),
        cell(main_b)
    );
    for (r, c) in [(a, main_a), (b, main_b)] {
        if let Crossing::At { k_lo, k_hi, .. } = c {
            eprintln!("    {} 的穿越点夹在 k={k_lo:.2} 与 k={k_hi:.2} 两行之间，按上表线性插值", r.stem);
        }
    }

    eprintln!("\n════════ 盲区倍数 ════════");
    match (main_a, main_b) {
        (Crossing::At { d: da, .. }, Crossing::At { d: db, .. }) => {
            eprintln!("  {}：偏差 d={db:.3}（k={:.3}）就把 ΔE00 中位数推高 {TARGET_MEDIAN:.1}", b.stem, 1.0 - db);
            eprintln!("  {}：同样要 d={da:.3}（k={:.3}）", a.stem, 1.0 - da);
            let factor = if db > 1e-9 { da / db } else { f64::INFINITY };
            eprintln!("  → 盲区倍数 = **{factor:.1}×**");
            eprintln!(
                "    即：同一处饱和度缺陷，在 {} 这类近中性样本上表现出的 ΔE00 恶化只有 {} 的 1/{factor:.1}。",
                a.stem, b.stem
            );
            if factor > 1.0 {
                eprintln!(
                    "    判读：**假设成立**——原来那 15 张竖构图样本对色度缺陷确实存在 {factor:.1} 倍盲区，"
                );
                eprintln!("          用它们验收'色度达标'不可信。");
            } else {
                eprintln!(
                    "    ⚠ 判读：**与预期相反**——近中性样本并不比彩色样本钝（{factor:.1}×），注入/拟合口径需要复核。"
                );
            }
        }
        _ => {
            eprintln!("  两张图的穿越点没能同时落在扫描区间内，无法给出倍数：");
            eprintln!("    {} → {}", a.stem, cell(main_a));
            eprintln!("    {} → {}", b.stem, cell(main_b));
        }
    }

    if let (Some(sa), Some(sb)) = (secant(a), secant(b)) {
        eprintln!("\n  补充口径（近基线灵敏度，d ∈ [0, 0.10]）：");
        eprintln!("    {} ΔE00中位 +{sa:.3} / 单位偏差", a.stem);
        eprintln!("    {} ΔE00中位 +{sb:.3} / 单位偏差", b.stem);
        if sa.abs() > 1e-9 {
            eprintln!("    → 灵敏度之比 = {:.1}×", sb / sa);
        }
    }

    // ---- 对照：重拟合到底补偿掉了多少（任务书警告的那个坑，这里把它量化）----
    eprintln!("\n════════ 对照：重拟合补偿掉多少 ════════");
    eprintln!("  （不重拟合 = 把 k=1.0 那份「出厂标定」直接套在注入后的数据上）");
    for d in [0.20_f64, 0.40, 0.60] {
        for r in [a, b] {
            let Some(row) = r.rows.iter().find(|x| (x.deviation() - d).abs() < 1e-6) else {
                continue;
            };
            let absorbed = if row.median_frozen > 1e-9 {
                1.0 - row.median / row.median_frozen
            } else {
                0.0
            };
            eprintln!(
                "   d={d:.2}（k={:.2}）{:>9}：重拟合 {:>6.3}  不重拟合 {:>6.3}  → 拟合吃掉 {:>5.1}%",
                row.k,
                r.stem,
                row.median,
                row.median_frozen,
                absorbed * 100.0
            );
        }
    }
    let frozen_a = crossing(&a.rows, a.rows[0].median_frozen + TARGET_MEDIAN, |r| {
        r.median_frozen
    });
    let frozen_b = crossing(&b.rows, b.rows[0].median_frozen + TARGET_MEDIAN, |r| {
        r.median_frozen
    });
    if let (Crossing::At { d: da, .. }, Crossing::At { d: db, .. }) = (frozen_a, frozen_b) {
        let f = if db > 1e-9 { da / db } else { f64::INFINITY };
        eprintln!("\n  若改按「不重拟合」口径，同一判据下：");
        eprintln!(
            "    {} 要 d={da:.3}，{} 要 d={db:.3} → 倍数为 {f:.1}×（{}）",
            a.stem,
            b.stem,
            if f > 1.0 { "彩色样本更敏感" } else { "仍旧不更敏感" }
        );
        eprintln!("    但验收测试会**重新标定**，所以主口径仍取重拟合那一列。");
    }

    eprintln!("\n  注：绝对阈值 {TARGET_MEDIAN:.1} 因两张图的基线都已在它之上而退化，故不作为主口径；");
    eprintln!("      主口径为「相对各自基线再高 {TARGET_MEDIAN:.1}」，两边定义一致、可比。");
}

// ---------------------------------------------------------------------------
// 测试
// ---------------------------------------------------------------------------

#[test]
fn chroma_blindspot_of_neutral_samples() {
    eprintln!("════════ 诊断 F：色度盲区量化 ════════");
    eprintln!("判据：ΔE00 中位数（ProPhoto 线性 → CIELAB D50）；达标线 {TARGET_MEDIAN:.1}");
    eprintln!("口径：曲线 {CURVE_BINS} 箱 + 3D LUT {LUT_EDGE}³（同脚手架）；**每个 k 都重新拟合**");
    eprintln!("注入：v' = lum + k·(v − lum)，lum 取 Rec.709 权重（亮度守恒，只缩放色度）");
    eprintln!(
        "扫描：k ∈ {:?}（任务书要求的 1.00/0.95/0.90/0.85/0.80/0.70/0.60 全含）",
        K_SWEEP
    );

    let reports: Vec<StemReport> = STEMS
        .iter()
        .filter_map(|(stem, role)| {
            let r = run_stem(stem, role);
            if r.is_none() {
                eprintln!("跳过 {stem}：样本缺席，无法做注入实验");
            }
            r
        })
        .collect();

    if reports.is_empty() {
        eprintln!("跳过本诊断：两张图的样本都不在位。");
        return;
    }

    for r in &reports {
        print_stem(r);
    }

    if reports.len() == 2 {
        print_conclusion(&reports[0], &reports[1]);
    } else {
        eprintln!("\n⚠ 只到位一张图，无法给出盲区倍数（需要 DSC_0001 与 DSC_0141 同时在场）。");
    }
}
