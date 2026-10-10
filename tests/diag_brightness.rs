//! 诊断 B：ΔE00 残差随**亮度**的分布。
//!
//! 5.5 的实测给出「未达标」的总体数字（中位数 2.177、P95 5.971），但总体分位数
//! 看不出残差是**处处都差**还是**只在某一段差**——这两者的修法完全不同：前者要换
//! 模型，后者只需给那一段更多自由度。本测试把总体数字按亮度拆成 16 段。
//!
//! **分箱键走 sRGB 编码域**：线性域里整个阴影挤在前两箱（线性 0.023 对应 L*≈18），
//! 暗部的问题会被高光样本淹掉；编码域等距的箱对应的感知跨度大致相当。这也正是
//! `fit::fit_curve` 的定义域（`fit::encode_srgb`），两者对齐后「哪一段欠拟合」
//! 可以直接与曲线的分箱分辨率对照。
//!
//! 拟合集与评估集由脚手架分开（`fitted_neutral_fixture`），因此下面是**留出误差**，
//! 不是训练误差。样本缺席时跳过而非失败，与仓库内其余依赖真实样张的测试一致。
//!
//! 最后一段是 **17³ LUT 网格占用**：曲线已经被搬到 sRGB 编码域拟合，但 `fit::fit_lut`
//! 仍用**线性**值索引网格（`src/fit.rs` 的 `idx`），`lut_edge = 17` → 线性间距 1/16。
//! 本段量化「17³ 的分辨率有多少真的落在样本上」，用来验证暗部是否被压进同一个格子。

mod common;

use common::*;
use nikonrawview::fit::{encode_srgb, forward_curve};
use std::cmp::{Ordering, Reverse};

/// 分箱数：编码域 16 等分，每箱宽 0.0625。
const BINS: usize = 16;
/// 达标线（见 `tasks.md` 5.5）。
const TARGET_MEDIAN: f64 = 1.0;
const TARGET_P95: f64 = 3.0;

// ---------------------------------------------------------------------------
// 工具
// ---------------------------------------------------------------------------

/// Rec.709 亮度。
///
/// `DiagSample::lum_ours()` 作用于**变换前**的值，所以「变换后亮度」要在这里另算：
/// 带符号偏差的两个端点必须用同一把尺子量。
fn lum(v: [f32; 3]) -> f64 {
    (0.2126 * v[0] + 0.7152 * v[1] + 0.0722 * v[2]) as f64
}

/// 中位数。
///
/// 刻意复用 `common::quantiles` 的**最近秩**实现，而不是另写一套：仓库里已经因为
/// `round` 与 `ceil` 的差别报错过一次中位数（见 5.5），两套百分位定义会互相打架。
fn median(v: &[f64]) -> f64 {
    quantiles(&mut v.to_vec()).0
}

/// 一组值的中位数；空集返回 `None`（空箱没有中位数可言，不能拿 NaN 冒充）。
fn med_of(v: &[f64]) -> Option<f64> {
    if v.is_empty() {
        None
    } else {
        Some(median(v))
    }
}

/// 一组值的 P95；空集返回 `None`。
fn p95_of(v: &[f64]) -> Option<f64> {
    if v.is_empty() {
        None
    } else {
        Some(quantiles(&mut v.to_vec()).1)
    }
}

/// 线性亮度落在哪一个编码域箱。
fn bin_of(linear: f64) -> usize {
    ((encode_srgb(linear as f32) * BINS as f32) as usize).min(BINS - 1)
}

/// 定宽数值列；空箱打 `—`。
fn col(opt: Option<f64>, prec: usize) -> String {
    let t = match opt {
        Some(v) => format!("{:.*}", prec, v),
        None => "—".to_string(),
    };
    format!("{t:>10}")
}

/// 复刻 `fit::fit_lut` 的网格索引（**就近取整**）。
///
/// 与 `src/fit.rs` 里 `fit_lut` 的 `idx` 逐字一致：拟合时每个样本投给最近的网格节点，
/// 节点的 `cnt` 决定它有没有数据。复刻一份是为了回答「17³ 的分辨率有多少真的被用到」
/// ——空节点只能靠邻域迭代从别处填，而那正是暗部被整体抬亮的一个来源。
fn node_index(v: [f32; 3], edge: usize) -> usize {
    let q = |x: f32| ((x.clamp(0.0, 1.0) * (edge - 1) as f32).round() as usize).min(edge - 1);
    (q(v[0]) * edge + q(v[1])) * edge + q(v[2])
}

/// 复刻 `fit::apply_lut` 的**插值格**索引（**向下取整**）。
///
/// 口径与上面不同：`apply_lut` 把样本所在的那一格用 8 个角做三线性插值，所以
/// 「误差被线性化在多小的范围内」要按这个口径数。
fn cell_index(v: [f32; 3], edge: usize) -> usize {
    let q = |x: f32| ((x.clamp(0.0, 1.0) * (edge - 1) as f32).floor() as usize).min(edge - 2);
    (q(v[0]) * edge + q(v[1])) * edge + q(v[2])
}

// ---------------------------------------------------------------------------
// 分箱
// ---------------------------------------------------------------------------

/// 一箱的汇总。
struct Bin {
    count: usize,
    /// 该箱全部样本的 ΔE00。
    de: Vec<f64>,
    /// 该箱全部样本的「变换后亮度 − 参考亮度」。
    bias: Vec<f64>,
    /// 该箱全部样本的参考亮度。
    ref_lum: Vec<f64>,
    /// 该箱内 ΔE00 > P95 目标线的样本数。
    over: usize,
}

impl Bin {
    fn new() -> Bin {
        Bin {
            count: 0,
            de: Vec::new(),
            bias: Vec::new(),
            ref_lum: Vec::new(),
            over: 0,
        }
    }
}

/// 误差最大的样本，便于回看具体像素。
struct Extreme {
    de: f64,
    xy: (usize, usize),
    /// 变换**前**的亮度（喂给基准变换的输入）。
    ours: f64,
    theirs: f64,
}

/// 一次分箱的全部汇总。
struct Tally {
    bins: Vec<Bin>,
    /// 全部样本的 ΔE00；顺序与 `keys` 一致。
    de: Vec<f64>,
    /// 全部样本的带符号亮度偏差。
    bias: Vec<f64>,
    /// 全部样本的分箱键值。
    keys: Vec<f64>,
    /// 尾部样本（ΔE00 > P95 目标线）的分箱键值。
    tail_keys: Vec<f64>,
    worst: Vec<Extreme>,
}

/// 按 `key` 给出的亮度分箱，并顺带收集整体统计量。
fn tally(samples: &[DiagSample], key: impl Fn(&DiagSample) -> f64) -> Tally {
    let mut t = Tally {
        bins: (0..BINS).map(|_| Bin::new()).collect(),
        de: Vec::with_capacity(samples.len()),
        bias: Vec::with_capacity(samples.len()),
        keys: Vec::with_capacity(samples.len()),
        tail_keys: Vec::new(),
        worst: Vec::with_capacity(samples.len()),
    };

    for s in samples {
        let de = s.delta_e();
        let ours = s.lum_ours() as f64;
        let theirs = s.lum_theirs() as f64;
        // 偏差 = 我们 − 参考：负值表示我们偏暗。
        let bias = lum(s.ours_transformed) - theirs;
        let k = key(s);

        let bin = &mut t.bins[bin_of(k)];
        bin.count += 1;
        bin.de.push(de);
        bin.bias.push(bias);
        bin.ref_lum.push(theirs);
        if de > TARGET_P95 {
            bin.over += 1;
            t.tail_keys.push(k);
        }

        t.de.push(de);
        t.bias.push(bias);
        t.keys.push(k);
        t.worst.push(Extreme {
            de,
            xy: s.xy,
            ours,
            theirs,
        });
    }
    t
}

// ---------------------------------------------------------------------------
// 打印
// ---------------------------------------------------------------------------

/// 打印一张 16 箱分布表。
///
/// 表头用 ASCII 是为了列对齐：中文是双宽字符，混在数字列里必然错位，因此列义写在
/// 表的上下方，表格本体保持纯数字。
fn print_table(t: &Tally, title: &str) {
    let n = t.de.len();
    let over_total: usize = t.bins.iter().map(|b| b.over).sum();
    println!();
    println!("{title}");
    println!("bin  enc_lo  enc_hi    count   share%     dE_med    dE_p95   bias_med    ref_med   tail%");
    for (i, b) in t.bins.iter().enumerate() {
        let lo = i as f32 / BINS as f32;
        let hi = (i + 1) as f32 / BINS as f32;
        let share = b.count as f64 / n as f64 * 100.0;
        let tail = if over_total == 0 {
            0.0
        } else {
            b.over as f64 / over_total as f64 * 100.0
        };
        println!(
            "{i:>3}  {lo:>6.4}  {hi:>6.4}  {:>7}  {:>7.2}%  {}  {}  {}  {}  {tail:>6.2}%",
            b.count,
            share,
            col(med_of(&b.de), 4),
            col(p95_of(&b.de), 4),
            col(med_of(&b.bias), 5),
            col(med_of(&b.ref_lum), 5),
        );
    }
    println!(
        "     合计 {n} 个样本 · ΔE00 > {TARGET_P95:.1} 的 {over_total} 个（{:.2}%）",
        over_total as f64 / n as f64 * 100.0
    );
    println!("     列义：enc_lo/enc_hi = 分箱键的编码域边界（左闭右开）");
    println!("           count/share% = 样本数 / 占本次评估集");
    println!("           dE_med/dE_p95 = 该箱 ΔE00 的中位数与 95 分位");
    println!("           bias_med = 该箱「变换后亮度 − 参考亮度」的中位数（Rec.709 线性域，负 = 我们偏暗）");
    println!("           ref_med  = 该箱参考亮度中位数；tail% = 该箱占全部超 P95 目标线样本的比例");
}

/// 打印三个关键问题的答案。
fn print_verdict(t: &mut Tally, label: &str) {
    let n = t.de.len();
    let mut all_de = t.de.clone();
    let (med, p95, max) = quantiles(&mut all_de);

    println!();
    println!("── 结论[{label}] ──");
    println!(
        "整体 ΔE00：中位数 {med:.4} · P95 {p95:.4} · 最大 {max:.4}（目标 ≤{TARGET_MEDIAN:.1} / ≤{TARGET_P95:.1}）→ {}",
        if med <= TARGET_MEDIAN && p95 <= TARGET_P95 {
            "达标"
        } else {
            "未达标"
        }
    );
    let n_over_1 = t.de.iter().filter(|&&d| d > TARGET_MEDIAN).count();
    let n_over_3 = t.de.iter().filter(|&&d| d > TARGET_P95).count();
    println!(
        "超阈样本：ΔE00 > {TARGET_MEDIAN:.1} 的 {n_over_1} 个（{:.1}%）；> {TARGET_P95:.1} 的 {n_over_3} 个（{:.1}%）",
        n_over_1 as f64 / n as f64 * 100.0,
        n_over_3 as f64 / n as f64 * 100.0
    );

    // 问题 1/2：残差是全局的还是集中在某段——用「箱中位数是否超目标线」界定「哪一段」。
    let over_bins: Vec<usize> = (0..BINS)
        .filter(|&i| t.bins[i].count > 0 && median(&t.bins[i].de) > TARGET_MEDIAN)
        .collect();
    let over_share = over_bins
        .iter()
        .map(|&i| t.bins[i].count)
        .sum::<usize>() as f64
        / n as f64
        * 100.0;
    println!();
    println!("1) 残差是全局的，还是集中在某个亮度段？");
    println!(
        "   16 箱中 {} 箱的 ΔE00 中位数 > {TARGET_MEDIAN:.1}；这些箱合计覆盖样本 {over_share:.2}%。",
        over_bins.len()
    );
    if over_bins.is_empty() {
        println!("   即：没有任何一段的中位数超标，残差是全局均匀的。");
    } else {
        for &i in &over_bins {
            let lo = i as f32 / BINS as f32;
            let hi = (i + 1) as f32 / BINS as f32;
            let b = &t.bins[i];
            let share = b.count as f64 / n as f64 * 100.0;
            println!(
                "   箱 {i:>2} 编码域 [{lo:.4}, {hi:.4})：样本 {:>7}（{share:.2}%）· ΔE00 中位 {:.4} · P95 {:.4}",
                b.count,
                median(&b.de),
                quantiles(&mut b.de.clone()).1
            );
        }
    }

    // 问题 2 的第二个角度：尾部样本（> P95 目标线）压在哪一段。
    let over_total: usize = t.bins.iter().map(|b| b.over).sum();
    if over_total > 0 {
        let mut by_tail: Vec<usize> = (0..BINS).filter(|&i| t.bins[i].over > 0).collect();
        by_tail.sort_by(|&a, &b| t.bins[b].over.cmp(&t.bins[a].over));
        println!();
        println!("   尾部样本（ΔE00 > {TARGET_P95:.1}，共 {over_total} 个）最集中的箱：");
        for &i in by_tail.iter().take(3) {
            let b = &t.bins[i];
            let share = b.count as f64 / n as f64;
            let tail_share = b.over as f64 / over_total as f64;
            let lo = i as f32 / BINS as f32;
            let hi = (i + 1) as f32 / BINS as f32;
            println!(
                "   箱 {i:>2} 编码域 [{lo:.4}, {hi:.4})：占样本 {:.2}% 却占尾部 {:.2}%（集中倍数 {:.2}×）· 该箱内超阈率 {:.1}%",
                share * 100.0,
                tail_share * 100.0,
                tail_share / share,
                b.over as f64 / b.count as f64 * 100.0
            );
        }
        let key_med_bin = bin_of(median(&t.keys));
        let tail_med_bin = bin_of(median(&t.tail_keys));
        println!(
            "   全体样本的分箱键中位数落在箱 {key_med_bin}，尾部样本落在箱 {tail_med_bin}——两者相距 {tail_med_bin} 与 {key_med_bin} 的箱差即为集中程度。"
        );
    }

    // 问题 3：带符号的亮度偏差。负 = 我们偏暗。
    let bmed = median(&t.bias);
    let dark = (0..BINS)
        .filter(|&i| t.bins[i].count > 0 && med_of(&t.bins[i].bias).is_some_and(|v| v < 0.0))
        .count();
    let bright = (0..BINS)
        .filter(|&i| t.bins[i].count > 0 && med_of(&t.bins[i].bias).is_some_and(|v| v > 0.0))
        .count();
    let dir = if bmed < 0.0 {
        "我们整体偏暗"
    } else if bmed > 0.0 {
        "我们整体偏亮"
    } else {
        "整体无偏"
    };
    println!();
    println!("3) 带符号的亮度偏差（我们 − 参考）：");
    println!("   整体中位数 {bmed:+.6}（线性域）→ {dir}；16 箱中 {dark} 箱中位偏暗、{bright} 箱中位偏亮。");

    // 最差样本：与 5.5 的「误差最大样本」对齐，便于回看像素。
    t.worst
        .sort_by(|a, b| b.de.partial_cmp(&a.de).unwrap_or(Ordering::Equal));
    println!();
    println!("误差最大的 5 个样本（坐标 = 裁剪后 x,y；比值 = 参考亮度 / 变换前亮度）：");
    for e in t.worst.iter().take(5) {
        let (x, y) = e.xy;
        let ratio = if e.ours > 0.0 {
            format!("{:.2}×", e.theirs / e.ours)
        } else {
            "—".to_string()
        };
        println!(
            "   ΔE00 {:>7.3}  位置 ({x:>4},{y:>4})  变换前亮度 {:.5}  参考亮度 {:.5}  {ratio}",
            e.de, e.ours, e.theirs
        );
    }
}

/// 17³ LUT 网格占用：名义分辨率有多少真的落在样本上，以及"坏"是不是坏在被压进同一格。
///
/// 动机：曲线已经搬到 sRGB 编码域拟合，而 `fit::fit_lut` 仍用**线性**值索引网格
/// （`src/fit.rs` 的 `idx`），`lut_edge = 17` → 线性间距 1/16 = 0.0625。线性域里
/// 整个阴影都落在第一个格子内，这正是 5.5 已经为曲线诊断过一次、但没在 3D 上修的病。
fn print_grid_occupancy(fx: &FittedFixture) {
    let edge = fx.transform.lut_edge;
    let curve = &fx.transform.curve;
    let n = fx.eval.len();
    let n_nodes = edge * edge * edge;

    let mut node_cnt = vec![0usize; n_nodes];
    let mut node_de: Vec<Vec<f64>> = (0..n_nodes).map(|_| Vec::new()).collect();
    let mut node_over = vec![0usize; n_nodes];
    let mut cell_cnt = vec![0usize; n_nodes];
    let mut cell_over = vec![0usize; n_nodes];
    let mut cell0_de: Vec<f64> = Vec::new();
    let mut all_de: Vec<f64> = Vec::with_capacity(n);
    let mut total_over = 0usize;

    for s in &fx.eval {
        // 与 fit_lut 一致：网格的输入是**曲线之后**的线性值，不是原始值。
        let c = forward_curve(curve, s.ours_linear);
        let de = s.delta_e();
        let over = de > TARGET_P95;

        let i = node_index(c, edge);
        node_cnt[i] += 1;
        node_de[i].push(de);
        if over {
            node_over[i] += 1;
        }

        let ci = cell_index(c, edge);
        cell_cnt[ci] += 1;
        if over {
            cell_over[ci] += 1;
        }
        if ci == 0 {
            cell0_de.push(de);
        }
        all_de.push(de);
        if over {
            total_over += 1;
        }
    }

    let empty = node_cnt.iter().filter(|&&c| c == 0).count();
    let tiny = node_cnt.iter().filter(|&&c| (1..10).contains(&c)).count();
    let small = node_cnt.iter().filter(|&&c| (10..100).contains(&c)).count();
    let mid = node_cnt.iter().filter(|&&c| (100..1000).contains(&c)).count();
    let big = node_cnt.iter().filter(|&&c| c >= 1000).count();
    assert_eq!(
        empty + tiny + small + mid + big,
        n_nodes,
        "占用分档必须覆盖全部节点"
    );

    let mut used: Vec<usize> = (0..n_nodes).filter(|&i| node_cnt[i] > 0).collect();
    used.sort_by_key(|&i| Reverse(node_cnt[i]));

    let global_med = median(&all_de);
    let spend = 1.0 / (edge - 1) as f64;
    let tail_pct = |k: usize| {
        if total_over == 0 {
            0.0
        } else {
            k as f64 / total_over as f64 * 100.0
        }
    };

    println!();
    println!(
        "── 17³ LUT 网格占用（lut_edge={edge}，**线性**等距，间距 1/{} = {spend:.4}）──",
        edge - 1
    );
    println!(
        "口径：索引 = (q(r)·edge + q(g))·edge + q(b)，q(x) = round(x.clamp(0,1)·(edge-1))——与 src/fit.rs 的 fit_lut 逐字一致"
    );
    println!("样本：评估集 {n} 个（与拟合集同样本量、同分布，占用图景可直接迁移到 fit_lut 的 cnt）");
    println!(
        "总格数 {n_nodes}；有样本的格 {} 个（{:.2}%）；空格 {empty} 个（{:.2}%）",
        n_nodes - empty,
        (n_nodes - empty) as f64 / n_nodes as f64 * 100.0,
        empty as f64 / n_nodes as f64 * 100.0
    );
    println!(
        "格内样本数分档：0 个 {empty} · 1~9 个 {tiny} · 10~99 个 {small} · 100~999 个 {mid} · ≥1000 个 {big}"
    );
    println!();
    println!("占用最多的前 5 格（索引 0 = 格坐标 (0,0,0)，即全暗角）：");
    println!("rank   index    r    g    b    count   share%     dE_med   tail%");
    for (k, &i) in used.iter().take(5).enumerate() {
        let r = i / (edge * edge);
        let g = (i / edge) % edge;
        let b = i % edge;
        let rank = k + 1;
        println!(
            "{rank:>4}  {i:>6}  {r:>3}  {g:>3}  {b:>3}  {:>7}  {:>6.2}%  {}  {:>6.2}%",
            node_cnt[i],
            node_cnt[i] as f64 / n as f64 * 100.0,
            col(med_of(&node_de[i]), 4),
            tail_pct(node_over[i]),
        );
    }

    let c0 = node_cnt[0];
    let m0 = match med_of(&node_de[0]) {
        Some(m) => format!("{m:.4}"),
        None => "—".to_string(),
    };
    let ratio0 = match med_of(&node_de[0]) {
        Some(m) if global_med > 0.0 => format!("{:.2}×", m / global_med),
        _ => "—".to_string(),
    };
    println!();
    println!(
        "首格（拟合口径，索引 0）：样本 {c0} 个（{:.2}%）· ΔE00 中位 {m0} vs 全局 {global_med:.4} → {ratio0} · 占全部 ΔE00 > {TARGET_P95:.1} 尾部的 {:.2}%",
        c0 as f64 / n as f64 * 100.0,
        tail_pct(node_over[0])
    );

    let cell0 = cell_cnt[0];
    let cm0 = match med_of(&cell0_de) {
        Some(m) => format!("{m:.4}"),
        None => "—".to_string(),
    };
    println!(
        "首格（插值口径，= apply_lut 的 floor 格 (0,0,0)，覆盖线性 [0,{spend:.4})³ = 编码域 [0,{:.4})³）：样本 {cell0} 个（{:.2}%）· ΔE00 中位 {cm0} · 占全部 ΔE00 > {TARGET_P95:.1} 尾部的 {:.2}%",
        encode_srgb(spend as f32),
        cell0 as f64 / n as f64 * 100.0,
        tail_pct(cell_over[0])
    );
    println!(
        "含义：{:.2}% 的格子完全没有样本，占用前 5 格就吃掉 {:.2}% 的样本——17³ 的名义分辨率几乎全花在没人去的亮部；",
        empty as f64 / n_nodes as f64 * 100.0,
        used.iter().take(5).map(|&i| node_cnt[i]).sum::<usize>() as f64 / n as f64 * 100.0
    );
    println!(
        "      而被压进同一格的样本只能用 8 个角做三线性插值，暗部的强非线性修正因此被抹平。"
    );
}

// ---------------------------------------------------------------------------
// 测试
// ---------------------------------------------------------------------------

/// ΔE00 残差随亮度的分布。
///
/// 本测试只做**诊断**：它断言分箱自洽（覆盖全部样本、无 NaN），但**不**断言达标与否
/// ——当前已知未达标，把结论写死成断言只会让用例变红而不会带来信息。达标判定留给
/// 5.6（拟合/验证分离后的验收）。
#[test]
fn delta_e_residual_by_luminance_bin() {
    let Some(fx) = fitted_neutral_fixture() else {
        // 脚手架已打印「跳过 …」提示：样本缺席时跳过而非失败。
        return;
    };
    assert!(!fx.eval.is_empty(), "评估集不应为空");

    // 主表：按下标（变换前）亮度分箱——这是曲线的定义域。
    let mut by_ours = tally(&fx.eval, |s| s.lum_ours() as f64);
    // 对照表：按参考亮度分箱。若两张表的集中段落在同一条亮度带上，说明集中段是
    // 真实的亮度段，而不是分箱键选错造成的假象。
    let mut by_ref = tally(&fx.eval, |s| s.lum_theirs() as f64);

    let n = fx.eval.len();
    for (name, t) in [("主表", &by_ours), ("对照表", &by_ref)] {
        assert_eq!(
            t.bins.iter().map(|b| b.count).sum::<usize>(),
            n,
            "{name}的分箱必须覆盖全部样本（漏样本会让所有占比失去意义）"
        );
        assert!(
            t.de.iter().all(|d| d.is_finite()),
            "{name}的 ΔE00 不应出现 NaN/inf"
        );
    }

    println!();
    println!("── 诊断 B：ΔE00 残差随亮度的分布 ──");
    println!(
        "样本：评估集 {n} 个 · 拟合集 {} 个 · 参考尺寸 {}×{}（两者用互不相交的像素子集）",
        fx.fit_count, fx.size.0, fx.size.1
    );
    println!(
        "分箱键：sRGB 编码域等距 {BINS} 箱（与 fit::fit_curve 的定义域一致）；重采样步长 13px，占全图 1/169"
    );

    print_table(&by_ours, "表 1（主表）：键 = sRGB 编码(变换前亮度)");
    print_verdict(&mut by_ours, "主表，键 = 变换前亮度");

    // 病根量化：暗部残差是不是"因为被压进同一个 LUT 格子"。
    print_grid_occupancy(&fx);

    print_table(&by_ref, "表 2（对照）：键 = sRGB 编码(参考亮度)");
    print_verdict(&mut by_ref, "对照表，键 = 参考亮度");
}
