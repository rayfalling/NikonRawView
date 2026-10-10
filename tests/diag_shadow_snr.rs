//! 诊断 C：暗部的大误差像素，是"落在噪声里"还是"系统性偏差"？
//!
//! # 要回答的问题
//!
//! 留出评估集的 ΔE00 未达标（中位数≈2.1、P95≈5.9），而误差最大的样本全在深阴影：
//! 我们读出接近中性，尼康读出明显偏暖。若这些像素的信号本身就在噪声底噪附近，那这个
//! "差异"就不可匹配——继续拟合只会把噪声学进去，应当排除而不是硬拟合。
//!
//! # 采样方式（一处必须写清的取舍）
//!
//! 脚手架 `common::fitted_neutral_fixture()` 的评估样本是**抽样过**的：以 13 像素为
//! 步长取点、再取其中的奇数序号，样本之间最近也隔着 26 像素——「xy 周围的 3×3 邻域」
//! 在 eval 里根本取不到。这是任务里点名的陷阱。
//!
//! 因此本测试**完整沿用脚手架的拟合流程**（色调曲线、3D LUT、留出评估集全部来自
//! `fitted_neutral_fixture()`），另外用**逐字相同的解码参数**把同一张 NEF 再解码一次，
//! 只为取每个评估样本真实的 3×3 邻域。解码已被验证逐位可复现（任务 2.7），所以两次
//! 解码的同一坐标必然是同一个像素——第 [1] 节对此下了断言：坐标映射一旦错了会立刻
//! 失败，而不是悄悄给出一份错的统计。
//!
//! # 两个噪声口径
//!
//! 规格要的是「邻域亮度均值 / 邻域亮度标准差」。但 3×3 的标准差里**混着真实纹理**，
//! 会把噪声估高。因此另给一个高频底噪口径：4 邻域 Laplacian（它抵消线性梯度，白噪声
//! 下方差为 20σ²）的稳健估计。两者都打印；判断"差异是否在噪声之上"以高频底噪为准。
//!
//! # 暗部 10% 的统计口径
//!
//! 「我们亮度最低的 10%」取的是**评估样本**按 `lum_ours` 排序后的最低 10%——评估样本
//! 是像素的均匀抽样，所以这与"最暗的 10% 像素"在分布上等价（要逐像素做符号统计就得把
//! 272 MB 的参考 TIF 再读一遍并转换，本诊断不需要那么多像素来定方向）。
//!
//! # 附带：拍摄参数
//!
//! 第 [9] 节顺带把 `simple/` 里各张 NEF 的 ISO / 曝光 / 光圈 / 拍摄时间读出来。它服务于
//! 一个交叉判读：若"暗部大误差来自噪声"，不同照片的噪声水平不同，ΔE00 的 P95 就不该在
//! 四张图上几乎一模一样。ISO 是否一致，正是这个解释能否成立的前提之一。
//!
//! # 附带：分阶段归因
//!
//! 第 [11] 节把同一评估集分别过「恒等」「仅曲线」「曲线+LUT」，看误差是哪一级引入的——
//! 排除噪声之后，这一步才回答"那到底是什么造成的"。

mod common;

use common::{fitted_neutral_fixture, samples_dir};
use nikonrawview::libraw::{decode_working_space, Decoded, Demosaic, Options};
use nikonrawview::tiff::Tiff;
use nikonrawview::transform::BaseTransform;
use std::cmp::Ordering;
use std::path::{Path, PathBuf};

/// 局部窗口边长（3×3）。
const WIN: usize = 3;
/// 高斯下 E|Z|。
const E_ABS_Z: f64 = 0.797_884_560_802_865;
/// 高斯下 median|Z| / E|Z|。
///
/// 本诊断在**成百上千个像素上取 |L| 的中位数**，而 median|L| = 0.6745·E|L|；若还用
/// E|L| 的常数换算 σ，底噪会被系统性低估约 1/3（比值则被高估同样比例）。
const MEDIAN_ABS_Z: f64 = 0.674_489_750_196_081_7;
/// 4 邻域 Laplacian 在白噪声下的方差系数：Var(L) = 20σ²。
const LAPLACIAN_VAR: f64 = 20.0;
/// "我们读出中性"的通道极差阈值（线性域）。
const NEUTRAL_SPREAD: f64 = 0.005;
/// "参考读出有色"的通道极差阈值（线性域）。
const TINTED_SPREAD: f64 = 0.02;
/// SNR 分箱。
const SNR_BINS: &[(f64, f64, &str)] = &[
    (0.0, 1.0, "<1"),
    (1.0, 2.0, "1-2"),
    (2.0, 4.0, "2-4"),
    (4.0, 8.0, "4-8"),
    (8.0, 16.0, "8-16"),
    (16.0, 32.0, "16-32"),
    (32.0, f64::INFINITY, ">=32"),
];

// ---------------------------------------------------------------------------
// 局部统计
// ---------------------------------------------------------------------------

/// 一个评估样本的 3×3 邻域统计。
#[derive(Debug, Clone, Copy)]
struct Local {
    /// 邻域亮度均值（信号）。
    lum_mean: f64,
    /// **中心像素**的亮度——用于和脚手架样本对齐做坐标自检。
    lum_center: f64,
    /// 邻域亮度标准差（噪声的**上界**：里面还含真实纹理）。
    lum_sd: f64,
    /// 高频底噪：Laplacian 稳健估计，受纹理影响远小于标准差。
    hf_sd: f64,
    /// 亮度信噪比 = 邻域亮度均值 / 邻域亮度标准差。
    snr: f64,
    /// 邻域色度信号模长 |均值(r−g, b−g)|。
    chroma_signal: f64,
    /// 邻域色度噪声：两个色度分量的标准差合并。
    chroma_sd: f64,
    /// 色度信噪比。
    snr_chroma: f64,
}

/// 一个评估样本 + 它的局部统计。
struct Row {
    /// 在留出评估集中的序号。
    idx: usize,
    xy: (usize, usize),
    de: f64,
    ours_linear: [f32; 3],
    ours_transformed: [f32; 3],
    theirs_linear: [f32; 3],
    lum_ours: f64,
    lum_theirs: f64,
    local: Local,
}

fn cmp_f64(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

/// 安全比值：分母为 0 时返回 `+∞`（而不是 NaN），便于进"最高"那一箱。
fn ratio_or_inf(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        a / b
    } else {
        f64::INFINITY
    }
}

/// 按 Rec.709 权重取 16 位像素的亮度（0..1）。
fn lum(p: [u16; 3]) -> f64 {
    (0.2126 * p[0] as f64 + 0.7152 * p[1] as f64 + 0.0722 * p[2] as f64) / 65535.0
}

/// 均值与**样本**标准差（n−1）。
fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0);
    (mean, var.max(0.0).sqrt())
}

fn median_of(mut v: Vec<f64>) -> f64 {
    common::quantiles(&mut v).0
}

/// 按 ΔE00 从大到小排序的下标。
fn order_by_de(rows: &[Row]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..rows.len()).collect();
    idx.sort_by(|&a, &b| cmp_f64(rows[b].de, rows[a].de));
    idx
}

/// 每组的中位"明度项占比"（ΔE00 的平方分解，其余为色度侧）。
fn lightness_share_median(v: &[&Row]) -> f64 {
    median_of(
        v.iter()
            .map(|r| {
                nikonrawview::deltae::delta_e_prophoto_decomposed(r.ours_transformed, r.theirs_linear)
                    .lightness_share()
            })
            .collect(),
    )
}

/// 噪声足迹：把 ours 沿"参考−我们"的方向移动 `k·σ_hf`，经**同一个变换**后 ΔE00 变多少。
///
/// 方向已经对准误差本身，所以这是"这么大的噪声最多能造成多大 ΔE00"的**上限**式估计。
fn noise_footprint(t: &BaseTransform, r: &Row, k: f64) -> f64 {
    let o = r.ours_linear;
    let e = [
        r.theirs_linear[0] as f64 - o[0] as f64,
        r.theirs_linear[1] as f64 - o[1] as f64,
        r.theirs_linear[2] as f64 - o[2] as f64,
    ];
    let len = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
    if len <= 0.0 {
        return 0.0;
    }
    let s = k * r.local.hf_sd / len;
    let moved = [
        (o[0] as f64 + e[0] * s) as f32,
        (o[1] as f64 + e[1] * s) as f32,
        (o[2] as f64 + e[2] * s) as f32,
    ];
    nikonrawview::deltae::delta_e_prophoto(t.apply(moved), r.ours_transformed)
}

/// 取 `(cx, cy)` 周围 3×3 的局部统计；邻域越出解码图时返回 `None`。
fn local_stats(d: &Decoded, cx: usize, cy: usize) -> Option<Local> {
    let r = WIN / 2;
    if cx < r || cy < r || cx + r >= d.width || cy + r >= d.height {
        return None;
    }
    let inv = 1.0 / 65535.0;
    let mut l = [0.0f64; WIN * WIN];
    let mut rg = [0.0f64; WIN * WIN];
    let mut bg = [0.0f64; WIN * WIN];
    let mut i = 0usize;
    for dy in 0..WIN {
        for dx in 0..WIN {
            let p = d.at(cx + dx - r, cy + dy - r)?;
            l[i] = lum(p);
            rg[i] = (p[0] as f64 - p[1] as f64) * inv;
            bg[i] = (p[2] as f64 - p[1] as f64) * inv;
            i += 1;
        }
    }
    let (lum_mean, lum_sd) = mean_sd(&l);
    let (m_rg, sd_rg) = mean_sd(&rg);
    let (m_bg, sd_bg) = mean_sd(&bg);
    // 行主序 3×3：中心是下标 4，上下左右是 1/7/3/5。
    // 白噪声下 E|L| = E|Z|·√20·σ；本诊断取的是 |L| 的中位数，故用 median|Z| 的常数。
    let lap = (l[1] + l[7] + l[3] + l[5] - 4.0 * l[4]).abs();
    let hf_sd = lap / (MEDIAN_ABS_Z * E_ABS_Z * LAPLACIAN_VAR.sqrt());
    let chroma_signal = (m_rg * m_rg + m_bg * m_bg).sqrt();
    let chroma_sd = ((sd_rg * sd_rg + sd_bg * sd_bg) / 2.0).sqrt();
    Some(Local {
        lum_mean,
        lum_center: l[4],
        lum_sd,
        hf_sd,
        snr: ratio_or_inf(lum_mean, lum_sd),
        chroma_signal,
        chroma_sd,
        snr_chroma: ratio_or_inf(chroma_signal, chroma_sd),
    })
}

/// 三通道极差（线性域）。
fn spread(v: [f32; 3]) -> f64 {
    let mx = v.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let mn = v.iter().copied().fold(f32::INFINITY, f32::min);
    (mx - mn) as f64
}

/// SNR 落在哪一箱。
fn snr_bin(snr: f64) -> Option<usize> {
    SNR_BINS
        .iter()
        .position(|&(lo, hi, _)| snr >= lo && (snr < hi || hi.is_infinite()))
}

/// Spearman 秩相关（平均秩处理并列）。
fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let rx = ranks(x);
    let ry = ranks(y);
    let n = rx.len() as f64;
    if rx.len() != ry.len() || rx.len() < 3 {
        return f64::NAN;
    }
    let mx = rx.iter().sum::<f64>() / n;
    let my = ry.iter().sum::<f64>() / n;
    let cov: f64 = rx.iter().zip(&ry).map(|(a, b)| (a - mx) * (b - my)).sum();
    let vx: f64 = rx.iter().map(|a| (a - mx) * (a - mx)).sum();
    let vy: f64 = ry.iter().map(|b| (b - my) * (b - my)).sum();
    if vx > 0.0 && vy > 0.0 {
        cov / (vx.sqrt() * vy.sqrt())
    } else {
        f64::NAN
    }
}

fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| cmp_f64(v[a], v[b]));
    let mut out = vec![0.0f64; v.len()];
    let mut i = 0usize;
    while i < idx.len() {
        let mut j = i + 1;
        while j < idx.len() && cmp_f64(v[idx[j]], v[idx[i]]) == Ordering::Equal {
            j += 1;
        }
        // 并列的一段占据秩 i+1..=j，取平均秩。
        let avg = (i + 1 + j) as f64 / 2.0;
        for k in &idx[i..j] {
            out[*k] = avg;
        }
        i = j;
    }
    out
}

// ---------------------------------------------------------------------------
// 报告
// ---------------------------------------------------------------------------

/// [2] 按 SNR 分箱，给出每箱的 ΔE00 中位数、P95 与样本占比。
fn report_bins(rows: &[Row]) {
    eprintln!("\n[2] 按局部信噪比（3×3 亮度 均值/标准差）分箱");
    eprintln!("    箱       样本数    占比      ΔE00中位  ΔE00_P95  ΔE00>5占比  我们亮度中位  SNR中位");
    let total = rows.len() as f64;
    let inf_n = rows.iter().filter(|r| r.local.snr.is_infinite()).count();
    for (b, &(_, _, label)) in SNR_BINS.iter().enumerate() {
        let sel: Vec<&Row> = rows.iter().filter(|r| snr_bin(r.local.snr) == Some(b)).collect();
        if sel.is_empty() {
            eprintln!("    {label:>6}         0      0.0%   (空箱)");
            continue;
        }
        let mut des: Vec<f64> = sel.iter().map(|r| r.de).collect();
        let (med, p95, _) = common::quantiles(&mut des);
        let share = 100.0 * sel.len() as f64 / total;
        let over5 = 100.0 * sel.iter().filter(|r| r.de > 5.0).count() as f64 / sel.len() as f64;
        let lum_med = median_of(sel.iter().map(|r| r.lum_ours).collect());
        let snr_med = median_of(sel.iter().map(|r| r.local.snr).collect());
        eprintln!(
            "    {label:>6}  {:>8}  {share:>6.1}%  {med:>8.3}  {p95:>8.3}  {over5:>9.1}%  {lum_med:>11.4}  {snr_med:>6.2}",
            sel.len()
        );
    }
    if inf_n > 0 {
        eprintln!("    （其中 {inf_n} 个样本的 3×3 完全平直，SNR 记为 +∞，落在最高箱）");
    }
}

/// [3] 相关性：整体、暗部内部、中间调内部。
fn report_correlation(rows: &[Row], dark: &[&Row], mid: &[&Row]) {
    let pairs = |v: &[&Row]| -> (Vec<f64>, Vec<f64>) {
        v.iter()
            .filter(|r| r.local.snr.is_finite())
            .map(|r| (r.local.snr, r.de))
            .unzip()
    };
    let all: Vec<&Row> = rows.iter().collect();
    let (snr_all, de_all) = pairs(&all);
    let (snr_dark, de_dark) = pairs(dark);
    let (snr_mid, de_mid) = pairs(mid);

    let lums: Vec<f64> = rows.iter().map(|r| r.lum_ours).collect();
    let des: Vec<f64> = rows.iter().map(|r| r.de).collect();

    eprintln!("\n[3] 相关性（Spearman 秩相关：−1 完全反向 / 0 无关 / +1 完全同向）");
    eprintln!(
        "    全体 {} 个样本：ρ(SNR, ΔE00) = {:.3}",
        snr_all.len(),
        spearman(&snr_all, &de_all)
    );
    eprintln!("    对照：ρ(我们亮度, ΔE00) = {:.3}   ← 亮度本身与误差的关系", spearman(&lums, &des));
    // 局部方差本身（不除以均值）与 ΔE00 的关系：SNR 里含了均值，这里单独看"噪声大小"。
    let sds: Vec<f64> = rows.iter().map(|r| r.local.lum_sd).collect();
    let hfs: Vec<f64> = rows.iter().map(|r| r.local.hf_sd).collect();
    eprintln!(
        "    对照：ρ(3×3 σ, ΔE00) = {:.3}，ρ(σ_hf, ΔE00) = {:.3}   ← 局部方差大是否伴随误差大",
        spearman(&sds, &des),
        spearman(&hfs, &des)
    );
    eprintln!(
        "    暗部 10% 内部（{} 个）：ρ(SNR, ΔE00) = {:.3}   ← **控制住亮度之后** SNR 还有没有解释力",
        snr_dark.len(),
        spearman(&snr_dark, &de_dark)
    );
    eprintln!(
        "    中间调内部（{} 个）：ρ(SNR, ΔE00) = {:.3}",
        snr_mid.len(),
        spearman(&snr_mid, &de_mid)
    );

    if dark.len() >= 20 {
        let cut = median_of(dark.iter().map(|r| r.local.snr).collect());
        let lo_half: Vec<f64> = dark
            .iter()
            .filter(|r| r.local.snr <= cut)
            .map(|r| r.de)
            .collect();
        let hi_half: Vec<f64> = dark
            .iter()
            .filter(|r| r.local.snr > cut)
            .map(|r| r.de)
            .collect();
        eprintln!(
            "    暗部内部按 SNR 中位（{cut:.2}）对半切：低 SNR 半 ΔE00 中位 {:.3}（n={}），高 SNR 半 {:.3}（n={}）",
            median_of(lo_half),
            dark.iter().filter(|r| r.local.snr <= cut).count(),
            median_of(hi_half),
            dark.iter().filter(|r| r.local.snr > cut).count(),
        );
    }
}

/// [4] 误差最大那一档的画像。
fn report_top_profile(rows: &[Row], dark: &[&Row]) {
    let mut by_de: Vec<usize> = (0..rows.len()).collect();
    by_de.sort_by(|&a, &b| cmp_f64(rows[b].de, rows[a].de));
    let top = (rows.len() / 10).max(1);
    let worst: Vec<&Row> = by_de[..top].iter().map(|&i| &rows[i]).collect();
    let rest: Vec<&Row> = by_de[top..].iter().map(|&i| &rows[i]).collect();

    let profile = |v: &[&Row]| -> (f64, f64, f64, f64) {
        (
            median_of(v.iter().map(|r| r.de).collect()),
            median_of(v.iter().map(|r| r.local.snr).collect()),
            100.0 * v.iter().filter(|r| r.local.snr < 2.0).count() as f64 / v.len() as f64,
            median_of(v.iter().map(|r| r.lum_ours).collect()),
        )
    };
    let (de_w, snr_w, low_w, lum_w) = profile(&worst);
    let (de_r, snr_r, low_r, lum_r) = profile(&rest);
    let low_all = 100.0 * rows.iter().filter(|r| r.local.snr < 2.0).count() as f64 / rows.len() as f64;

    eprintln!("\n[4] 误差最大的一档（ΔE00 最高的 10%，n={top}）");
    eprintln!("                     ΔE00中位   SNR中位   SNR<2占比   我们亮度中位");
    eprintln!("    ΔE00 最高 10%     {de_w:>7.3}   {snr_w:>7.2}   {low_w:>8.1}%   {lum_w:>11.4}");
    eprintln!("    其余 90%          {de_r:>7.3}   {snr_r:>7.2}   {low_r:>8.1}%   {lum_r:>11.4}");
    eprintln!("    全体              {:>7.3}   {:>7.2}   {low_all:>8.1}%   {:>11.4}",
        median_of(rows.iter().map(|r| r.de).collect()),
        median_of(rows.iter().map(|r| r.local.snr).collect()),
        median_of(rows.iter().map(|r| r.lum_ours).collect()),
    );

    let dark_ids: std::collections::BTreeSet<usize> = dark.iter().map(|r| r.idx).collect();
    let over5: Vec<&Row> = rows.iter().filter(|r| r.de > 5.0).collect();
    if !over5.is_empty() {
        let n_low = over5.iter().filter(|r| r.local.snr < 2.0).count();
        let n_dark = over5.iter().filter(|r| dark_ids.contains(&r.idx)).count();
        eprintln!(
            "    ΔE00 > 5 的 {} 个样本：SNR<2 占 {:.1}%，落在最暗 10% 占 {:.1}%（全体 SNR<2 只占 {low_all:.1}%）",
            over5.len(),
            100.0 * n_low as f64 / over5.len() as f64,
            100.0 * n_dark as f64 / over5.len() as f64,
        );
    }
}

/// [5] 差异到底在不在噪声之上。
///
/// 三条互补的证据：
/// 1. 差异量 / 局部底噪：|Δ亮度|/σ_hf、|Δ色度|/色度σ；
/// 2. **噪声足迹**：把 ours 沿"误差方向"整体移动 1σ（以及保守的 3σ），经过同一变换后
///    ΔE00 会变多少——这是"这么大的噪声最多能造成多大 ΔE00"的直接估计；
/// 3. 实际 ΔE00 与足迹的比值。
fn report_noise_vs_signal(rows: &[Row], dark: &[&Row], t: &BaseTransform) {
    let dnr_lum = |r: &Row| ratio_or_inf((r.lum_theirs - r.lum_ours).abs(), r.local.hf_sd);
    let dnr_lum3 = |r: &Row| ratio_or_inf((r.lum_theirs - r.lum_ours).abs(), r.local.lum_sd);
    // 色度差异：参考的 (r−g, b−g) 与我们的之差，除以我们的局部色度 σ。
    let dnr_chroma = |r: &Row| {
        let (o, v) = (r.ours_linear, r.theirs_linear);
        let dg = (v[0] as f64 - v[1] as f64) - (o[0] as f64 - o[1] as f64);
        let db = (v[2] as f64 - v[1] as f64) - (o[2] as f64 - o[1] as f64);
        ratio_or_inf((dg * dg + db * db).sqrt(), r.local.chroma_sd)
    };
    // 噪声足迹：沿"参考−我们"的方向移动 k·σ_hf，再走同一变换看 ΔE00 变多少。
    // 这是**上限**式估计——方向已经对准了误差本身，噪声若真能解释误差，足迹就该接近实际。
    let col = |f: &dyn Fn(&Row) -> f64| -> (f64, f64) {
        (
            median_of(dark.iter().map(|r| f(r)).collect()),
            median_of(rows.iter().map(f).collect()),
        )
    };
    let (d1, a1) = col(&dnr_lum);
    let (d2, a2) = col(&dnr_lum3);
    let (d3, a3) = col(&dnr_chroma);
    let (d4, a4) = col(&|r: &Row| noise_footprint(t, r, 1.0));
    let (d5, a5) = col(&|r: &Row| noise_footprint(t, r, 3.0));
    let (d6, a6) = col(&|r: &Row| r.local.snr_chroma);
    let (d7, a7) = col(&|r: &Row| r.local.lum_mean);
    let (d8, a8) = col(&|r: &Row| r.local.chroma_signal);
    let de_dark = median_of(dark.iter().map(|r| r.de).collect());
    let de_all = median_of(rows.iter().map(|r| r.de).collect());

    eprintln!("\n[5] 差异在不在噪声之上（中位数）");
    eprintln!("    量                                    暗部 10%     全体");
    eprintln!("    邻域亮度均值（信号）                     {d7:>8.4}   {a7:>8.4}");
    eprintln!("    |Δ亮度| / 高频底噪 σ_hf                 {d1:>8.2}   {a1:>8.2}");
    eprintln!("    |Δ亮度| / 3×3 σ（含纹理，偏保守）        {d2:>8.2}   {a2:>8.2}");
    eprintln!("    |Δ色度| / 我们局部色度 σ                {d3:>8.2}   {a3:>8.2}");
    eprintln!("    我们局部色度信号模长                     {d8:>8.4}   {a8:>8.4}");
    eprintln!("    我们色度信噪比（色度信号/色度σ）          {d6:>8.2}   {a6:>8.2}");
    eprintln!("    噪声足迹：沿误差方向移 1σ 的 ΔE00        {d4:>8.2}   {a4:>8.2}");
    eprintln!("    噪声足迹：沿误差方向移 3σ 的 ΔE00        {d5:>8.2}   {a5:>8.2}");
    eprintln!("    实际 ΔE00                               {de_dark:>8.2}   {de_all:>8.2}");
    eprintln!(
        "    ⇒ 实际误差 / 3σ 噪声足迹 = {:.1}（暗部）；1σ = {:.1}",
        ratio_or_inf(de_dark, d5),
        ratio_or_inf(de_dark, d4)
    );
}

/// 暗部的三通道方向统计。
#[derive(Default)]
struct SignStats {
    n: usize,
    opp: [usize; 3],
    all_opp: usize,
    ge2_opp: usize,
    axis_opp: usize,
    neutral_vs_tinted: usize,
    ours_neutral: usize,
    theirs_tinted: usize,
    theirs_warm: usize,
    ours_warm: usize,
    warm_o: Vec<f64>,
    warm_t: Vec<f64>,
    spread_o: Vec<f64>,
    spread_t: Vec<f64>,
}

fn sign_stats(rows: &[&Row]) -> SignStats {
    let mut s = SignStats::default();
    for r in rows {
        let o = r.ours_linear;
        let t = r.theirs_linear;
        let dev_o = [
            o[0] as f64 - r.lum_ours,
            o[1] as f64 - r.lum_ours,
            o[2] as f64 - r.lum_ours,
        ];
        let dev_t = [
            t[0] as f64 - r.lum_theirs,
            t[1] as f64 - r.lum_theirs,
            t[2] as f64 - r.lum_theirs,
        ];
        let mut n_opp = 0usize;
        for ((slot, a), b) in s.opp.iter_mut().zip(dev_o).zip(dev_t) {
            if a * b < 0.0 {
                *slot += 1;
                n_opp += 1;
            }
        }
        if n_opp == 3 {
            s.all_opp += 1;
        }
        if n_opp >= 2 {
            s.ge2_opp += 1;
        }
        let warm_o = o[0] as f64 - o[2] as f64;
        let warm_t = t[0] as f64 - t[2] as f64;
        if warm_o * warm_t < 0.0 {
            s.axis_opp += 1;
        }
        let so = spread(o);
        let st = spread(t);
        if so < NEUTRAL_SPREAD {
            s.ours_neutral += 1;
        }
        if st > TINTED_SPREAD {
            s.theirs_tinted += 1;
        }
        if so < NEUTRAL_SPREAD && st > TINTED_SPREAD {
            s.neutral_vs_tinted += 1;
        }
        if warm_t > 0.0 {
            s.theirs_warm += 1;
        }
        if warm_o > 0.0 {
            s.ours_warm += 1;
        }
        s.warm_o.push(warm_o);
        s.warm_t.push(warm_t);
        s.spread_o.push(so);
        s.spread_t.push(st);
        s.n += 1;
    }
    s
}

fn pct(x: usize, n: usize) -> f64 {
    100.0 * x as f64 / n.max(1) as f64
}

/// [6] 分组画像：三个关键分组的中位 Lab 与 ΔE00 的明度/色度构成。
///
/// Lab 取「参考」与「**变换后**的我们」——与 ΔE00 同一口径。深阴影里 L* 对 Y 的斜率
/// 极大（Y≈0.007 处 dL*/dY ≈ 1000），且 ProPhoto 的蓝原色对 Y 的贡献几乎为 0，
/// 所以**一个纯色度差异也会落在 L* 轴上**：必须 L*/a*/b* 三项一起看，不能只看 L*。
fn report_groups(named: &[(&str, Vec<&Row>)]) {
    eprintln!("\n[6] 分组画像（中位数；Lab 口径与 ΔE00 一致）");
    for (name, v) in named {
        if v.is_empty() {
            continue;
        }
        let lab_of = |f: &dyn Fn(&Row) -> [f32; 3]| -> [f64; 3] {
            let mut l = Vec::new();
            let mut a = Vec::new();
            let mut b = Vec::new();
            for r in v {
                let lab = f(r);
                l.push(lab[0] as f64);
                a.push(lab[1] as f64);
                b.push(lab[2] as f64);
            }
            [median_of(l), median_of(a), median_of(b)]
        };
        let lo = lab_of(&|r| nikonrawview::deltae::lab_from_prophoto_linear(r.ours_transformed));
        let lt = lab_of(&|r| nikonrawview::deltae::lab_from_prophoto_linear(r.theirs_linear));
        let chroma = |lab: [f64; 3]| (lab[1] * lab[1] + lab[2] * lab[2]).sqrt();
        let de = median_of(v.iter().map(|r| r.de).collect());
        let light = lightness_share_median(v);
        let snr = median_of(v.iter().map(|r| r.local.snr).collect());
        let rel = |f: &dyn Fn(&Row) -> [f32; 3]| -> f64 {
            median_of(
                v.iter()
                    .map(|r| {
                        let lab = f(r);
                        let c = ((lab[1] * lab[1] + lab[2] * lab[2]) as f64).sqrt();
                        c / lab[0].max(1e-6) as f64
                    })
                    .collect(),
            )
        };
        let rel_o = rel(&|r| nikonrawview::deltae::lab_from_prophoto_linear(r.ours_transformed));
        let rel_t = rel(&|r| nikonrawview::deltae::lab_from_prophoto_linear(r.theirs_linear));
        eprintln!("    --- {name}（n={}）---", v.len());
        eprintln!(
            "        我们 L* {:>6.2}  a* {:>+7.2}  b* {:>+7.2}  C* {:>6.2}",
            lo[0],
            lo[1],
            lo[2],
            chroma(lo)
        );
        eprintln!(
            "        参考 L* {:>6.2}  a* {:>+7.2}  b* {:>+7.2}  C* {:>6.2}",
            lt[0],
            lt[1],
            lt[2],
            chroma(lt)
        );
        eprintln!(
            "        ΔL* {:>+6.2}   ΔC* {:>+6.2}   ΔE00 {de:>5.2}（明度项占 ΔE00² 的 {:.0}%）  SNR 中位 {snr:>5.2}  相对彩度 C*/L* 我们/参考 {rel_o:>5.2}/{rel_t:>5.2}",
            lt[0] - lo[0],
            chroma(lt) - chroma(lo),
            100.0 * light
        );
    }
}

/// [7] 方向统计：逐组回答"我们读出中性而参考偏色吗"。
///
/// 注意单通道符号在**近中性**像素上是噪声主导的：我们的通道偏差本来就≈0，符号随机。
/// 判读要看"三通道全相反""暖冷轴相反"与"中位通道极差"三组数，而不是单通道百分比。
fn report_signs(named: &[(&str, Vec<&Row>)]) {
    eprintln!("\n[7] 方向统计：暗部到底是不是「我们中性、参考偏色」");
    for (name, v) in named {
        let s = sign_stats(v);
        let n = s.n.max(1);
        eprintln!(
            "    --- {name}（n={}，中位亮度 我们 {:.4} / 参考 {:.4}）---",
            s.n,
            median_of(v.iter().map(|r| r.lum_ours).collect()),
            median_of(v.iter().map(|r| r.lum_theirs).collect())
        );
        eprintln!(
            "        三通道全相反 {:.1}%   至少两通道相反 {:.1}%   暖冷轴(R−B)符号相反 {:.1}%",
            pct(s.all_opp, n),
            pct(s.ge2_opp, n),
            pct(s.axis_opp, n)
        );
        eprintln!(
            "        我们读出中性(spread<{NEUTRAL_SPREAD}) {:.1}%   参考读出有色(spread>{TINTED_SPREAD}) {:.1}%   我们中性且参考有色 {:.1}%",
            pct(s.ours_neutral, n),
            pct(s.theirs_tinted, n),
            pct(s.neutral_vs_tinted, n)
        );
        eprintln!(
            "        参考偏暖(R>B) {:.1}%   我们偏暖(R>B) {:.1}%   中位 R−B：我们 {:+.4} / 参考 {:+.4}",
            pct(s.theirs_warm, n),
            pct(s.ours_warm, n),
            median_of(s.warm_o),
            median_of(s.warm_t)
        );
        eprintln!(
            "        中位通道极差：我们 {:.4} / 参考 {:.4}；单通道符号相反 R {:.1}% / G {:.1}% / B {:.1}%",
            median_of(s.spread_o),
            median_of(s.spread_t),
            pct(s.opp[0], n),
            pct(s.opp[1], n),
            pct(s.opp[2], n)
        );
    }
}

/// ProPhoto 的 Y（Lab 的明度轴就是它的立方根）。
///
/// **不要**用 Rec.709 权重去近似它：ProPhoto 的蓝原色对 Y 的贡献几乎为 0，而它占
/// Rec.709 亮度的 7.2%。深阴影里色度差一大，两者就显著分叉。
fn pp_y(v: [f32; 3]) -> f64 {
    (0.288_040_2 * v[0] as f64 + 0.711_874_1 * v[1] as f64 + 0.000_085_7 * v[2] as f64).max(0.0)
}

/// [11] 分阶段诊断：恒等 / 仅曲线 / 曲线+LUT 在同一评估集上的 ΔE00 与中位 Y。
///
/// 它回答"噪声之外，误差是哪一级造成的"——深阴影若在**仅曲线**时就已经对不上，
/// 问题在曲线的定义域/分箱；若是加上 LUT 之后才崩，问题在 3D LUT 的网格分辨率。
fn report_stages(named: &[(&str, Vec<&Row>)], t: &BaseTransform) {
    let curve_only = BaseTransform::from_fit(0, "仅曲线", Vec::new(), t.curve.clone(), 0, Vec::new());
    eprintln!("\n[11] 变换分阶段：误差是哪一级造成的（同一评估集，中位数）");
    eprintln!("    组                恒等ΔE00  仅曲线ΔE00  曲线+LUTΔE00 |  中位 Y：原始 / 仅曲线 / 曲线+LUT / 参考");
    for (name, v) in named {
        if v.is_empty() {
            continue;
        }
        let med = |f: &dyn Fn(&Row) -> f64| median_of(v.iter().map(|r| f(r)).collect());
        let de0 = med(&|r| nikonrawview::deltae::delta_e_prophoto(r.ours_linear, r.theirs_linear));
        let dec = med(&|r| {
            nikonrawview::deltae::delta_e_prophoto(curve_only.apply(r.ours_linear), r.theirs_linear)
        });
        let def = med(&|r| r.de);
        let y0 = med(&|r| pp_y(r.ours_linear));
        let yc = med(&|r| pp_y(curve_only.apply(r.ours_linear)));
        let yf = med(&|r| pp_y(r.ours_transformed));
        let yt = med(&|r| pp_y(r.theirs_linear));
        eprintln!(
            "    {name:<13}  {de0:>8.2}  {dec:>10.2}  {def:>11.2} |  {y0:.4} / {yc:.4} / {yf:.4} / {yt:.4}"
        );
    }
}

/// [12] 暗部传递曲线：按**原始输入 Y** 分箱，看曲线、曲线+LUT、参考各自把 Y 送到哪里。
///
/// 这是判"暗部对比度是否被压平"的直接证据：把每箱的四个中位 Y 排在一起，
/// 一眼就能看出传递曲线在暗部是抬了、压了还是平了。
fn report_shadow_transfer(rows: &[&Row], t: &BaseTransform) {
    let curve_only = BaseTransform::from_fit(0, "仅曲线", Vec::new(), t.curve.clone(), 0, Vec::new());
    const EDGES: &[f64] = &[0.0, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128];
    eprintln!("\n[12] 暗部传递曲线（按原始输入 Y 分箱，均为中位数）");
    eprintln!("    输入 Y 区间            n     Y原始    Y仅曲线  Y曲线+LUT   Y参考    曲线+LUT/参考");
    for w in EDGES.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        let sel: Vec<&Row> = rows
            .iter()
            .copied()
            .filter(|r| (lo..hi).contains(&pp_y(r.ours_linear)))
            .collect();
        if sel.is_empty() {
            continue;
        }
        let med = |f: &dyn Fn(&Row) -> f64| median_of(sel.iter().map(|r| f(r)).collect());
        let y0 = med(&|r| pp_y(r.ours_linear));
        let yc = med(&|r| pp_y(curve_only.apply(r.ours_linear)));
        let yf = med(&|r| pp_y(r.ours_transformed));
        let yt = med(&|r| pp_y(r.theirs_linear));
        eprintln!(
            "    [{lo:.4}, {hi:.4})  {:>7}  {y0:.5}  {yc:.5}  {yf:.5}  {yt:.5}  {:>8.2}×",
            sel.len(),
            ratio_or_inf(yf, yt)
        );
    }
}

/// [8] 误差最大的若干样本。
fn report_worst(rows: &[Row], n: usize) {
    let mut by_de: Vec<usize> = (0..rows.len()).collect();
    by_de.sort_by(|&a, &b| cmp_f64(rows[b].de, rows[a].de));
    eprintln!("\n[8] 误差最大的 {n} 个样本：它们的局部信噪比是多少");
    eprintln!("    序号     xy          ΔE00    SNR   σ_hf      |Δ亮度|/σ_hf  我们线性 → 参考线性");
    for &i in by_de.iter().take(n) {
        let r = &rows[i];
        let dnr = ratio_or_inf((r.lum_theirs - r.lum_ours).abs(), r.local.hf_sd);
        eprintln!(
            "    #{:<6} ({:>4},{:>4})  {:>6.2}  {:>5.2}  {:.4}   {:>9.1}     [{:.3} {:.3} {:.3}] → [{:.3} {:.3} {:.3}]",
            r.idx,
            r.xy.0,
            r.xy.1,
            r.de,
            r.local.snr,
            r.local.hf_sd,
            dnr,
            r.ours_linear[0],
            r.ours_linear[1],
            r.ours_linear[2],
            r.theirs_linear[0],
            r.theirs_linear[1],
            r.theirs_linear[2],
        );
    }
}

// ---------------------------------------------------------------------------
// [9] 拍摄参数（尽力读取）：用于判断"噪声由拍摄条件决定"这一解释是否成立
// ---------------------------------------------------------------------------

/// Exif 前缀读取窗口——ExifIFD 与 MakerNote 都远在前 4 MiB 内。
const EXIF_PREFIX: u64 = 4 * 1024 * 1024;

fn rational_of(t: &Tiff, ifd: usize, tag: u16) -> Option<f64> {
    let e = t.find(ifd, tag).ok()??;
    let b = t.bytes(&e).ok()?;
    if b.len() < 8 {
        return None;
    }
    let num = t.order.u32(&b[0..4]) as f64;
    let den = t.order.u32(&b[4..8]) as f64;
    (den != 0.0).then_some(num / den)
}

fn ascii_of(t: &Tiff, ifd: usize, tag: u16) -> Option<String> {
    let e = t.find(ifd, tag).ok()??;
    let b = t.bytes(&e).ok()?;
    let s = String::from_utf8_lossy(b).trim_end_matches('\0').trim().to_string();
    (!s.is_empty()).then_some(s)
}

/// 从一份（可能是文件前缀的）NEF 字节里读 ISO / 曝光时间 / 光圈 / 拍摄时间。
fn parse_exposure(data: &[u8]) -> Option<(u32, f64, f64, String)> {
    let t = Tiff::locate(data).ok()?;
    let ifd0 = t.ifd0_offset().ok()?;
    let exif = t.find(ifd0, 0x8769).ok()??;
    let exif_abs = t.deref_offset(&exif).ok()?;
    let iso = t.find(exif_abs, 0x8827).ok().flatten().and_then(|e| t.u16_value(&e).ok())? as u32;
    Some((
        iso,
        rational_of(&t, exif_abs, 0x829A).unwrap_or(f64::NAN),
        rational_of(&t, exif_abs, 0x829D).unwrap_or(f64::NAN),
        ascii_of(&t, exif_abs, 0x9003).unwrap_or_else(|| "?".into()),
    ))
}

/// 先读前缀，读不到再退回整份文件（50 MB 一次，只在这几个样本上发生）。
fn exposure_of(nef: &Path) -> Option<(u32, f64, f64, String)> {
    use std::io::Read;
    if let Ok(f) = std::fs::File::open(nef) {
        let mut src = f.take(EXIF_PREFIX);
        let mut head = Vec::new();
        if src.read_to_end(&mut head).is_ok() {
            if let Some(v) = parse_exposure(&head) {
                return Some(v);
            }
        }
    }
    std::fs::read(nef).ok().and_then(|d| parse_exposure(&d))
}

fn report_exposure() {
    let dir = samples_dir();
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x.eq_ignore_ascii_case("nef")))
        .collect();
    files.sort();
    if files.is_empty() {
        return;
    }
    eprintln!("\n[9] 拍摄参数对照（Exif 尽力读取）——用于判断「噪声随拍摄条件变化」能否解释 P95 恒定");
    let mut isos: Vec<u32> = Vec::new();
    let mut unreadable = 0usize;
    let shown = files.len().min(20);
    for p in files.iter().take(shown) {
        let name = p.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        match exposure_of(p) {
            Some((iso, et, fnum, when)) => {
                isos.push(iso);
                let expo = if et.is_finite() && et > 0.0 && et < 1.0 {
                    format!("1/{:.0}s", 1.0 / et)
                } else {
                    format!("{et:.3}s")
                };
                eprintln!("    {name:<14} ISO {iso:<5} 曝光 {expo:>9}  光圈 f/{fnum:<4.1}  {when}");
            }
            None => {
                unreadable += 1;
                eprintln!("    {name:<14} （Exif 读不到）");
            }
        }
    }
    isos.sort_unstable();
    isos.dedup();
    eprintln!(
        "    共 {} 个 NEF（列出前 {shown} 个），可读 {} 个、读不到 {unreadable} 个；ISO 取值 {isos:?}",
        files.len(),
        shown - unreadable
    );
    if isos.len() <= 1 {
        eprintln!(
            "    ⇒ ISO 一致，「不同照片噪声水平不同」不成立——P95 恒定与「噪声决定误差」并不矛盾"
        );
    } else {
        let stops = (*isos.last().unwrap() as f64 / isos[0] as f64).log2();
        eprintln!(
            "    ⇒ ISO 跨度 {}–{}（约 {stops:.1} 档）⇒ 噪声水平并不相同；若大误差主要来自噪声，四张图的 P95 不该几乎恒定",
            isos[0],
            isos.last().unwrap()
        );
    }
}

// ---------------------------------------------------------------------------
// 用例
// ---------------------------------------------------------------------------

#[test]
fn shadow_residual_noise_or_bias() {
    let Some(fx) = fitted_neutral_fixture() else {
        return;
    };
    let nef = samples_dir().join("DSC_0001.NEF");
    eprintln!("=== 诊断 C：暗部大误差像素 —— 噪声还是系统性偏差？ ===");
    eprintln!(
        "有效区 {}×{}，拟合样本 {}，留出评估样本 {}",
        fx.size.0,
        fx.size.1,
        fx.fit_count,
        fx.eval.len()
    );

    // ---- [0] 基线：先复现未达标的结论 ----
    let mut des: Vec<f64> = fx.eval.iter().map(|s| s.delta_e()).collect();
    let (med, p95, max) = common::quantiles(&mut des);
    eprintln!("\n[0] 基线：留出评估集 ΔE00 中位 {med:.3}（目标 ≤1.0）、P95 {p95:.3}（目标 ≤3.0）、最大 {max:.3}");

    // ---- [1] 再解码一次，取每个样本真实的 3×3 邻域 ----
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let dec = decode_working_space(&nef, &opts).expect("密集解码应成功（脚手架刚用同一参数解过）");
    let model = nikonrawview::camera::read_model(&nef).expect("应能读出机型");
    let entry = nikonrawview::camera::lookup(&model).expect("机型应在边距表中");
    let mg = if dec.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = mg.effective(dec.width, dec.height).expect("边距应能裁出有效区");
    assert_eq!((cw, ch), fx.size, "两次解码的有效区尺寸应一致");

    let mut rows: Vec<Row> = Vec::with_capacity(fx.eval.len());
    let mut missing = 0usize;
    let mut coord_err = 0.0f64;
    for (idx, s) in fx.eval.iter().enumerate() {
        let (x, y) = s.xy;
        let Some(local) = local_stats(&dec, mg.left + x, mg.top + y) else {
            missing += 1;
            continue;
        };
        // 坐标自检：密集解码里中心像素的亮度，必须与脚手架给的样本逐位对上。
        coord_err = coord_err.max((local.lum_center - s.lum_ours() as f64).abs());
        rows.push(Row {
            idx,
            xy: s.xy,
            de: s.delta_e(),
            ours_linear: s.ours_linear,
            ours_transformed: s.ours_transformed,
            theirs_linear: s.theirs_linear,
            lum_ours: s.lum_ours() as f64,
            lum_theirs: s.lum_theirs() as f64,
            local,
        });
    }
    eprintln!(
        "\n[1] 密集解码 {}×{}（旋转={}），边距 左{} 上{} 右{} 下{}；{} 个样本取到 3×3 邻域，{} 个越界被排除",
        dec.width, dec.height, dec.rotated, mg.left, mg.top, mg.right, mg.bottom, rows.len(), missing
    );
    eprintln!("    坐标自检：中心像素亮度与脚手架样本的最大偏差 = {coord_err:.3e}（应为 0）");
    assert!(
        rows.len() * 100 >= fx.eval.len() * 99,
        "3×3 邻域应当几乎处处可得，实际只取到 {}/{}",
        rows.len(),
        fx.eval.len()
    );
    assert!(rows.len() > 10_000, "样本太少（{}），统计不可信", rows.len());
    assert!(
        coord_err < 1e-6,
        "坐标映射不一致（最大偏差 {coord_err:.3e}）：密集解码与脚手架取的不是同一批像素，本诊断的邻域统计无效"
    );

    // 暗部 10% 与中间调：先按我们的亮度排序。
    let dark_n = (rows.len() / 10).max(1);
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| cmp_f64(rows[a].lum_ours, rows[b].lum_ours));
    let dark: Vec<&Row> = order[..dark_n].iter().map(|&i| &rows[i]).collect();
    let mid: Vec<&Row> = rows
        .iter()
        .filter(|r| r.lum_ours > 0.10 && r.lum_ours <= 0.60)
        .collect();
    eprintln!(
        "    暗部 10% = 我们亮度最低的 {dark_n} 个样本，亮度区间 [{:.4}, {:.4}]；中间调（0.10<亮度≤0.60）{} 个",
        dark[0].lum_ours,
        dark[dark.len() - 1].lum_ours,
        mid.len()
    );

    report_bins(&rows);
    report_correlation(&rows, &dark, &mid);
    report_top_profile(&rows, &dark);
    report_noise_vs_signal(&rows, &dark, &fx.transform);

    // 三个分组对照：最暗 10%、ΔE00 最高 10%、全体。
    let by_de = order_by_de(&rows);
    let top_n = (rows.len() / 10).max(1);
    let worst10: Vec<&Row> = by_de[..top_n].iter().map(|&i| &rows[i]).collect();
    let all: Vec<&Row> = rows.iter().collect();
    let groups: Vec<(&str, Vec<&Row>)> = vec![
        ("暗部最暗10%", dark.clone()),
        ("ΔE00最高10%", worst10),
        ("全体", all),
    ];
    report_groups(&groups);
    report_signs(&groups);
    report_stages(&groups, &fx.transform);
    report_shadow_transfer(&groups[2].1, &fx.transform);
    report_worst(&rows, 10);
    report_exposure();

    // ---- [10] 结论（全部由上面的数字算出来，不是写死的） ----
    let worst_snr = median_of(by_de[..top_n].iter().map(|&i| rows[i].local.snr).collect());
    let rest_snr = median_of(by_de[top_n..].iter().map(|&i| rows[i].local.snr).collect());
    let snrs: Vec<f64> = rows.iter().filter(|r| r.local.snr.is_finite()).map(|r| r.local.snr).collect();
    let des_all: Vec<f64> = rows.iter().filter(|r| r.local.snr.is_finite()).map(|r| r.de).collect();
    let rho = spearman(&snrs, &des_all);
    let d_snr: Vec<f64> = dark.iter().filter(|r| r.local.snr.is_finite()).map(|r| r.local.snr).collect();
    let d_de: Vec<f64> = dark.iter().filter(|r| r.local.snr.is_finite()).map(|r| r.de).collect();
    let rho_dark = spearman(&d_snr, &d_de);
    let dnr_dark = median_of(
        dark.iter()
            .map(|r| ratio_or_inf((r.lum_theirs - r.lum_ours).abs(), r.local.hf_sd))
            .collect(),
    );
    let foot_dark = median_of(dark.iter().map(|r| noise_footprint(&fx.transform, r, 3.0)).collect());
    let de_dark = median_of(dark.iter().map(|r| r.de).collect());
    let lab_l = |f: &dyn Fn(&Row) -> [f32; 3]| {
        median_of(dark.iter().map(|r| f(r)[0] as f64).collect())
    };
    let l_ours = lab_l(&|r| nikonrawview::deltae::lab_from_prophoto_linear(r.ours_transformed));
    let l_theirs = lab_l(&|r| nikonrawview::deltae::lab_from_prophoto_linear(r.theirs_linear));
    let light_dark = lightness_share_median(&dark);
    let sd = sign_stats(&dark);

    eprintln!("\n[10] 结论");
    eprintln!(
        "    1) 高 ΔE00 是否伴随低 SNR：最高 10% 的 SNR 中位 {worst_snr:.2}，其余 {rest_snr:.2}；\
         全体 ρ(SNR, ΔE00) = {rho:.3}，暗部内部 ρ = {rho_dark:.3}"
    );
    eprintln!(
        "       ⇒ {}",
        if rho_dark < -0.3 {
            "是，即使在暗部内部，SNR 越低误差也越大"
        } else if rho_dark < -0.1 {
            "整体上有负相关，但**控制住亮度之后**只剩弱相关——高误差主要是「暗」而不是「噪」"
        } else {
            "不是：暗部内部 SNR 与 ΔE00 基本无关，高误差不是「低信噪比」的直接后果"
        }
    );
    eprintln!(
        "    2) 差异是否在噪声之上：暗部 |Δ亮度|/σ_hf 中位 = {dnr_dark:.2}；\
         沿误差方向移 3σ 的噪声足迹 ΔE00 = {foot_dark:.2}，而实际 ΔE00 中位 = {de_dark:.2} → 实际/3σ = {:.1} 倍",
        ratio_or_inf(de_dark, foot_dark)
    );
    eprintln!(
        "       ⇒ {}",
        if ratio_or_inf(de_dark, foot_dark) > 3.0 {
            "远超噪声：这批暗部像素的误差不可能由局部噪声解释，必须被匹配"
        } else {
            "与噪声同量级：不能排除不可匹配"
        }
    );
    let dl = l_ours - l_theirs; // 正 = 变换后我们偏亮
    eprintln!(
        "    3) 误差构成：暗部（**变换后**）L* 我们 {l_ours:.2} / 参考 {l_theirs:.2}（ΔL* {dl:+.2}，正=我们偏亮），\
         明度项占 ΔE00² 的 {:.0}%",
        100.0 * light_dark
    );
    eprintln!(
        "       ⇒ {}",
        if dl > 3.0 {
            "深阴影被**抬得比参考亮**，是变换在深阴影端的系统性偏差，不是噪声；看 [11] 判断该归咎曲线还是 LUT"
        } else if dl < -3.0 {
            "深阴影比参考**暗**，同样是系统性偏差，不是噪声"
        } else {
            "明暗差不大，误差主体在色度侧"
        }
    );
    eprintln!(
        "    4) 方向是否随机：暗部三通道全相反 {:.1}%、暖冷轴相反 {:.1}%、参考偏暖 {:.1}%",
        pct(sd.all_opp, sd.n),
        pct(sd.axis_opp, sd.n),
        pct(sd.theirs_warm, sd.n)
    );

    assert!(med.is_finite() && p95.is_finite() && max.is_finite(), "基线 ΔE00 应可计算");
    assert!(!worst_snr.is_nan() && !rest_snr.is_nan(), "SNR 中位不应为 NaN");
    assert_eq!(sd.n, dark_n, "符号统计应覆盖整个暗部十分位");
}
