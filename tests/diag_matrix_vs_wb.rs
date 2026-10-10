//! 诊断 G：**逐图最优 3×3 是否随拍摄白平衡系统变化**。
//!
//! # 要检验的假设
//!
//! 「我们用的相机矩阵（`camera::color_matrix("Nikon Z 8")`）与尼康的等效色彩变换不是同一个，
//! 且偏离程度随拍摄白平衡的极端程度增大。」
//!
//! 机制上的说法是：Adobe DNG 与尼康按「两个已知光源各一个矩阵」存标定、依拍摄白平衡插值，
//! 而 LibRaw 的 `adobe_coeff` 表对每台相机只存一个矩阵、不插值。于是白平衡离标定光源越远，
//! LibRaw 越偏。
//!
//! 本诊断**不检验该机制本身**（那要读解码器源码），只检验它可观测的推论：逐图解出的
//! 「相机线性 → 参考 ProPhoto 线性」最优变换，与我们的固定矩阵之间的**功能**偏离，
//! 是否随拍摄白平衡（`B/G`）变化。
//!
//! # 采集
//!
//! 1. 白平衡取 `libraw::read_wb` 的 `cam_mul`，以 `B/G`、`R/G` 表征——这是**相机记录的**
//!    原始白平衡，不是 NX Studio 边车里可能改过的那份。
//! 2. 自变量 `c` 取**原始解码**（`OutputColor::Camera`，管线里的 ①），参考 `t` 取脚手架
//!    `fitted_fixture_for` 给的 `theirs_linear`。注意 `c` 里**已经含白平衡**
//!    （① 是 `scale_colors` 之后的值），所以逐图最优变换会把「还需要多少白平衡型修正」
//!    一并吸收进去——这正是要量的东西。
//! 3. **空间留出**：画面左半拟合、右半评估；再反向拟合一次，两半之差即同图噪声底。
//!    样本排除任一侧的饱和/压黑点（`c` 或 `t` 的任一分量 ≤0 或 ≥1）——那里解码层与参考
//!    都做过裁切，会把裁切行为误当成色彩关系。
//!
//! # 两把尺子，以及为什么**只有第二把能用**
//!
//! - **系数尺（`‖ΔM‖F/‖M0‖F`）—— 不可用，本诊断只用它做反面例证。**
//!   单张画面在相机空间的样本分布**近乎退化**（近中性画面 → `c ∝ (1,1,1)` 占主导，
//!   正规矩阵 `Σccᵀ` 病态），因此「最优 3×3」的 9 个系数**不是可识别的量**，实测系数能到
//!   `2.34 / −3.39 / 2.80` 这种非物理值。实测证据：系数离 `M0` 5.5 倍远，功能上却只改善
//!   14%；同一张图两半各自解出的系数也能差 0.75。**「矩阵差多少」在这种数据上量不出来，
//!   只能量「函数差多少」。**
//! - **功能尺（尺度不变色度误差）—— 判据用它。** 对每个样本算
//!   `Δlog2(R/G)`、`Δlog2(B/G)`（预测 vs 参考），取向量长度的中位数。两图各自乘任意正常数
//!   都不改变它，因此**整体增益与曝光/色调的尺度差被完全排除**；逐通道的差由对角模型
//!   `M0·diag(w)` 吸收。这正是诊断 E 用的那把「尺度不变、无亮度污染」的尺子。
//!   两个被拟合的模型（`M0·diag(w)`、`M_opt`）都是按**线性 rms** 拟合的，色度尺只是给它们的
//!   留出预测换一把尺子打分；因此「`M_opt` 不比对角模型好」的含义是**多出来的形状自由度没有**
//!   **换来色度精度**，而不是「色度意义上最优的 3×3 也只能做到这样」。
//!
//! # 四模型阶梯（同一份留出集，都只用左半拟合）
//!
//! | 模型 | 参数数 | 含义 |
//! |---|---|---|
//! | `M0` | 0 | 我们的固定矩阵原样 |
//! | `k·M0` | 1 | 只允许整体增益（色调/曝光尺度） |
//! | `M0·diag(w)` | 3 | 只允许「矩阵之前的三通道增益」= 白平衡型修正 |
//! | `M_opt` | 9 | 无约束 3×3（逐图最优） |
//!
//! 四者**互相嵌套**，所以留出读数的落差可直接读成「哪一级自由度吃掉了偏离」。
//! 解码器里 `scale_colors` 就在 `cam_xyz` 之前，因此 3 参数那一级正是「白平衡型修正」
//! 在本管线里的对应物。在尺度不变的色度尺上，`M0` 与 `k·M0` 必然给出同一个读数
//! （比值对整体增益不变），所以色度表只列 `M0 → M0·diag(w) → M_opt` 三级。
//!
//! # 白平衡型修正 `w` 与跨图迁移
//!
//! `w` 只有 3 个参数、良态、且直接可读（`w_B/w_G` 就是「参考相对我们还需要多少蓝通道增益」）。
//! 因此除了逐图看它，还做一件更贴近假设的事：**把 A 图解出的 `w_A` 用到 B 图的留出半幅上**，
//! 比较「A、B 同档白平衡」与「A、B 不同档」两种迁移的色度误差。若所需修正真的随拍摄白平衡
//! 变化，跨档迁移应当明显差于同档迁移。
//!
//! # 一条必须随结论一起说的限制
//!
//! `simple/` 里那批连拍白平衡几乎相同，真正不同的白平衡只有几档（见 [0] 普查）。
//! **三个点不足以确认「单调趋势」**——三点取值在随机顺序下恰好单调的先验概率是 1/3。
//! 本诊断只回答「是否随 `B/G` 变化」，并且把**档内离散**（同档白平衡的多张图 × 两个半幅估计）
//! 当作噪声底：档间落差没有超出噪声底时，如实报「不可辨」。读数门槛先定后看
//! （见 `SNR_WEAK`/`SNR_STRONG`），不在看到数字之后挑口径。
//!
//! # 其它限制
//!
//! - 只对 7 张图做完整配对（每张要读 260 MB 参考 TIF + 两次解码，约 1 分钟）；其余 NEF
//!   只做白平衡普查（`unpack`，不解马赛克），用于确认「只有几档白平衡」。
//! - `M_opt` 是**逐图**解出来的，会吸收该图的一切线性差异，因此它**不是物理相机矩阵**，
//!   只能当「线性上界」用。
//! - 色度尺对**逐通道**的曲线差不是免疫的（只对整体增益免疫）：参考导出带 NEUTRAL 的逐通道
//!   色调曲线，那一部分会被算进「对角吃不掉的残差」里。故「形状份额」是**上界**。
//! - 环境变量：`NRV_DIAG_STEM=DSC_0001,DSC_0141` 只跑指定图（迭代用）；
//!   `NRV_DIAG_CENSUS=0` 跳过白平衡普查。

mod common;

use common::{fitted_fixture_for, quantiles, samples_dir};
use nikonrawview::camera;
use nikonrawview::color;
use nikonrawview::deltae;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::mat3::{self, Mat3};
use std::time::Instant;

// ---------------------------------------------------------------------------
// 口径与门槛（先定后看）
// ---------------------------------------------------------------------------

/// 默认参与矩阵求解的图。
///
/// 前 5 张同属那批 42 秒连拍（白平衡几乎相同，用来量档内离散），后 2 张是另外两档白平衡。
const DEFAULT_STEMS: [&str; 7] = [
    "DSC_0001", "DSC_0002", "DSC_0003", "DSC_0004", "DSC_0010", "DSC_0141", "DSC_8562",
];

/// 饱和/压黑门槛（开区间）：任一侧的任一分量触到就丢弃该样本。
const GATE_LO: f32 = 0.0;
const GATE_HI: f32 = 1.0;

/// 色度比值的可用门槛：任一通道低于它就谈不了 `R/G`、`B/G`（16 位量化下 1 个计数是 1/65535）。
const RATIO_GATE: f32 = 1e-4;

/// 单个半幅可用的最少配对数，低于此值认为样本不足。
const MIN_PAIRS: usize = 1000;

/// 视作同一档白平衡的 `B/G` 容差。
const WB_GROUP_TOL: f64 = 2e-3;

/// 信噪比读数门槛：档间极差 ÷ 档内合并标准差。
const SNR_WEAK: f64 = 2.0;
const SNR_STRONG: f64 = 5.0;

/// 四模型阶梯的名字，供打印。
const MODELS: [&str; 4] = ["M0", "k·M0", "M0·diag(w)", "M_opt(3×3)"];

// ---------------------------------------------------------------------------
// 基础工具
// ---------------------------------------------------------------------------

/// 解码层的 16 位输出 → 0..1。与脚手架 `fitted_fixture_for` 的换算保持一致。
fn norm96(p: [u16; 3]) -> [f32; 3] {
    [p[0] as f32 / 65535.0, p[1] as f32 / 65535.0, p[2] as f32 / 65535.0]
}

/// 逐分量最大绝对差。
fn max_abs_diff(a: Mat3, b: Mat3) -> f64 {
    let mut m = 0f64;
    for (ra, rb) in a.iter().zip(b.iter()) {
        for (x, y) in ra.iter().zip(rb.iter()) {
            m = m.max((x - y).abs() as f64);
        }
    }
    m
}

fn frob(m: Mat3) -> f64 {
    let mut s = 0f64;
    for row in m.iter() {
        for v in row.iter() {
            s += (*v as f64).powi(2);
        }
    }
    s.sqrt()
}

/// 整体缩放。
fn scaled(m: Mat3, k: f32) -> Mat3 {
    let mut o = m;
    for row in o.iter_mut() {
        for v in row.iter_mut() {
            *v *= k;
        }
    }
    o
}

fn fmt_rows(m: Mat3) -> [String; 3] {
    [
        format!("[{:>9.5} {:>9.5} {:>9.5}]", m[0][0], m[0][1], m[0][2]),
        format!("[{:>9.5} {:>9.5} {:>9.5}]", m[1][0], m[1][1], m[1][2]),
        format!("[{:>9.5} {:>9.5} {:>9.5}]", m[2][0], m[2][1], m[2][2]),
    ]
}

/// 线性残差读数（尺子一：与参考同为 0..1 刻度，受整体尺度主导）。
#[derive(Debug, Clone, Copy)]
struct Score {
    rms: f64,
    median: f64,
    p95: f64,
}

/// 尺度不变的色度读数（尺子二，**判据**）。
#[derive(Debug, Clone, Copy)]
struct RatioErr {
    /// `|(Δlog2(R/G), Δlog2(B/G))|` 的中位数。
    median: f64,
    p95: f64,
    /// 两个分量的**带符号**中位数，用来看偏差方向是否系统性。
    med_rg: f64,
    med_bg: f64,
    n: usize,
}

impl RatioErr {
    fn empty() -> Self {
        Self { median: f64::NAN, p95: f64::NAN, med_rg: f64::NAN, med_bg: f64::NAN, n: 0 }
    }
}

/// 用给定的前端 `f` 在一批样本上算 rms 与 ΔE00 分位数。
fn score<F: Fn([f32; 3]) -> [f32; 3]>(pairs: &[([f32; 3], [f32; 3])], f: F) -> Score {
    let mut sum_sq = 0f64;
    let mut n = 0usize;
    let mut des = Vec::with_capacity(pairs.len());
    for (c, t) in pairs {
        let p = f(*c);
        for k in 0..3 {
            let e = (p[k] - t[k]) as f64;
            sum_sq += e * e;
            n += 1;
        }
        des.push(deltae::delta_e_prophoto(p, *t));
    }
    let (median, p95, _max) = quantiles(&mut des);
    Score { rms: (sum_sq / n.max(1) as f64).sqrt(), median, p95 }
}

/// 三个通道都有限且超过 [`RATIO_GATE`]——只有这样才能谈 `R/G`、`B/G`。
fn usable_ratio(v: &[f32; 3]) -> bool {
    v.iter().all(|x| x.is_finite() && *x > RATIO_GATE)
}

/// 尺度不变的色度误差：`(log2 R/G, log2 B/G)` 的「预测 − 参考」，取向量长度的中位数。
///
/// 预测与参考各自乘任意正常数都不改变这两个分量，因此**整体增益、曝光尺度被完全排除**。
/// 只统计两边的三个通道都超过 [`RATIO_GATE`] 的样本。
fn ratio_error<F: Fn([f32; 3]) -> [f32; 3]>(pairs: &[([f32; 3], [f32; 3])], f: F) -> RatioErr {
    let mut lens = Vec::new();
    let mut drg = Vec::new();
    let mut dbg = Vec::new();
    for (c, t) in pairs {
        let p = f(*c);
        if !usable_ratio(&p) || !usable_ratio(t) {
            continue;
        }
        let a = ((p[0] / p[1]) as f64).log2() - ((t[0] / t[1]) as f64).log2();
        let b = ((p[2] / p[1]) as f64).log2() - ((t[2] / t[1]) as f64).log2();
        drg.push(a);
        dbg.push(b);
        lens.push((a * a + b * b).sqrt());
    }
    if lens.is_empty() {
        return RatioErr::empty();
    }
    let n = lens.len();
    let (median, p95, _max) = quantiles(&mut lens);
    let (med_rg, _, _) = quantiles(&mut drg);
    let (med_bg, _, _) = quantiles(&mut dbg);
    RatioErr { median, p95, med_rg, med_bg, n }
}

/// 静态的矩阵比较读数：`(相对 Frobenius 偏离, 矩阵空间最优增益, 去增益后的相对形状残差)`。
///
/// **注意**：在近中性的画面样本上这些系数不可识别，见模块文档「两把尺子」。本函数只用于
/// 展示「系数尺为什么不能用」。
fn matrix_metrics(m: Mat3, m0: Mat3) -> (f64, f64, f64) {
    let (mut num, mut den) = (0f64, 0f64);
    for i in 0..3 {
        for j in 0..3 {
            num += m[i][j] as f64 * m0[i][j] as f64;
            den += (m0[i][j] as f64).powi(2);
        }
    }
    let k = if den > 0.0 { num / den } else { 1.0 };
    let mut diff = 0f64;
    let mut shape = 0f64;
    for i in 0..3 {
        for j in 0..3 {
            diff += (m[i][j] as f64 - m0[i][j] as f64).powi(2);
            shape += (m[i][j] as f64 - k * m0[i][j] as f64).powi(2);
        }
    }
    let base = den.sqrt().max(1e-12);
    (diff.sqrt() / base, k, shape.sqrt() / base)
}

/// 在 `pairs` 上求让 `k·M0·c` 最接近 `t` 的标量 `k`（一元最小二乘）。
fn best_gain(m0: Mat3, pairs: &[([f32; 3], [f32; 3])]) -> f64 {
    let (mut num, mut den) = (0f64, 0f64);
    for (c, t) in pairs {
        let p = mat3::mul_vec(m0, *c);
        for k in 0..3 {
            num += p[k] as f64 * t[k] as f64;
            den += (p[k] as f64).powi(2);
        }
    }
    if den > 0.0 {
        num / den
    } else {
        1.0
    }
}

fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// 三元线性方程组，Cramer 法则；退化时返回 `None`。
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-18 {
        return None;
    }
    let mut out = [0f64; 3];
    for i in 0..3 {
        let mut m = a;
        for r in 0..3 {
            m[r][i] = b[r];
        }
        out[i] = det3(m) / d;
    }
    Some(out)
}

/// 在 `pairs` 上求 `w` 使 `M0·diag(w)·c` 最接近 `t`。
///
/// 预测对 `w` 是线性的：`M0·diag(w)·c = Σ_j w_j · u_j`，其中 `u_j = M0 的第 j 列 × c_j`。
/// 于是仍是三元最小二乘，只需一次 3×3 求解——**良态**，且 `w` 直接可读。
///
/// 这就是「白平衡型修正」的形式：解码器里 `scale_colors` 正是在 `cam_xyz` 之前按通道乘系数。
fn fit_diag_gain(m0: Mat3, pairs: &[([f32; 3], [f32; 3])]) -> Option<[f64; 3]> {
    let mut a = [[0f64; 3]; 3];
    let mut b = [0f64; 3];
    for (c, t) in pairs {
        // u[j][k] = M0[k][j] · c_j
        let mut u = [[0f64; 3]; 3];
        for j in 0..3 {
            for k in 0..3 {
                u[j][k] = m0[k][j] as f64 * c[j] as f64;
            }
        }
        for r in 0..3 {
            for s in 0..3 {
                a[r][s] += (0..3).map(|k| u[r][k] * u[s][k]).sum::<f64>();
            }
            b[r] += (0..3).map(|k| u[r][k] * t[k] as f64).sum::<f64>();
        }
    }
    solve3(a, b)
}

/// 按 `w` 施加白平衡型修正后过 `M0`。
fn apply_diag(m0: Mat3, w: [f64; 3], c: [f32; 3]) -> [f32; 3] {
    mat3::mul_vec(
        m0,
        [
            c[0] * w[0] as f32,
            c[1] * w[1] as f32,
            c[2] * w[2] as f32,
        ],
    )
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().min(ys.len());
    if n < 2 {
        return f64::NAN;
    }
    let mx = xs[..n].iter().sum::<f64>() / n as f64;
    let my = ys[..n].iter().sum::<f64>() / n as f64;
    let (mut sxy, mut sxx, mut syy) = (0f64, 0f64, 0f64);
    for (x, y) in xs[..n].iter().zip(ys[..n].iter()) {
        let dx = x - mx;
        let dy = y - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return f64::NAN;
    }
    sxy / (sxx * syy).sqrt()
}

/// 读数：信噪比 → 判定词。
fn snr_verdict(snr: f64) -> &'static str {
    if !snr.is_finite() {
        "样本不足"
    } else if snr < SNR_WEAK {
        "不可辨（未超出噪声底）"
    } else if snr < SNR_STRONG {
        "弱证据"
    } else {
        "明确变化"
    }
}

fn median_of(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[v.len() / 2]
}

/// `(中位数, 均值)`，不改动调用方的数据。
fn med_mean(v: &[f64]) -> (f64, f64) {
    let mut c = v.to_vec();
    (median_of(&mut c), mean_of(v))
}

fn mean_of(v: &[f64]) -> f64 {
    if v.is_empty() {
        f64::NAN
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

// ---------------------------------------------------------------------------
// 配对采集
// ---------------------------------------------------------------------------

type Pairs = Vec<([f32; 3], [f32; 3])>;

/// 配对采集的丢弃统计。
#[derive(Debug, Default, Clone, Copy)]
struct DropCount {
    missing: usize,
    cam_sat: usize,
    ref_sat: usize,
}

/// 把脚手架的评估样本（裁剪后坐标 + 参考值）与**我们自己的**相机空间解码配对，
/// 并按画面左右一分为二（左半拟合、右半评估）。
fn collect_pairs(
    fx: &common::FittedFixture,
    cam: &libraw::Decoded,
    margins: camera::Margins,
    cw: usize,
) -> (Pairs, Pairs, DropCount) {
    let mut left: Pairs = Vec::new();
    let mut right: Pairs = Vec::new();
    let mut drop = DropCount::default();
    for s in &fx.eval {
        let (x, y) = s.xy;
        let Some(p) = cam.at(x + margins.left, y + margins.top) else {
            drop.missing += 1;
            continue;
        };
        let c = norm96(p);
        let t = s.theirs_linear;
        if c.iter().any(|v| *v <= GATE_LO || *v >= GATE_HI) {
            drop.cam_sat += 1;
            continue;
        }
        if t.iter().any(|v| *v <= GATE_LO || *v >= GATE_HI) {
            drop.ref_sat += 1;
            continue;
        }
        if x < cw / 2 {
            left.push((c, t));
        } else {
            right.push((c, t));
        }
    }
    (left, right, drop)
}

// ---------------------------------------------------------------------------
// 一张图的结果
// ---------------------------------------------------------------------------

struct ImageRow {
    stem: String,
    b_over_g: f64,
    r_over_g: f64,
    cam_mul: [f32; 3],
    n_train: usize,
    n_eval: usize,
    drop: DropCount,
    /// 左半解出的矩阵（主估计；**系数不可识别**，只作反面例证）。
    m_left: Mat3,
    /// 两半矩阵逐分量最大差（系数不可识别性的证据）。
    half_diff: f64,
    /// 两半矩阵的相对 Frobenius 差（同上）。
    coef_instab: f64,
    /// 训练集上的 rms（乐观，不用于判定）。
    train_rms: f64,
    /// 留出 rms：左半解 → 右半评估 / 右半解 → 左半评估。
    rms_opt: f64,
    rms_rev: f64,
    /// 系数尺读数（左半 / 右半）。
    dist_left: f64,
    dist_right: f64,
    gain_left: f64,
    gain_right: f64,
    shape_left: f64,
    shape_right: f64,
    /// 线性尺的四模型留出读数：rms 与 ΔE00 中位/P95，序为 [`MODELS`]。
    rms: [f64; 4],
    de: [f64; 4],
    de_p95: [f64; 4],
    /// 拟合出的整体增益与三通道增益（左半 / 右半），及按绿通道归一化的后者。
    k_fit: f64,
    w_fit: [f64; 3],
    w_norm: [f64; 3],
    /// 反向（右半拟合）的 `w`，用作同图噪声底。
    w_norm_rev: [f64; 3],
    /// 线性尺的三段份额：整体增益 / 对角增量 / 形状增量。
    share_global: f64,
    share_diag: f64,
    share_shape: f64,
    /// 色度尺读数：`M0` / `M0·diag(w)` / `M_opt`，以及反向自拟合（噪声底）。
    re: [RatioErr; 3],
    re_rev: RatioErr,
    /// 色度尺的三段份额：对角增量 / 形状增量（整体增益在色度尺上恒为 0）。
    re_share_diag: f64,
    re_share_shape: f64,
    /// 留出半幅的配对，供跨图迁移使用。
    right: Pairs,
    fixture_secs: f64,
    decode_secs: f64,
}

/// 解一张图。
fn analyze_image(stem: &str, m0: Mat3) -> Option<ImageRow> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        common::skip("诊断 G", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }

    // ---- 白平衡：相机记录的 cam_mul ----
    let cam_mul = match libraw::read_wb(&nef) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("跳过 {stem}：读白平衡失败（{e}）");
            return None;
        }
    };
    if cam_mul[1] <= 0.0 {
        eprintln!("跳过 {stem}：cam_mul 的 G 为 {}，无法表征白平衡", cam_mul[1]);
        return None;
    }
    let b_over_g = (cam_mul[2] / cam_mul[1]) as f64;
    let r_over_g = (cam_mul[0] / cam_mul[1]) as f64;

    // ---- 参考值走脚手架：与其它诊断共用同一份 ICC 转换与同一批评估样本 ----
    let t_fx = Instant::now();
    let Some(fx) = fitted_fixture_for(stem, false) else {
        eprintln!("跳过 {stem}：脚手架未能构建夹具（样本缺失或读取失败）");
        return None;
    };
    let fixture_secs = t_fx.elapsed().as_secs_f64();

    // ---- ① 原始解码（相机空间）：本诊断的自变量，必须自己解一次 ----
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let t_dec = Instant::now();
    let cam = match libraw::decode_with_output_for_test(&nef, &opts, OutputColor::Camera) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("跳过 {stem}：相机空间解码失败（{e}）");
            return None;
        }
    };
    let decode_secs = t_dec.elapsed().as_secs_f64();

    // ---- 与脚手架**同一套**裁切/朝向逻辑，把裁剪坐标映射回解码像素 ----
    let Some(model) = camera::read_model(&nef) else {
        eprintln!("跳过 {stem}：读不出机型");
        return None;
    };
    let Some(entry) = camera::lookup(&model) else {
        eprintln!("跳过 {stem}：机型 {model} 不在边距表中");
        return None;
    };
    let margins = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let Some((cw, ch)) = margins.effective(cam.width, cam.height) else {
        eprintln!("跳过 {stem}：边距与解码尺寸 {}×{} 不自洽", cam.width, cam.height);
        return None;
    };
    if (cw, ch) != fx.size {
        eprintln!(
            "跳过 {stem}：我们的裁切尺寸 {cw}×{ch} 与脚手架用过的 {}×{} 不一致，坐标无法对齐",
            fx.size.0, fx.size.1
        );
        return None;
    }

    let (left, right, drop) = collect_pairs(&fx, &cam, margins, cw);
    if left.len() < MIN_PAIRS || right.len() < MIN_PAIRS {
        eprintln!(
            "跳过 {stem}：可用配对不足（左半 {}、右半 {}，门槛 {MIN_PAIRS}；\
             丢弃——越界 {}、相机侧饱和 {}、参考侧饱和 {}）",
            left.len(),
            right.len(),
            drop.missing,
            drop.cam_sat,
            drop.ref_sat
        );
        return None;
    }

    let Some(dl) = color::derive_matrix(&left) else {
        eprintln!("跳过 {stem}：左半样本退化，解不出矩阵");
        return None;
    };
    let Some(dr) = color::derive_matrix(&right) else {
        eprintln!("跳过 {stem}：右半样本退化，解不出矩阵");
        return None;
    };
    let m_left = dl.matrix;
    let m_right = dr.matrix;

    let (dist_left, gain_left, shape_left) = matrix_metrics(m_left, m0);
    let (dist_right, gain_right, shape_right) = matrix_metrics(m_right, m0);
    let mut diff_sum = 0f64;
    for i in 0..3 {
        for j in 0..3 {
            diff_sum += (m_left[i][j] as f64 - m_right[i][j] as f64).powi(2);
        }
    }
    let coef_instab = diff_sum.sqrt() / frob(m_left).max(1e-12);

    // ---- 线性尺：四模型阶梯（都在左半拟合，右半评估）----
    let k_fit = best_gain(m0, &left);
    let Some(w_fit) = fit_diag_gain(m0, &left) else {
        eprintln!("跳过 {stem}：对角增益的正规方程退化");
        return None;
    };
    let Some(w_rev) = fit_diag_gain(m0, &right) else {
        eprintln!("跳过 {stem}：反向对角增益退化");
        return None;
    };
    if !w_fit[1].is_finite() || !w_rev[1].is_finite() || w_fit[1].abs() < 1e-9 || w_rev[1].abs() < 1e-9 {
        eprintln!("跳过 {stem}：绿通道增益接近 0，w 无法归一化");
        return None;
    }
    let s0 = score(&right, |c| mat3::mul_vec(m0, c));
    let sk = score(&right, |c| mat3::mul_vec(scaled(m0, k_fit as f32), c));
    let sd = score(&right, |c| apply_diag(m0, w_fit, c));
    let so = score(&right, |c| mat3::mul_vec(m_left, c));

    let total = s0.rms - so.rms;
    let share_of = |x: f64| if total > 1e-12 { x / total } else { f64::NAN };
    let share_global = share_of(s0.rms - sk.rms);
    let share_diag = share_of(sk.rms - sd.rms);
    let share_shape = share_of(sd.rms - so.rms);

    // ---- 色度尺：尺度不变，整体增益在这一把尺子上恒等 ----
    let r0 = ratio_error(&right, |c| mat3::mul_vec(m0, c));
    let rd = ratio_error(&right, |c| apply_diag(m0, w_fit, c));
    let ro = ratio_error(&right, |c| mat3::mul_vec(m_left, c));
    let rev = ratio_error(&left, |c| mat3::mul_vec(m_right, c));
    let re_total = r0.median - ro.median;
    let re_share_of = |x: f64| if re_total > 1e-12 { x / re_total } else { f64::NAN };
    let re_share_diag = re_share_of(r0.median - rd.median);
    let re_share_shape = re_share_of(rd.median - ro.median);

    let w_norm = [w_fit[0] / w_fit[1], 1.0, w_fit[2] / w_fit[1]];
    let w_norm_rev = [w_rev[0] / w_rev[1], 1.0, w_rev[2] / w_rev[1]];

    Some(ImageRow {
        stem: stem.to_string(),
        b_over_g,
        r_over_g,
        cam_mul,
        n_train: left.len(),
        n_eval: right.len(),
        drop,
        m_left,
        half_diff: max_abs_diff(m_left, m_right),
        coef_instab,
        train_rms: dl.rms,
        rms_opt: so.rms,
        rms_rev: score(&left, |c| mat3::mul_vec(m_right, c)).rms,
        dist_left,
        dist_right,
        gain_left,
        gain_right,
        shape_left,
        shape_right,
        rms: [s0.rms, sk.rms, sd.rms, so.rms],
        de: [s0.median, sk.median, sd.median, so.median],
        de_p95: [s0.p95, sk.p95, sd.p95, so.p95],
        k_fit,
        w_fit,
        w_norm,
        w_norm_rev,
        share_global,
        share_diag,
        share_shape,
        re: [r0, rd, ro],
        re_rev: rev,
        re_share_diag,
        re_share_shape,
        right,
        fixture_secs,
        decode_secs,
    })
}

// ---------------------------------------------------------------------------
// 分组与统计
// ---------------------------------------------------------------------------

/// 按 `B/G` 把图分档（`rows` 必须已按 `B/G` 升序）。
fn wb_groups(rows: &[ImageRow]) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        match out.last_mut() {
            Some(g) if r.b_over_g - rows[g[0]].b_over_g <= WB_GROUP_TOL => g.push(i),
            _ => out.push(vec![i]),
        }
    }
    out
}

/// 各档均值、档内合并标准差、档间极差。
///
/// `rev` 为 `Some` 时，同图的反向估计也进档内统计——那样「同图两半之差」也算进噪声底，
/// 噪声底更保守。
fn group_stats(rows: &[ImageRow], vals: &[f64], rev: Option<&[f64]>) -> (Vec<f64>, f64, f64) {
    let groups = wb_groups(rows);
    let mut means = Vec::new();
    let (mut ss, mut df) = (0f64, 0usize);
    for g in &groups {
        let mut xs: Vec<f64> = Vec::new();
        for &i in g {
            xs.push(vals[i]);
            if let Some(rv) = rev {
                xs.push(rv[i]);
            }
        }
        let m = xs.iter().sum::<f64>() / xs.len() as f64;
        for v in &xs {
            ss += (v - m).powi(2);
        }
        df += xs.len() - 1;
        means.push(m);
    }
    let sd = if df > 0 { (ss / df as f64).sqrt() } else { f64::NAN };
    let hi = means.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let lo = means.iter().cloned().fold(f64::INFINITY, f64::min);
    (means, sd, hi - lo)
}

fn fmt_means(means: &[f64]) -> String {
    means.iter().map(|m| format!("{m:.5}")).collect::<Vec<_>>().join(" / ")
}

/// 一行「标量 → 相关系数 + 档间/档内落差」。
fn print_scalar_line(rows: &[ImageRow], xs: &[f64], name: &str, vals: &[f64], rev: Option<&[f64]>) {
    let r = pearson(xs, vals);
    let (means, sd, range) = group_stats(rows, vals, rev);
    let snr = if sd > 0.0 { range / sd } else { f64::NAN };
    eprintln!(
        "    {:<26} r {:+6.3}  档均值 {:<34} 档间极差 {:.5}  档内 SD {:.5}  a/b {:>5.2}  {}",
        name,
        r,
        fmt_means(&means),
        range,
        sd,
        snr,
        snr_verdict(snr)
    );
}

// ---------------------------------------------------------------------------
// 打印
// ---------------------------------------------------------------------------

/// [0] 白平衡普查：`simple/` 下**全部** NEF（只 unpack，不解码）。
fn print_census() -> usize {
    let dir = samples_dir();
    let Ok(entries) = std::fs::read_dir(&dir) else {
        eprintln!("[0] 白平衡普查：读不到样本目录 {}，跳过", dir.display());
        return 0;
    };
    let mut items: Vec<(String, f64, f64)> = Vec::new();
    for e in entries.flatten() {
        let p = e.path();
        let is_nef = p.extension().map(|x| x.eq_ignore_ascii_case("nef")).unwrap_or(false);
        if !is_nef {
            continue;
        }
        let stem = p.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
        match libraw::read_wb(&p) {
            Ok(w) if w[1] > 0.0 => items.push((stem, (w[0] / w[1]) as f64, (w[2] / w[1]) as f64)),
            Ok(_) => eprintln!("  {stem}：cam_mul 的 G 非正，跳过"),
            Err(e) => eprintln!("  {stem}：读白平衡失败（{e}）"),
        }
    }
    if items.is_empty() {
        eprintln!("[0] 白平衡普查：目录里没有可用 NEF");
        return 0;
    }
    items.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..items.len() {
        match groups.last_mut() {
            Some(g) if items[i].2 - items[g[0]].2 <= WB_GROUP_TOL => g.push(i),
            _ => groups.push(vec![i]),
        }
    }

    eprintln!(
        "\n[0] 白平衡普查（simple/ 下全部 NEF，只 unpack 不解马赛克）：{} 个 NEF，按 B/G 分档（容差 {WB_GROUP_TOL}）",
        items.len()
    );
    for (gi, g) in groups.iter().enumerate() {
        let mut names: Vec<&str> = g.iter().map(|&i| items[i].0.as_str()).collect();
        names.sort_unstable();
        let bgs: Vec<f64> = g.iter().map(|&i| items[i].2).collect();
        let lo = bgs.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = bgs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        eprintln!(
            "  档 {}  B/G {:.4}（档内极差 {:.4}，占 {:.3}%）  R/G {:.4}  n = {}",
            gi + 1,
            items[g[0]].2,
            hi - lo,
            if items[g[0]].2 > 0.0 { 100.0 * (hi - lo) / items[g[0]].2 } else { 0.0 },
            items[g[0]].1,
            g.len()
        );
        eprintln!("        {}", names.join(" "));
    }
    eprintln!("  ⇒ 真正不同的白平衡只有 {} 档。", groups.len());
    eprintln!("     注：普查覆盖全部 NEF，只有同时有 .TIF 参考导出的图才进入求解（见下面的图列表）。");
    groups.len()
}

/// [1] 逐图细节。
fn print_detail(rows: &[ImageRow]) {
    eprintln!("\n[1] 逐图细节（左半拟合 → 右半评估；两半各解一次用于估计同图噪声底）");
    for r in rows {
        eprintln!(
            "\n  ── {}  B/G {:.4}  R/G {:.4}  cam_mul (R {:.4}, G {:.4}, B {:.4}) ──",
            r.stem, r.b_over_g, r.r_over_g, r.cam_mul[0], r.cam_mul[1], r.cam_mul[2]
        );
        eprintln!(
            "     样本：训练 {} / 评估 {}；丢弃——越界 {}、相机侧饱和 {}、参考侧饱和 {}",
            r.n_train, r.n_eval, r.drop.missing, r.drop.cam_sat, r.drop.ref_sat
        );
        let ml = fmt_rows(r.m_left);
        eprintln!("     逐图最优 3×3（**系数不可识别，仅存档**）M(左半) = {}", ml[0]);
        eprintln!("                                                 {}", ml[1]);
        eprintln!("                                                 {}", ml[2]);
        eprintln!(
            "       两半逐分量最大差 {:.5}；两半相对 Frobenius 差 {:.4}；训练 rms {:.6}（乐观）",
            r.half_diff, r.coef_instab, r.train_rms
        );
        eprintln!(
            "       留出 rms：左半解→右半 {:.6}、右半解→左半 {:.6}（同图噪声底）",
            r.rms_opt, r.rms_rev
        );
        eprintln!(
            "       系数尺（**不能**读成「矩阵差这么多」，见模块文档）：‖M−M0‖F/‖M0‖F {:.4}（反向 {:.4}）；矩阵空间增益 {:.4}（反向 {:.4}）；去增益形状 {:.4}（反向 {:.4}）",
            r.dist_left, r.dist_right, r.gain_left, r.gain_right, r.shape_left, r.shape_right
        );
        eprintln!(
            "     白平衡型修正 w（左半拟合，已按绿通道归一化）：R {:.5}  G 1.00000  B {:.5}；反向 {:.5} / {:.5}",
            r.w_norm[0], r.w_norm[2], r.w_norm_rev[0], r.w_norm_rev[2]
        );
        eprintln!(
            "     线性尺四模型留出 rms（{}）：{}",
            MODELS.join(" → "),
            r.rms.iter().map(|v| format!("{v:.5}")).collect::<Vec<_>>().join(" → ")
        );
        eprintln!(
            "       份额：整体增益 {:.1}% · 对角增量 {:.1}% · 形状增量 {:.1}%（拟合 k = {:.4}）",
            100.0 * r.share_global,
            100.0 * r.share_diag,
            100.0 * r.share_shape,
            r.k_fit
        );
        eprintln!(
            "       留出 ΔE00 中位：{}；P95：{}",
            r.de.iter().map(|v| format!("{v:.3}")).collect::<Vec<_>>().join(" → "),
            r.de_p95.iter().map(|v| format!("{v:.3}")).collect::<Vec<_>>().join(" → ")
        );
        eprintln!(
            "     色度尺（尺度不变）|Δlog2(R/G,B/G)| 中位：M0 {:.5} → M0·diag(w) {:.5} → M_opt {:.5}；反向自拟合 {:.5}（噪声底）",
            r.re[0].median, r.re[1].median, r.re[2].median, r.re_rev.median
        );
        eprintln!(
            "       同一读数的 P95：M0 {:.5} → M0·diag(w) {:.5} → M_opt {:.5}；可用样本 {} 个",
            r.re[0].p95, r.re[1].p95, r.re[2].p95, r.re[0].n
        );
        eprintln!(
            "       带符号中位：M0 的 Δlog2R/G {:+.5}、Δlog2B/G {:+.5}；色度落差里对角增量 {:.1}% 、形状增量 {:.1}%",
            r.re[0].med_rg, r.re[0].med_bg, 100.0 * r.re_share_diag, 100.0 * r.re_share_shape
        );
        eprintln!(
            "     耗时：脚手架 {:.1} s（含 260 MB 参考 TIF 与一次 ProPhoto 解码）、① 解码 {:.1} s",
            r.fixture_secs, r.decode_secs
        );
    }
}

/// [2] 系数尺（不可用）与线性尺的表。
fn print_compact(rows: &[ImageRow]) {
    eprintln!("\n[2] 摆在一起的对照表（按 B/G 升序）");
    eprintln!("  (a) 系数尺 —— **不可作判据**，列在这里是为了说明「矩阵差多少」量不出来：");
    eprintln!(
        "  {:<10} {:>7} {:>11} {:>11} {:>10} {:>10}",
        "图", "B/G", "‖ΔM‖F/‖M0‖", "去增益形状", "两半差", "两半相对差"
    );
    for r in rows {
        eprintln!(
            "  {:<10} {:>7.4} {:>11.4} {:>11.4} {:>10.5} {:>10.4}",
            r.stem, r.b_over_g, r.dist_left, r.shape_left, r.half_diff, r.coef_instab
        );
    }
    let max_half = rows.iter().map(|r| r.half_diff).fold(0.0f64, f64::max);
    let min_dist = rows.iter().map(|r| r.dist_left).fold(f64::INFINITY, f64::min);
    eprintln!(
        "      同一个 3×3 的系数能差到 {max_half:.2} 而功能几乎不变；系数离 M0 有 {min_dist:.1} 倍远时功能只改善十几个百分点。"
    );
    eprintln!("\n  (b) 线性尺 —— 受整体尺度主导（参考带 NEUTRAL 色调曲线，我们没有），只作参考：");
    eprintln!(
        "  {:<10} {:>7}  {:<46} {:>9} {:>9} {:>9}",
        "图", "B/G", format!("留出 rms：{}", MODELS.join(" → ")), "增益份额", "对角份额", "形状份额"
    );
    for r in rows {
        eprintln!(
            "  {:<10} {:>7.4}  {:.5} → {:.5} → {:.5} → {:.5}  {:>8.1}% {:>8.1}% {:>8.1}%",
            r.stem,
            r.b_over_g,
            r.rms[0],
            r.rms[1],
            r.rms[2],
            r.rms[3],
            100.0 * r.share_global,
            100.0 * r.share_diag,
            100.0 * r.share_shape
        );
    }
}

/// [3] 色度尺（判据）。
fn print_chroma(rows: &[ImageRow]) {
    eprintln!("\n[3] 色度尺（**判据**）：尺度不变的 |Δlog2(R/G), Δlog2(B/G)| 中位数");
    eprintln!("  这一把尺子上 M0 与 k·M0 读数相同（比值对整体增益不变），故只比三级。");
    eprintln!("  注意：`M0·diag(w)` 与 `M_opt` 都是按**线性 rms** 拟合的，这里只是把它们的留出预测");
    eprintln!("  换一把尺子打分，因此「M_opt 不比对角模型好」说明的是**多出来的形状自由度没换来色度精度**，");
    eprintln!("  不等于「色度意义上最优的 3×3 做不到更好」。");
    eprintln!(
        "  {:<10} {:>7} {:>10} {:>10} {:>10} {:>10}   {:<22} {:>9} {:>9}",
        "图", "B/G", "M0", "M0·diag(w)", "M_opt", "自拟合反向", "带符号 Δlog2R/G,Δlog2B/G", "对角份额", "形状份额"
    );
    for r in rows {
        eprintln!(
            "  {:<10} {:>7.4} {:>10.5} {:>10.5} {:>10.5} {:>10.5}   {:+8.5},{:+.5}      {:>8.1}% {:>8.1}%",
            r.stem,
            r.b_over_g,
            r.re[0].median,
            r.re[1].median,
            r.re[2].median,
            r.re_rev.median,
            r.re[0].med_rg,
            r.re[0].med_bg,
            100.0 * r.re_share_diag,
            100.0 * r.re_share_shape
        );
    }

    let xs: Vec<f64> = rows.iter().map(|r| r.b_over_g).collect();
    eprintln!("\n  色度尺标量的档间/档内落差（a/b = 档间极差 ÷ 档内合并标准差）：");
    let m0: Vec<f64> = rows.iter().map(|r| r.re[0].median).collect();
    let md: Vec<f64> = rows.iter().map(|r| r.re[1].median).collect();
    let mo: Vec<f64> = rows.iter().map(|r| r.re[2].median).collect();
    let mrev: Vec<f64> = rows.iter().map(|r| r.re_rev.median).collect();
    print_scalar_line(rows, &xs, "M0 的色度误差", &m0, None);
    print_scalar_line(rows, &xs, "M0·diag(w) 的色度误差", &md, None);
    print_scalar_line(rows, &xs, "M_opt 的色度误差", &mo, Some(&mrev));
    let sd: Vec<f64> = rows.iter().map(|r| r.re_share_diag).collect();
    let sh: Vec<f64> = rows.iter().map(|r| r.re_share_shape).collect();
    print_scalar_line(rows, &xs, "对角份额（色度）", &sd, None);
    print_scalar_line(rows, &xs, "形状份额（色度）", &sh, None);
}

/// [4] 白平衡型修正 `w`（3 参数、良态、直接可读）。
fn print_w(rows: &[ImageRow]) {
    eprintln!("\n[4] 白平衡型修正 w = M0·diag(w) 里那个 w（已按绿通道归一化）");
    eprintln!("  w_R 是「参考相对我们还需要多少红通道增益」，w_B 同理。1.0 = 我们的解码已经对上。");
    eprintln!(
        "  {:<10} {:>7} {:>10} {:>10} {:>10}   反向的 w_R / w_B",
        "图", "B/G", "w_R", "w_G", "w_B"
    );
    for r in rows {
        eprintln!(
            "  {:<10} {:>7.4} {:>10.5} {:>10.5} {:>10.5}   {:.5} / {:.5}",
            r.stem, r.b_over_g, r.w_norm[0], r.w_norm[1], r.w_norm[2], r.w_norm_rev[0], r.w_norm_rev[2]
        );
    }
    let xs: Vec<f64> = rows.iter().map(|r| r.b_over_g).collect();
    let wr: Vec<f64> = rows.iter().map(|r| r.w_norm[0]).collect();
    let wb: Vec<f64> = rows.iter().map(|r| r.w_norm[2]).collect();
    let wr_rev: Vec<f64> = rows.iter().map(|r| r.w_norm_rev[0]).collect();
    let wb_rev: Vec<f64> = rows.iter().map(|r| r.w_norm_rev[2]).collect();
    print_scalar_line(rows, &xs, "w_R", &wr, Some(&wr_rev));
    print_scalar_line(rows, &xs, "w_B", &wb, Some(&wb_rev));
}

/// [5] 跨图迁移：把 A 图解出的 `w_A` 用到 B 图的留出半幅上（**配 M0**，不是配 A 自己的矩阵）。
fn print_transfer(rows: &[ImageRow], m0: Mat3) {
    let groups = wb_groups(rows);
    let mut gid = vec![0usize; rows.len()];
    for (gi, g) in groups.iter().enumerate() {
        for &i in g {
            gid[i] = gi;
        }
    }

    eprintln!("\n[5] 跨图迁移：把 A 图左半解出的 w_A 用到 B 图右半，读尺度不变色度误差");
    eprintln!("  若「所需修正」随拍摄白平衡变化，跨档迁移应明显差于同档迁移。");
    let mut same: Vec<f64> = Vec::new();
    let mut cross: Vec<f64> = Vec::new();
    let mut selfv: Vec<f64> = Vec::new();
    for (b, rb) in rows.iter().enumerate() {
        selfv.push(rb.re[1].median);
        for (a, ra) in rows.iter().enumerate() {
            if a == b {
                continue;
            }
            let e = ratio_error(&rb.right, |c| apply_diag(m0, ra.w_fit, c)).median;
            if gid[a] == gid[b] {
                same.push(e);
            } else {
                cross.push(e);
            }
        }
    }
    let (sm, sa) = med_mean(&same);
    let (cm, ca) = med_mean(&cross);
    let (selfm, selfa) = med_mean(&selfv);
    eprintln!(
        "  同档迁移（A、B 同一 B/G 档）：{:>3} 对   中位 {:.5}   均值 {:.5}",
        same.len(),
        sm,
        sa
    );
    eprintln!(
        "  跨档迁移（A、B 不同 B/G 档）：{:>3} 对   中位 {:.5}   均值 {:.5}",
        cross.len(),
        cm,
        ca
    );
    eprintln!(
        "  自用（A = B，左半拟合 → 右半评估）：{:>3} 张   中位 {:.5}   均值 {:.5}",
        selfv.len(),
        selfm,
        selfa
    );
    let ratio = if sm > 0.0 { cm / sm } else { f64::NAN };
    eprintln!(
        "  ⇒ 跨档 ÷ 同档 = {ratio:.3}。（比值明显大于 1 才说明所需修正与拍摄白平衡有关；\
         接近 1 说明「换成哪张图算出来的修正」比「那张图什么白平衡」更重要。）"
    );

    eprintln!("\n  逐图对照（每一行是一张 B 图）：");
    eprintln!(
        "  {:<10} {:>7}  {:>12} {:>14} {:>14}",
        "B 图", "B/G", "自用", "同档 A 迁移中位", "跨档 A 迁移中位"
    );
    for (b, rb) in rows.iter().enumerate() {
        let mut s: Vec<f64> = Vec::new();
        let mut c: Vec<f64> = Vec::new();
        for (a, ra) in rows.iter().enumerate() {
            if a == b {
                continue;
            }
            let e = ratio_error(&rb.right, |cc| apply_diag(m0, ra.w_fit, cc)).median;
            if gid[a] == gid[b] {
                s.push(e);
            } else {
                c.push(e);
            }
        }
        eprintln!(
            "  {:<10} {:>7.4}  {:>12.5} {:>14} {:>14}",
            rb.stem,
            rb.b_over_g,
            rb.re[1].median,
            if s.is_empty() { "—".to_string() } else { format!("{:.5}", med_mean(&s).0) },
            if c.is_empty() { "—".to_string() } else { format!("{:.5}", med_mean(&c).0) }
        );
    }
}

/// [6] 判定。
fn print_verdict(rows: &[ImageRow], census_groups: usize) {
    let xs: Vec<f64> = rows.iter().map(|r| r.b_over_g).collect();
    let groups = wb_groups(rows);
    let mo: Vec<f64> = rows.iter().map(|r| r.re[2].median).collect();
    let mrev: Vec<f64> = rows.iter().map(|r| r.re_rev.median).collect();
    let (means, sd, range) = group_stats(rows, &mo, Some(&mrev));
    let snr = if sd > 0.0 { range / sd } else { f64::NAN };
    let r_opt = pearson(&xs, &mo);
    let monotone = means.windows(2).all(|w| w[1] > w[0]);

    let wb: Vec<f64> = rows.iter().map(|r| r.w_norm[2]).collect();
    let wb_rev: Vec<f64> = rows.iter().map(|r| r.w_norm_rev[2]).collect();
    let (w_means, w_sd, w_range) = group_stats(rows, &wb, Some(&wb_rev));
    let w_snr = if w_sd > 0.0 { w_range / w_sd } else { f64::NAN };

    let diag_share = mean_of(&rows.iter().map(|r| r.re_share_diag).collect::<Vec<_>>());
    let shape_share = mean_of(&rows.iter().map(|r| r.re_share_shape).collect::<Vec<_>>());
    let coef_ratio = median_of(&mut rows.iter().map(|r| r.dist_left).collect::<Vec<_>>());
    let coef_instab = median_of(&mut rows.iter().map(|r| r.coef_instab).collect::<Vec<_>>());

    eprintln!("\n[6] 判定");
    eprintln!("  样本：{} 张图参与求解；白平衡 {} 档（普查共 {} 档 NEF）。", rows.len(), groups.len(), census_groups.max(groups.len()));
    eprintln!(
        "  主标量（色度尺，M_opt 的留出色度误差）：均值 {:.5}；档均值（按 B/G 升序）{}",
        mean_of(&mo),
        fmt_means(&means)
    );
    eprintln!(
        "    r(B/G) = {r_opt:+.3}；档间极差 {range:.5} ÷ 档内合并 SD {sd:.5} = a/b {snr:.2} → {}",
        snr_verdict(snr)
    );
    eprintln!(
        "    档均值序列单调递增：{}（但这**不是**趋势证据，见第 1 条）",
        if monotone { "是" } else { "否" }
    );
    eprintln!(
        "  w_B（白平衡型修正）：档均值 {}；档间极差 {w_range:.5} ÷ 档内 SD {w_sd:.5} = a/b {w_snr:.2} → {}",
        fmt_means(&w_means),
        snr_verdict(w_snr)
    );
    eprintln!();
    eprintln!("  四条必须分开说（本项目已因「结论说得比证据满」返工两次）：");
    eprintln!(
        "  1) 单调趋势：本批白平衡只有 {} 档，**本诊断不判定单调性**——三个点恰好单调的",
        groups.len()
    );
    eprintln!("     先验概率是 1/3。任何「随白平衡越极端而越差」的说法都超出本批证据。");
    eprintln!("  2) 尺度：「矩阵差多少」在本数据上**量不出来**。系数尺显示 {coef_ratio:.1} 倍的偏离，");
    eprintln!("     但同一张图两半的系数相对差就有 {coef_instab:.2}，且那些非物理系数（2.3 / −3.4 / 2.8 一类）");
    eprintln!("     功能上只比 M0 好十几个百分点。⇒ 任何基于系数逐分量比较的结论都不可信；");
    eprintln!("     本诊断的判读全部走色度尺与 3 参数的 w。");
    eprintln!(
        "  3) 白平衡型修正是否随 B/G 变化：w_B 的读数「{}」（a/b {:.2}）；",
        snr_verdict(w_snr),
        w_snr
    );
    eprintln!(
        "     色度尺上 M_opt 的读数「{}」（a/b {:.2}）。判读依据是**档间落差是否超出档内离散**，",
        snr_verdict(snr),
        snr
    );
    eprintln!("     不是 r 的大小——r 在样本少、且含同档重复图时会被重复点主导。");
    eprintln!(
        "  4) 偏离的形态：色度落差里，白平衡型对角修正平均吃掉 {:.1}%，余下 {:.1}% 要靠完整 3×3 形状。",
        100.0 * diag_share,
        100.0 * shape_share
    );
    eprintln!("     （份额可以超过 100% 甚至为负：那说明 9 参数的 3×3 在留出集上并不比 3 参数的对角");
    eprintln!("     模型更准，即多出来的形状自由度没有换来色度精度。见 [3] 的拟合口径说明。）");
    if diag_share > 0.8 {
        eprintln!("     ⇒ 参考端与我们的色度差基本是「三通道增益」形态，不是矩阵形状。");
    } else if diag_share < 0.2 {
        eprintln!("     ⇒ 色度差基本吃不掉：形态是**完整 3×3**，与白平衡型（对角）修正不同一类。");
    } else {
        eprintln!("     ⇒ 对角与形状都占相当份额，不能只归给其中一类（形状份额还是**上界**，见下）。");
    }
    eprintln!();
    eprintln!("  限制（必须随结论一起转述）：");
    eprintln!("  · 色度尺只对**整体**增益免疫，对**逐通道**的曲线差不免疫；参考端 NEUTRAL 的逐通道");
    eprintln!("    色调曲线会被算进「形状份额」里，故形状份额是上界、对角份额是下界。");
    eprintln!("  · M_opt 逐图独立，吸收了该图的一切线性差异（曝光、黑电平、参考端残差），");
    eprintln!("    因此它**不是物理相机矩阵**，只能用来看「线性上界」。");
    eprintln!("  · 参考 TIF 若在 NX Studio 里改过白平衡，该图测的就是「两份白平衡之差」而非「两种渲染之差」；");
    eprintln!("    [4] 的 w 表正是查这件事的地方（w 明显偏离 1 的图要单独存疑）。");
}

// ---------------------------------------------------------------------------
// 入口
// ---------------------------------------------------------------------------

fn selected_stems() -> Vec<String> {
    match std::env::var("NRV_DIAG_STEM") {
        Ok(v) if !v.trim().is_empty() => v
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect(),
        _ => DEFAULT_STEMS.iter().map(|s| s.to_string()).collect(),
    }
}

fn census_enabled() -> bool {
    !matches!(std::env::var("NRV_DIAG_CENSUS").as_deref(), Ok("0") | Ok("no") | Ok("false"))
}

#[test]
fn optimal_matrix_vs_white_balance() {
    eprintln!("=== 诊断 G：逐图最优变换是否随拍摄白平衡系统变化 ===");
    let Some(cm) = camera::color_matrix("Nikon Z 8") else {
        eprintln!("跳过：机型表里没有 Z8 的色彩矩阵，无从比较");
        return;
    };
    let m0 = cm.matrix;
    eprintln!("固定矩阵 M0（{}，来源：{}）", cm.output_space, cm.source);
    let rows0 = fmt_rows(m0);
    eprintln!("  {}", rows0[0]);
    eprintln!("  {}", rows0[1]);
    eprintln!("  {}", rows0[2]);

    let census_groups = if census_enabled() {
        print_census()
    } else {
        eprintln!("\n[0] 白平衡普查：按 NRV_DIAG_CENSUS 跳过");
        0
    };

    let stems = selected_stems();
    eprintln!("\n将要参与求解的图（{} 张）：{}", stems.len(), stems.join(" "));
    let mut rows: Vec<ImageRow> = Vec::new();
    for stem in &stems {
        if let Some(r) = analyze_image(stem, m0) {
            rows.push(r);
        }
    }
    if rows.is_empty() {
        eprintln!("跳过：没有一张图能同时取到 NEF、参考 TIF 与可用配对");
        return;
    }
    rows.sort_by(|a, b| a.b_over_g.partial_cmp(&b.b_over_g).unwrap_or(std::cmp::Ordering::Equal));

    print_detail(&rows);
    print_compact(&rows);
    print_chroma(&rows);
    print_w(&rows);
    print_transfer(&rows, m0);
    print_verdict(&rows, census_groups);

    eprintln!("\n口径提醒：线性尺上的 ΔE00 是在**没有色调曲线**的裸 3×3 上算的，");
    eprintln!("绝对值远大于 5.5/5.6 报告的 1.698 / 3.543 / 4.767，只能在同一批样本之间横向比较；");
    eprintln!("判定请以 [3] 的尺度不变色度尺与 [4] 的 3 参数 w 为准。");
}

/// 仪器自检：对角增益模型必须能解回已知的 `w`，且此时它与无约束 3×3 给出同一个矩阵。
///
/// 这条同时说明「`M0·diag(w)` 是无约束 3×3 的真子集」——[3]/[6] 里「对角份额」的读法依赖这一点。
#[test]
fn instrument_recovers_a_known_diagonal_gain() {
    let m0 = camera::color_matrix("Nikon Z 8").expect("Z8 应有固定矩阵").matrix;
    let w_true = [1.07f64, 1.0, 0.93];
    let mut pairs: Pairs = Vec::new();
    for i in 0..60 {
        for j in 0..60 {
            let c = [
                0.05 + i as f32 * 0.014,
                0.05 + j as f32 * 0.014,
                0.05 + ((i * 7 + j * 13) % 60) as f32 * 0.014,
            ];
            let cw = [
                c[0] * w_true[0] as f32,
                c[1] * w_true[1] as f32,
                c[2] * w_true[2] as f32,
            ];
            pairs.push((c, mat3::mul_vec(m0, cw)));
        }
    }

    let w = fit_diag_gain(m0, &pairs).expect("应能求解");
    for (k, (got, want)) in w.iter().zip(w_true.iter()).enumerate() {
        assert!((got - want).abs() < 1e-4, "w[{k}] = {got}，期望 {want}（全部 {w:?}）");
    }

    let full = color::derive_matrix(&pairs).expect("应能求解");
    let expect = mat3::mul(
        m0,
        [[w_true[0] as f32, 0.0, 0.0], [0.0, w_true[1] as f32, 0.0], [0.0, 0.0, w_true[2] as f32]],
    );
    for (i, row) in expect.iter().enumerate() {
        for (j, want) in row.iter().enumerate() {
            let got = full.matrix[i][j];
            assert!(
                (got - want).abs() < 1e-4,
                "[{i}][{j}] 无约束 3×3 应等于 M0·diag(w)：{got} vs {want}"
            );
        }
    }
    assert!(full.rms < 1e-5, "残差应可忽略，得到 {}", full.rms);
    assert!(best_gain(m0, &pairs) > 0.0, "标量增益应为正");

    // 色度尺对整体增益不敏感：k·M0 与 M0 必须给出同一个色度读数
    let e0 = ratio_error(&pairs, |c| mat3::mul_vec(m0, c));
    let ek = ratio_error(&pairs, |c| mat3::mul_vec(scaled(m0, 7.5), c));
    assert!(
        (e0.median - ek.median).abs() < 1e-9,
        "色度尺应对整体增益免疫：{:.9} vs {:.9}",
        e0.median,
        ek.median
    );
    assert!(e0.n > 0 && e0.median.is_finite(), "色度尺应有可用样本，得到 {e0:?}");
}
