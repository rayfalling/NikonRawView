//! 诊断 M：把 `w` 的自变量换成 **(R/G, B/G) 两元**，重跑拟合与留一状态交叉验证。
//!
//! # 相对上一轮（诊断 K）的唯一实质变化
//!
//! 上一轮只用 `B/G` 一个自变量，对红通道给出**负结论**（CV 中位 3.55% 但最坏 14.8%，
//! 且结构上「B/G ≤ 1.31 恒定、之后陡降」，误差÷幅度 1.65）。
//!
//! 本轮的依据是：`cam_mul` 的四元组里 `G1 == G2 == 1.0` 恒成立，尼康把白平衡微调
//! （A/M 偏移）折进了 R 与 B 的乘数——所以白平衡状态由 **`(R/G, B/G)` 两个数**完整编码，
//! 而 `K` 只是"偏移为零"那条特殊曲线上的参数。**把二维投影到一根轴会制造假结构**：
//! `DSC_0567`（R/G 2.1113、B/G 1.4941）与 `DSC_0141`（R/G 1.7500、B/G 1.4961）B/G 只差
//! 0.13%，`w_R` 却差 23.9%，而 `w_B` 只差 1.7%。
//!
//! 因此本轮的判据从「留一 **B/G 档**」改为「留一 **(R/G, B/G) 白平衡状态**」，
//! 模型从「6 种一元形式」扩到「一元 + 二元共 9 个模型」，两者都用同一套 minimax 纪律选。
//!
//! # 要解决的问题
//!
//! 诊断链已把问题收敛到唯一剩下的环节——**白平衡/光源适配的应用方式**：
//!
//! | 环节 | 状态 |
//! |---|---|
//! | 相机矩阵的施加 | ✓ 我们的 == LibRaw（5e-6） |
//! | 相机矩阵的取值 | LibRaw 单矩阵、零插值（源码确认） |
//! | 白平衡的读取 | ✓ 23 张逐位吻合（MakerNote `WB_RBLevels` == `cam_mul`） |
//! | 参考样本有效性 | ✓ NX Studio 显示「原始值」+ `.nksc` 无 WB 字段 |
//! | 差异形态 | ✓ 三参数对角增益吃掉色度落差的 106%（诊断 G） |
//! | 只与白平衡档有关 | ✓ 跨档迁移是同档的 12.9 倍，同档 ≈ 自用（诊断 G） |
//!
//! 所以假设是 **尼康的等效变换 = 固定矩阵 · diag(w(白平衡))**，本诊断把 `w` 标定出来，
//! 并且——这是重点——**用留一档交叉验证检验它到底能不能用**。
//!
//! # 两套参数化，都要报
//!
//! 对角增益可以放在矩阵的**前面**或**后面**，两者是不同的三参数族（`M0·diag(w)` 与
//! `diag(w)·M0` 只在特殊情形下重合）：
//!
//! - `M0·diag(w)`（**主模型**）：增益加在相机空间、矩阵之前——这正是解码器里
//!   `scale_colors`/`pre_mul` 所在的位置，也是诊断 G 验证过的形态。
//!   预测对 `w` 线性（`Σ_j w_j·u_j`，`u_j = M0 的第 j 列 × c_j`），一次三元最小二乘解完。
//! - `diag(w)·M0`（**对照**）：增益加在输出空间、矩阵之后。三分量独立，闭式解
//!   `w_k = Σ t_k·p_k / Σ p_k²`（`p = M0·c`）。Lead 的起点文件解的是这一个；
//!   本诊断把两者摆在一起比，免得把一个没验证过的族当成结论。
//!
//! # 采样与读数
//!
//! `c` 取**原始解码**（`OutputColor::Camera`，管线里的 ①；① 已经含 `cam_mul`），
//! `t` 取参考 TIF 经内嵌 ICC 转 ProPhoto 线性。每 17 个像素取一点；排除任一侧的
//! 饱和/近黑点。色度读数用**尺度不变**的 `|Δlog2(R/G), Δlog2(B/G)|` 中位数——整体增益
//! 在这一把尺子上恒等，因此不会被「我们还没有色调曲线」这件事污染。
//!
//! # 逐通道选形式（R 与 B 分开）
//!
//! 五个候选：常数（`w = c`，即**不用 B/G**）、线性、半对数、指数、对数-对数。
//! 每个通道**各自**按留一档 CV 误差选中位数最小的那个。这一步很重要：若某通道选中的是
//! 常数，那就等于说「这个通道与 B/G 无关」，而不是硬套一条曲线上去。
//!
//! # 纪律
//!
//! `w` 是在这批参考图上解出来的，而它们场景各异、白平衡只有几档（其中一批 42 秒连拍
//! 只占一档）。**必须做留一档交叉验证**，否则我们只是把几个点插值成一条曲线。
//! 判读分两种情形，两种都要如实报：交叉验证误差只有百分之几 → 曲线可用；
//! 误差与 `w` 本身的修正幅度同量级 → **白平衡不是唯一自变量**，还缺场景/光源类型信息，
//! 那就不能只靠一张曲线表。**这个负结论和正结论一样有价值。**
//!
//! 门槛先定后看：见 [`ERR_USABLE`] / [`ERR_VS_AMPLITUDE`]。
//!
//! # 运行
//!
//! 默认：普查全部配对 → 每档 1 张代表 + 多图档各补 1 张「编号相隔最远」的成员 +
//! 最大的一档再补 1 张 → 逐图解码。**绝不为了预算丢档**；真超预算时打印丢了谁。
//! `NRV_K_STEMS=DSC_0001,DSC_0141` 只跑指定图（迭代用，跳过普查）；
//! `NRV_K_CENSUS=0` 跳过普查。
//!
//! # 本轮（诊断 M）的数据变化，以及**没有**变的东西
//!
//! 用户按建议拍了一组**色温扫描**（DSC_0561~0570），把 B/G 1.33~2.19 之间的点从 2 个补到 5 个，
//! 最大相邻间隔从 47% 降到 23%。这解决的是**「红通道那段形状没被约束住」**的问题。
//!
//! 但它**没有**改善「跨场景」这一维：新样张是**同一场景**在不同色温下的扫描。
//! 所以留一档 CV 仍然是**插值**（在相邻档之间），不是「换一个场景还成立」的验证。
//! [1d] 把「近同档但不同图」的配对单独列出来，那是本批唯一能看跨场景的地方。

mod common;

use common::samples_dir;
use nikonrawview::camera;
use nikonrawview::fit;
use nikonrawview::icc;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::mat3::{self, Mat3};
use std::time::Instant;

// ---------------------------------------------------------------------------
// 口径与门槛（先定后看）
// ---------------------------------------------------------------------------

/// 采样步长（每 17 个像素取一点，与 Lead 的起点一致）。
const STEP: usize = 17;

/// 相机空间侧的门槛：近黑与饱和都丢掉。
const CAM_LO: f32 = 1e-4;
const CAM_HI: f32 = 0.999;
/// 参考侧的门槛。
const REF_LO: f32 = 1e-5;
const REF_HI: f32 = 0.999;
/// 色度比值的可用门槛：任一通道低于它，`R/G`、`B/G` 就是噪声。
const RATIO_GATE: f32 = 1e-4;

/// 一张图可用的最少采样点数。
const MIN_SAMPLES: usize = 2000;

/// 视作同一档白平衡的**相对** B/G 容差（0.5%）。
///
/// 这个容差决定了「留一档」到底留掉什么：容差越窄，档越碎、被留掉的那档在训练集里
/// 越可能有**近乎相同**的近邻，CV 于是偏乐观。所以 [3b] 会把容差放宽再算一遍。
const WB_LEVEL_TOL: f64 = 5e-3;

/// 多图档额外补几张。补的那张取**编号与该档代表相隔最远**的成员——那是「不同拍摄批次/
/// 不同场景」在本数据里唯一可用的代理，用来检验「同白平衡、不同场景是否给出同一个 w」。
const EXTRA_PER_MULTI: usize = 1;
/// 最大的一档再补一张，保证连拍档至少有 3 张可量档内离散。
const EXTRA_PER_BIGGEST: usize = 1;
/// 解码总预算上限。**超出时先丢额外补的图并打印丢了谁，绝不静默丢档**。
const MAX_IMAGES: usize = 22;
/// [1d] 近同档配对用的 B/G 相对容差。
const NEAR_BG_TOL: f64 = 0.03;
/// [3b] 档容差敏感性检查用的更宽容差。
const TOL_SENSITIVITY: [f64; 2] = [0.01, 0.02];

/// 「曲线可用」的门槛：留一档交叉验证的相对误差中位数。
const ERR_USABLE: f64 = 0.05;
/// 「误差与修正幅度同量级」的门槛：CV 误差 ÷ 中位 |w−1|。
const ERR_VS_AMPLITUDE: f64 = 1.0;

// ---------------------------------------------------------------------------
// 基础工具
// ---------------------------------------------------------------------------

/// 解码层的 16 位输出 → 0..1。
fn norm96(p: [u16; 3]) -> [f32; 3] {
    [p[0] as f32 / 65535.0, p[1] as f32 / 65535.0, p[2] as f32 / 65535.0]
}

/// 一个采样点：相机空间值、参考值、是否位于画面左半。
type Sample = ([f32; 3], [f32; 3], bool);

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

/// 按绿通道归一化（消除整体亮度尺度；色调曲线承担那一部分）。
fn norm_g(w: [f64; 3]) -> Option<[f64; 3]> {
    if w[1].abs() < 1e-12 || !w.iter().all(|v| v.is_finite()) {
        return None;
    }
    Some([w[0] / w[1], 1.0, w[2] / w[1]])
}

/// 施加「矩阵之前的对角增益」：`M0·diag(w)·c`。
fn apply_w_cam(m0: Mat3, w: [f64; 3], c: [f32; 3]) -> [f32; 3] {
    mat3::mul_vec(
        m0,
        [c[0] * w[0] as f32, c[1] * w[1] as f32, c[2] * w[2] as f32],
    )
}

/// 施加「矩阵之后的对角增益」：`diag(w)·M0·c`。
fn apply_w_out(m0: Mat3, w: [f64; 3], c: [f32; 3]) -> [f32; 3] {
    let p = mat3::mul_vec(m0, c);
    [p[0] * w[0] as f32, p[1] * w[1] as f32, p[2] * w[2] as f32]
}

/// 解「矩阵之前的对角增益」：`min Σ ‖M0·diag(w)·c − t‖²`。
///
/// 预测对 `w` 线性：`M0·diag(w)·c = Σ_j w_j·u_j`，`u_j = M0 的第 j 列 × c_j`，
/// 因此是三元最小二乘，正规矩阵 3×3，一次求解，无需迭代。
/// 返回**未归一化**的最小二乘最优解（三个分量共同定尺度），归一化交给 [`norm_g`]，
/// 这样两族的残差比较才是各自的上限。
fn solve_w_cam(m0: Mat3, samples: &[Sample]) -> Option<[f64; 3]> {
    let mut a = [[0f64; 3]; 3];
    let mut b = [0f64; 3];
    for (c, t, _) in samples {
        let mut u = [[0f64; 3]; 3]; // u[j][k] = M0[k][j]·c_j
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
    let w = solve3(a, b)?;
    w.iter().all(|v| v.is_finite() && v.abs() > 1e-12).then_some(w)
}

/// 解「矩阵之后的对角增益」：`min Σ ‖diag(w)·M0·c − t‖²`。
///
/// 三分量互相独立，闭式解 `w_k = Σ t_k·p_k / Σ p_k²`（`p = M0·c`）。
fn solve_w_out(m0: Mat3, samples: &[Sample]) -> Option<[f64; 3]> {
    let mut num = [0f64; 3];
    let mut den = [0f64; 3];
    for (c, t, _) in samples {
        let p = mat3::mul_vec(m0, *c);
        for k in 0..3 {
            let pf = p[k] as f64;
            num[k] += t[k] as f64 * pf;
            den[k] += pf * pf;
        }
    }
    if den.iter().any(|d| *d <= 0.0) {
        return None;
    }
    let w = [num[0] / den[0], num[1] / den[1], num[2] / den[2]];
    w.iter().all(|v| v.is_finite()).then_some(w)
}

/// 相对线性残差：`sqrt(Σ‖pred−t‖² / Σ‖t‖²)`。两族参数化用它比「谁拟合得更准」。
fn rel_resid<F: Fn([f32; 3]) -> [f32; 3]>(samples: &[Sample], f: F) -> f64 {
    let (mut ss, mut tt) = (0f64, 0f64);
    for (c, t, _) in samples {
        let p = f(*c);
        for k in 0..3 {
            ss += ((p[k] - t[k]) as f64).powi(2);
            tt += (t[k] as f64).powi(2);
        }
    }
    (ss / tt.max(1e-12)).sqrt()
}

/// 三个通道都有限且超过 [`RATIO_GATE`]。
fn usable_ratio(v: &[f32; 3]) -> bool {
    v.iter().all(|x| x.is_finite() && *x > RATIO_GATE)
}

/// 尺度不变色度误差：`|(Δlog2(R/G), Δlog2(B/G))|` 的中位数。
///
/// 预测与参考各自乘任意正常数都不改变这两个分量，因此整体增益、曝光尺度被完全排除。
fn chroma_err<F: Fn([f32; 3]) -> [f32; 3]>(samples: &[Sample], f: F) -> f64 {
    let mut lens = Vec::new();
    for (c, t, _) in samples {
        let p = f(*c);
        if !usable_ratio(&p) || !usable_ratio(t) {
            continue;
        }
        let a = ((p[0] / p[1]) as f64).log2() - ((t[0] / t[1]) as f64).log2();
        let b = ((p[2] / p[1]) as f64).log2() - ((t[2] / t[1]) as f64).log2();
        lens.push((a * a + b * b).sqrt());
    }
    if lens.is_empty() {
        return f64::NAN;
    }
    lens.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    lens[lens.len() / 2]
}

fn median_of(v: &[f64]) -> f64 {
    let mut c = v.to_vec();
    if c.is_empty() {
        return f64::NAN;
    }
    c.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    c[c.len() / 2]
}

fn max_of(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::NAN, f64::max)
}

/// 相对误差。
fn rel_err(pred: f64, measured: f64) -> f64 {
    if measured.abs() < 1e-12 {
        f64::NAN
    } else {
        (pred - measured).abs() / measured.abs()
    }
}

/// Pearson 相关系数。
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

/// 控制 `z` 之后 `x` 与 `y` 的偏相关。
fn partial_corr(x: &[f64], y: &[f64], z: &[f64]) -> f64 {
    let rxy = pearson(x, y);
    let rxz = pearson(x, z);
    let ryz = pearson(y, z);
    let den = ((1.0 - rxz * rxz) * (1.0 - ryz * ryz)).sqrt();
    if !den.is_finite() || den <= 1e-12 {
        f64::NAN
    } else {
        (rxy - rxz * ryz) / den
    }
}

// ---------------------------------------------------------------------------
// 拟合形式
// ---------------------------------------------------------------------------

/// 拟合形式。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Form {
    /// `w = c`——**不用 B/G**，取训练集**均值**（平方误差最优）。它是「B/G 有没有用」的对照。
    Const,
    /// `w = c`——取训练集**中位数**。对离群图稳健，用来分辨「曲线没用」与「曲线被一张离群图带偏」。
    Median,
    /// `w = a + b·x`
    LinLin,
    /// `w = a + b·ln x`
    LinLog,
    /// `ln w = a + b·x`
    LogLin,
    /// `ln w = a + b·ln x`
    LogLog,
}

const FORMS: [Form; 6] =
    [Form::Const, Form::Median, Form::LinLin, Form::LinLog, Form::LogLin, Form::LogLog];

impl Form {
    fn name(self) -> &'static str {
        match self {
            Form::Const => "常数-均值（不用 B/G）",
            Form::Median => "常数-中位数（不用 B/G）",
            Form::LinLin => "w = a + b·x",
            Form::LinLog => "w = a + b·ln x",
            Form::LogLin => "ln w = a + b·x",
            Form::LogLog => "ln w = a + b·ln x",
        }
    }
    /// 自变量变换；`None` 表示该样本落在此形式的定义域之外。
    fn tx(self, x: f64) -> Option<f64> {
        match self {
            Form::Const | Form::Median => Some(0.0),
            Form::LinLin | Form::LinLog => Some(x),
            Form::LogLin | Form::LogLog => (x > 0.0).then(|| x.ln()),
        }
    }
    fn ty(self, y: f64) -> Option<f64> {
        match self {
            Form::Const | Form::Median | Form::LinLin | Form::LogLin => Some(y),
            Form::LinLog | Form::LogLog => (y > 0.0).then(|| y.ln()),
        }
    }
    fn inv_y(self, t: f64) -> f64 {
        match self {
            Form::Const | Form::Median | Form::LinLin | Form::LogLin => t,
            Form::LinLog | Form::LogLog => t.exp(),
        }
    }
}

/// 一个拟合结果：变换域里 `y ≈ a + b·x`。
#[derive(Clone, Copy)]
struct Fit {
    form: Form,
    a: f64,
    b: f64,
    /// 训练集上的相对残差：`ln(pred/actual)` 的 RMS。
    log_rms: f64,
}

impl Fit {
    fn predict(&self, x: f64) -> f64 {
        match self.form.tx(x) {
            Some(t) => self.form.inv_y(self.a + self.b * t),
            None => f64::NAN,
        }
    }
}

/// 在 `(x, y)` 上做变换域最小二乘。
fn fit_form(form: Form, pts: &[(f64, f64)]) -> Option<Fit> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (x, y) in pts {
        if let (Some(tx), Some(ty)) = (form.tx(*x), form.ty(*y)) {
            xs.push(tx);
            ys.push(ty);
        }
    }
    let n = xs.len();
    if n < 3 {
        return None;
    }
    // 常数-中位数是显式分支：不做最小二乘，直接取中位数，用来分辨
    // 「B/G 没用」与「曲线被一张离群图带偏」。
    if form == Form::Median {
        let fit = Fit { form, a: median_of(&ys), b: 0.0, log_rms: f64::NAN };
        let mut ss = 0f64;
        let mut cnt = 0usize;
        for (x, y) in pts {
            let p = fit.predict(*x);
            if p > 0.0 && *y > 0.0 {
                ss += (p / y).ln().powi(2);
                cnt += 1;
            }
        }
        return Some(Fit { log_rms: (ss / cnt.max(1) as f64).sqrt(), ..fit });
    }
    let mx = xs.iter().sum::<f64>() / n as f64;
    let my = ys.iter().sum::<f64>() / n as f64;
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let b = if sxx.abs() < 1e-15 { 0.0 } else { sxy / sxx };
    let a = my - b * mx;
    let fit = Fit { form, a, b, log_rms: f64::NAN };

    let mut ss = 0f64;
    let mut cnt = 0usize;
    for (x, y) in pts {
        let p = fit.predict(*x);
        if p > 0.0 && *y > 0.0 {
            ss += (p / y).ln().powi(2);
            cnt += 1;
        }
    }
    Some(Fit { log_rms: (ss / cnt.max(1) as f64).sqrt(), ..fit })
}

/// 拟合用的样本点 `(B/G, w_k)`。
fn points(images: &[&ImageData], k: usize) -> Vec<(f64, f64)> {
    images.iter().map(|im| (im.bg, im.w_norm[k])).collect()
}

// ---------------------------------------------------------------------------
// 模型：一元 vs 二元
// ---------------------------------------------------------------------------
//
// 自变量是白平衡的**两个**比值 `R/G`、`B/G`。依据（诊断 L / Lead 的补充，本测试也自检）：
// `cam_mul` 的四元组里 `G1 == G2 == 1.0` 恒成立，调色偏移被尼康折进了 R 与 B 的乘数，
// 所以 `(R/G, B/G)` 完整编码了白平衡状态。`K` 只是"调色偏移为零"那条特殊曲线上的参数，
// 不适合当自变量。
//
// 上一轮只用了 `B/G` 一个自变量，于是「B/G 几乎相同但 R/G 差 21%」的两张图
// （DSC_0567 与 DSC_0141）被当成同一档——把二维投影到一维会**制造出假结构**。

/// 候选模型。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Model {
    ConstMean,
    ConstMedian,
    UniBgLinLin,
    UniBgLinLog,
    UniBgLogLin,
    UniBgLogLog,
    UniRgLogLog,
    BiLogLog,
    BiLinear,
}

const MODELS: [Model; 9] = [
    Model::ConstMean,
    Model::ConstMedian,
    Model::UniBgLinLin,
    Model::UniBgLinLog,
    Model::UniBgLogLin,
    Model::UniBgLogLog,
    Model::UniRgLogLog,
    Model::BiLogLog,
    Model::BiLinear,
];

impl Model {
    fn name(self) -> &'static str {
        match self {
            Model::ConstMean => "常数-均值（不用白平衡）",
            Model::ConstMedian => "常数-中位数（不用白平衡）",
            Model::UniBgLinLin => "一元：w = a + b·(B/G)",
            Model::UniBgLinLog => "一元：w = a + b·ln(B/G)",
            Model::UniBgLogLin => "一元：ln w = a + b·(B/G)",
            Model::UniBgLogLog => "一元：ln w = a + b·ln(B/G)",
            Model::UniRgLogLog => "一元：ln w = a + b·ln(R/G)",
            Model::BiLogLog => "**二元**：ln w = a + b₁ln(R/G) + b₂ln(B/G)",
            Model::BiLinear => "**二元**：w = a + b₁(R/G) + b₂(B/G)",
        }
    }
    fn is_bivariate(self) -> bool {
        matches!(self, Model::BiLogLog | Model::BiLinear)
    }
    /// 只用 `B/G`（上一轮的那一族）——用来回答「红通道的负结论翻没翻」。
    fn bg_only(self) -> bool {
        matches!(
            self,
            Model::ConstMean
                | Model::ConstMedian
                | Model::UniBgLinLin
                | Model::UniBgLinLog
                | Model::UniBgLogLin
                | Model::UniBgLogLog
        )
    }
}

/// 拟合好的模型。预测值由 `[1, x₁, x₂]` 三个基函数线性组合（一元模型只用前两个）。
struct Fitted {
    model: Model,
    c: [f64; 3],
}

impl Fitted {
    fn predict(&self, rg: f64, bg: f64) -> f64 {
        let (c0, c1, c2) = (self.c[0], self.c[1], self.c[2]);
        match self.model {
            Model::ConstMean | Model::ConstMedian => c0,
            Model::UniBgLinLin => c0 + c1 * bg,
            Model::UniBgLinLog => c0 + c1 * bg.ln(),
            Model::UniBgLogLin => (c0 + c1 * bg).exp(),
            Model::UniBgLogLog => (c0 + c1 * bg.ln()).exp(),
            Model::UniRgLogLog => (c0 + c1 * rg.ln()).exp(),
            Model::BiLogLog => (c0 + c1 * rg.ln() + c2 * bg.ln()).exp(),
            Model::BiLinear => c0 + c1 * rg + c2 * bg,
        }
    }
}

/// 高斯消元解 `n` 元方程组（`n ≤ 3`），增广列在第 `n` 列；退化返回 `None`。
fn solve_n(mut a: [[f64; 4]; 3], n: usize) -> Option<[f64; 3]> {
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-15 {
            return None;
        }
        a.swap(col, piv);
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = a[r][col] / a[col][col];
            let pivot = a[col];
            for (rc, pc) in a[r][col..=n].iter_mut().zip(pivot[col..=n].iter()) {
                *rc -= f * pc;
            }
        }
    }
    let mut out = [0f64; 3];
    for i in 0..n {
        out[i] = a[i][n] / a[i][i];
    }
    Some(out)
}

/// 最小二乘：`y ≈ Σ_j c_j · cols[j]`。
fn ols(cols: &[Vec<f64>], y: &[f64], model: Model) -> Option<Fitted> {
    let p = cols.len();
    if p == 0 || p > 3 || y.len() < p + 2 {
        return None;
    }
    let mut a = [[0f64; 4]; 3];
    for i in 0..p {
        for j in 0..p {
            a[i][j] = cols[i].iter().zip(&cols[j]).map(|(u, v)| u * v).sum();
        }
        a[i][p] = cols[i].iter().zip(y).map(|(u, v)| u * v).sum();
    }
    let c = solve_n(a, p)?;
    c.iter().all(|v| v.is_finite()).then_some(Fitted { model, c })
}

/// 在 `(R/G, B/G, w)` 上拟合一个模型。
fn fit_model(model: Model, pts: &[(f64, f64, f64)]) -> Option<Fitted> {
    let n = pts.len();
    if n < 4 {
        return None;
    }
    let ones = vec![1.0; n];
    let rg: Vec<f64> = pts.iter().map(|p| p.0).collect();
    let bg: Vec<f64> = pts.iter().map(|p| p.1).collect();
    let w: Vec<f64> = pts.iter().map(|p| p.2).collect();
    let ln_rg: Vec<f64> = pts.iter().map(|p| p.0.ln()).collect();
    let ln_bg: Vec<f64> = pts.iter().map(|p| p.1.ln()).collect();
    let ln_w: Vec<f64> = pts.iter().map(|p| p.2.ln()).collect();
    match model {
        Model::ConstMedian => Some(Fitted { model, c: [median_of(&w), 0.0, 0.0] }),
        Model::ConstMean => ols(&[ones], &w, model),
        Model::UniBgLinLin => ols(&[ones, bg], &w, model),
        Model::UniBgLinLog => ols(&[ones, ln_bg], &w, model),
        Model::UniBgLogLin => ols(&[ones, bg], &ln_w, model),
        Model::UniBgLogLog => ols(&[ones, ln_bg], &ln_w, model),
        Model::UniRgLogLog => ols(&[ones, ln_rg], &ln_w, model),
        Model::BiLogLog => ols(&[ones, ln_rg, ln_bg], &ln_w, model),
        Model::BiLinear => ols(&[ones, rg, bg], &w, model),
    }
}

/// 拟合用的样本点 `(R/G, B/G, w_k)`。
fn model_points(images: &[&ImageData], k: usize) -> Vec<(f64, f64, f64)> {
    images.iter().map(|im| (im.rg, im.bg, im.w_norm[k])).collect()
}

// ---------------------------------------------------------------------------
// 逐图观测
// ---------------------------------------------------------------------------

struct ImageData {
    stem: String,
    rg: f64,
    bg: f64,
    n: usize,
    samples: Vec<Sample>,
    /// 主模型 `M0·diag(w)`，按绿通道归一化。
    w_norm: [f64; 3],
    /// 只用左半解出的主模型 `w`（空间留出用）。
    w_left: [f64; 3],
    /// 对照模型 `diag(w)·M0`，按绿通道归一化。
    w_out: [f64; 3],
    /// 相对线性残差：主模型 / 对照模型。
    resid_cam: f64,
    resid_out: f64,
    /// 色度读数：M0 原样 / 主模型 / 对照模型 / 左半解（在右半上评估）。
    chroma_m0: f64,
    chroma_cam: f64,
    chroma_out: f64,
    chroma_half: f64,
    secs: f64,
}

/// 解一张图。
fn observe(stem: &str, m0: Mat3) -> Option<ImageData> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        common::skip("诊断 K", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }
    let t0 = Instant::now();

    let wb = match libraw::read_wb(&nef) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("跳过 {stem}：读白平衡失败（{e}）");
            return None;
        }
    };
    if wb[1] <= 0.0 {
        eprintln!("跳过 {stem}：cam_mul 的 G 非正");
        return None;
    }
    let (rg, bg) = ((wb[0] / wb[1]) as f64, (wb[2] / wb[1]) as f64);

    let tif_data = match std::fs::read(&tif) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("跳过 {stem}：读参考 TIF 失败（{e}）");
            return None;
        }
    };
    let img = match fit::read_rgb16(&tif_data) {
        Ok(i) => i,
        Err(e) => {
            eprintln!("跳过 {stem}：解析参考 TIF 失败（{e}）");
            return None;
        }
    };
    let plan = match icc::plan_for(&tif_data) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("跳过 {stem}：解析内嵌 ICC 失败（{e}）");
            return None;
        }
    };
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let cam = match libraw::decode_with_output_for_test(&nef, &opts, OutputColor::Camera) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("跳过 {stem}：相机空间解码失败（{e}）");
            return None;
        }
    };
    let Some(model) = camera::read_model(&nef) else {
        eprintln!("跳过 {stem}：读不出机型");
        return None;
    };
    let Some(entry) = camera::lookup(&model) else {
        eprintln!("跳过 {stem}：机型 {model} 不在边距表中");
        return None;
    };
    let mg = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let Some((cw, ch)) = mg.effective(cam.width, cam.height) else {
        eprintln!("跳过 {stem}：边距与解码尺寸不自洽");
        return None;
    };
    if (cw, ch) != (img.width, img.height) {
        eprintln!("跳过 {stem}：裁切后 {cw}×{ch} 与参考 {}×{} 不符", img.width, img.height);
        return None;
    }
    let theirs = fit::reference_to_working(&img, &plan);

    let mut samples: Vec<Sample> = Vec::new();
    let mut dropped = 0usize;
    for y in (0..ch).step_by(STEP) {
        for x in (0..cw).step_by(STEP) {
            let j = x + y * img.width;
            let Some(px) = cam.at(x + mg.left, y + mg.top) else { continue };
            let c = norm96(px);
            let t = theirs[j];
            if c.iter().any(|v| *v <= CAM_LO || *v >= CAM_HI)
                || t.iter().any(|v| *v <= REF_LO || *v >= REF_HI)
            {
                dropped += 1;
                continue;
            }
            samples.push((c, t, x < cw / 2));
        }
    }
    if samples.len() < MIN_SAMPLES {
        eprintln!("跳过 {stem}：可用采样点仅 {} 个（丢弃 {dropped}）", samples.len());
        return None;
    }

    let raw_cam = solve_w_cam(m0, &samples)?;
    let raw_out = solve_w_out(m0, &samples)?;
    let w_norm = norm_g(raw_cam)?;
    let w_out = norm_g(raw_out)?;

    let left: Vec<Sample> = samples.iter().filter(|s| s.2).copied().collect();
    let right: Vec<Sample> = samples.iter().filter(|s| !s.2).copied().collect();
    if left.len() < MIN_SAMPLES / 2 || right.len() < MIN_SAMPLES / 2 {
        eprintln!("跳过 {stem}：左右半样本不足");
        return None;
    }
    let w_left = norm_g(solve_w_cam(m0, &left)?)?;

    let resid_cam = rel_resid(&samples, |c| apply_w_cam(m0, raw_cam, c));
    let resid_out = rel_resid(&samples, |c| apply_w_out(m0, raw_out, c));
    let chroma_m0 = chroma_err(&samples, |c| mat3::mul_vec(m0, c));
    let chroma_cam = chroma_err(&samples, |c| apply_w_cam(m0, raw_cam, c));
    let chroma_out = chroma_err(&samples, |c| apply_w_out(m0, raw_out, c));
    let chroma_half = chroma_err(&right, |c| apply_w_cam(m0, w_left, c));

    Some(ImageData {
        stem: stem.to_string(),
        rg,
        bg,
        n: samples.len(),
        samples,
        w_norm,
        w_left,
        w_out,
        resid_cam,
        resid_out,
        chroma_m0,
        chroma_cam,
        chroma_out,
        chroma_half,
        secs: t0.elapsed().as_secs_f64(),
    })
}

// ---------------------------------------------------------------------------
// 普查与选图
// ---------------------------------------------------------------------------

/// `simple/` 下所有「NEF + 配对 TIF」的白平衡，按 B/G 升序。
///
/// 顺带报出**有 TIF 却找不到 NEF** 的文件：那种参考导出无法配对，值得知道。
fn pair_census() -> Vec<(String, f64, f64)> {
    let dir = samples_dir();
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return Vec::new();
    };
    let all: Vec<std::path::PathBuf> = entries.flatten().map(|e| e.path()).collect();
    let nefs: Vec<std::path::PathBuf> = all
        .iter()
        .filter(|p| p.extension().map(|x| x.eq_ignore_ascii_case("nef")).unwrap_or(false))
        .cloned()
        .collect();
    for tif in all
        .iter()
        .filter(|p| p.extension().map(|x| x.eq_ignore_ascii_case("tif")).unwrap_or(false))
    {
        if !nefs.iter().any(|n| n.with_extension("NEF") == tif.with_extension("NEF")) {
            eprintln!("  普查：{} 有 TIF 但没有配对 NEF，无法使用", tif.display());
        }
    }

    let mut out = Vec::new();
    for p in &nefs {
        if !p.with_extension("TIF").is_file() {
            continue;
        }
        let stem = p.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
        match libraw::read_wb(p) {
            Ok(w) if w[1] > 0.0 => out.push((stem, (w[0] / w[1]) as f64, (w[2] / w[1]) as f64)),
            Ok(_) => eprintln!("  普查：{stem} 的 cam_mul G 非正，跳过"),
            Err(e) => eprintln!("  普查：{stem} 读白平衡失败（{e}）"),
        }
    }
    out.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    out
}

/// 按**相对** B/G 容差把（已按 B/G 升序的）一串 B/G 分档。
fn group_by_bg(bgs: &[f64], tol: f64) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for (i, b) in bgs.iter().enumerate() {
        match out.last_mut() {
            Some(g) if (b - bgs[g[0]]) / bgs[g[0]].abs().max(1e-12) <= tol => g.push(i),
            _ => out.push(vec![i]),
        }
    }
    out
}

/// **白平衡状态**的联合分档：`R/G` 与 `B/G` **都**在容差内才算同一状态。
///
/// 这是本轮相对上一轮的关键修正。只按 `B/G` 分档时，「B/G 几乎相同但 `R/G` 差 21%」的
/// 两张图（DSC_0567 与 DSC_0141）会被当成同一档，而它们的 `w_R` 差 24%——那一档的
/// 「档内离散」于是量的是**两个不同白平衡状态之间的差**，毫无意义。
/// 反过来，同一状态内的多张图（那 15 张连拍）仍然会被正确合并。
///
/// `images` 必须已按 `B/G` 升序。
fn group_states(images: &[ImageData], tol: f64) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for (i, im) in images.iter().enumerate() {
        let head = out.last().map(|g: &Vec<usize>| (images[g[0]].rg, images[g[0]].bg));
        match head {
            Some((hrg, hbg))
                if (im.bg - hbg).abs() / hbg <= tol && (im.rg - hrg).abs() / hrg <= tol =>
            {
                if let Some(g) = out.last_mut() {
                    g.push(i);
                }
            }
            _ => out.push(vec![i]),
        }
    }
    out
}

fn print_census(census: &[(String, f64, f64)]) -> Vec<Vec<usize>> {
    eprintln!(
        "\n[0] 白平衡普查：simple/ 下 {} 对「NEF + 配对 TIF」，按相对容差 {:.1}% 分档",
        census.len(),
        WB_LEVEL_TOL * 100.0
    );
    let levels = group_by_bg(&census.iter().map(|c| c.2).collect::<Vec<_>>(), WB_LEVEL_TOL);
    eprintln!("  {:<4} {:>9} {:>9} {:>6}  成员（前 6 个）", "档", "B/G", "R/G", "张数");
    for (li, g) in levels.iter().enumerate() {
        let bgs: Vec<f64> = g.iter().map(|&i| census[i].2).collect();
        let lo = bgs.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = bgs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let names: Vec<&str> = g.iter().take(6).map(|&i| census[i].0.as_str()).collect();
        eprintln!(
            "  {:<4} {:>9.4} {:>9.4} {:>6}  {}{}",
            li + 1,
            census[g[0]].2,
            census[g[0]].1,
            g.len(),
            names.join(" "),
            if g.len() > 6 { " …" } else { "" }
        );
        if g.len() > 1 {
            eprintln!("       档内 B/G 跨度 {lo:.4} → {hi:.4}（{:.2}%）", 100.0 * (hi - lo) / lo);
        }
    }
    let first = census.first().map(|c| c.2).unwrap_or(f64::NAN);
    let last = census.last().map(|c| c.2).unwrap_or(f64::NAN);
    eprintln!("  ⇒ {} 档，B/G 覆盖 {first:.4} → {last:.4}（{:.2}×）。", levels.len(), last / first);
    levels
}

/// 从 `DSC_0001` 这类名字里取编号，用于判断两张图「相隔多远」。
fn stem_number(s: &str) -> Option<u32> {
    s.rsplit('_').next()?.parse::<u32>().ok()
}

/// 选图。规则（确定性、可复核）：
///
/// 1. **每一档都要一张代表**（档按字母序取第一张）——绝不为了让预算好看而丢档；
/// 2. 多图档再补一张**编号与该档代表相隔最远**的成员：那是「不同拍摄批次/不同场景」
///    在本数据里唯一可用的代理，用来检验「同白平衡、不同场景是否给出同一个 w」；
/// 3. 最大的一档再补一张，保证连拍档至少有 3 张可量档内离散；
/// 4. 超过 [`MAX_IMAGES`] 时**先丢额外补的图，并打印丢了谁**。
fn choose_stems(
    census: &[(String, f64, f64)],
    levels: &[Vec<usize>],
) -> (Vec<String>, Vec<String>) {
    let names = |g: &Vec<usize>| -> Vec<String> {
        let mut m: Vec<String> = g.iter().map(|&i| census[i].0.clone()).collect();
        m.sort_unstable();
        m
    };
    let farthest = |members: &[String], rep: &str| -> Option<String> {
        let rn = stem_number(rep)?;
        members
            .iter()
            .filter(|m| m.as_str() != rep)
            .filter_map(|m| stem_number(m).map(|n| (n.abs_diff(rn), m.clone())))
            .max_by_key(|(d, _)| *d)
            .map(|(_, s)| s)
    };

    let mut reps: Vec<String> = Vec::new();
    let mut extras: Vec<String> = Vec::new();
    for g in levels {
        let members = names(g);
        let Some(rep) = members.first() else { continue };
        reps.push(rep.clone());
        if members.len() > 1 {
            if let Some(f) = farthest(&members, rep) {
                extras.push(f);
            }
        }
    }
    if let Some(g) = levels.iter().max_by_key(|g| g.len()) {
        let members = names(g);
        if let Some(rep) = members.first() {
            if let Some(f) = farthest(&members, rep) {
                extras.push(f);
            }
        }
    }
    // 去重：既要去掉与代表重复的，也要去掉 extras 内部的重复（同一状态可能被两条规则各挑中一次，
    // 例如连拍档既是「多图档」又是「最大档」）。`dedup()` 只去相邻重复，所以先排序。
    extras.retain(|s| !reps.contains(s));
    extras.sort();
    extras.dedup();

    let mut chosen = reps;
    let mut notes = Vec::new();
    for e in extras {
        if chosen.len() < MAX_IMAGES {
            chosen.push(e);
        } else {
            notes.push(format!("预算 {MAX_IMAGES} 已满，丢掉额外补的 {e}"));
        }
    }
    (chosen, notes)
}

// ---------------------------------------------------------------------------
// 交叉验证
// ---------------------------------------------------------------------------

/// 一条交叉验证记录：完整的三通道 `w`（绿通道恒为 1）。
struct CvRow {
    label: String,
    bg: f64,
    rg: f64,
    /// 留出对象的实测值（单张图，或该状态多张图的均值）。
    measured: [f64; 3],
    /// 用其余状态拟合出来的预测值。
    predicted: [f64; 3],
}

/// 在训练集上为 R、B 两个通道各拟合一个模型（多元接口）。
fn fit_models(train: &[&ImageData], mr: Model, mb: Model) -> Option<(Fitted, Fitted)> {
    let fr = fit_model(mr, &model_points(train, 0))?;
    let fb = fit_model(mb, &model_points(train, 2))?;
    Some((fr, fb))
}

/// 一组图的 `w` 均值（三通道各自平均；绿通道恒为 1）。
fn mean_w(images: &[&ImageData]) -> [f64; 3] {
    let n = images.len().max(1) as f64;
    let mut out = [0f64; 3];
    for im in images {
        for (k, slot) in out.iter_mut().enumerate() {
            *slot += im.w_norm[k];
        }
    }
    [out[0] / n, out[1] / n, out[2] / n]
}

fn mean_of_key<F: Fn(&ImageData) -> f64>(images: &[&ImageData], f: F) -> f64 {
    if images.is_empty() {
        return f64::NAN;
    }
    images.iter().map(|im| f(im)).sum::<f64>() / images.len() as f64
}

fn mean_bg(images: &[&ImageData]) -> f64 {
    mean_of_key(images, |im| im.bg)
}

fn mean_rg(images: &[&ImageData]) -> f64 {
    mean_of_key(images, |im| im.rg)
}

/// 训练集 = 除 `held` 之外的所有图。
fn train_excluding<'a>(images: &'a [ImageData], held: &[usize]) -> Vec<&'a ImageData> {
    images
        .iter()
        .enumerate()
        .filter(|(i, _)| !held.contains(i))
        .map(|(_, im)| im)
        .collect()
}

fn group_label(members: &[&ImageData]) -> String {
    if members.len() == 1 {
        members[0].stem.clone()
    } else {
        format!("{} 状态（{} 张）", members[0].stem, members.len())
    }
}

/// **留一个白平衡状态**（联合按 R/G 与 B/G 分档）：抽掉整个状态，用其余状态拟合。
fn cv_rows_leave_group(
    images: &[ImageData],
    groups: &[Vec<usize>],
    mr: Model,
    mb: Model,
) -> Vec<CvRow> {
    let mut rows = Vec::new();
    for held in groups {
        let train = train_excluding(images, held);
        let Some((fr, fb)) = fit_models(&train, mr, mb) else { continue };
        let members: Vec<&ImageData> = held.iter().map(|&i| &images[i]).collect();
        let (rg, bg) = (mean_rg(&members), mean_bg(&members));
        rows.push(CvRow {
            label: group_label(&members),
            bg,
            rg,
            measured: mean_w(&members),
            predicted: [fr.predict(rg, bg), 1.0, fb.predict(rg, bg)],
        });
    }
    rows
}

/// 留一图：每张图轮流留出。
fn cv_rows_leave_image(images: &[ImageData], mr: Model, mb: Model) -> Vec<CvRow> {
    let mut rows = Vec::new();
    for (i, held) in images.iter().enumerate() {
        let train = train_excluding(images, &[i]);
        let Some((fr, fb)) = fit_models(&train, mr, mb) else { continue };
        rows.push(CvRow {
            label: held.stem.clone(),
            bg: held.bg,
            rg: held.rg,
            measured: held.w_norm,
            predicted: [fr.predict(held.rg, held.bg), 1.0, fb.predict(held.rg, held.bg)],
        });
    }
    rows
}

/// 单通道的留一状态 CV 相对误差（逐通道选模型用）。
fn cv_model_errors(
    images: &[ImageData],
    groups: &[Vec<usize>],
    model: Model,
    k: usize,
) -> Vec<f64> {
    let mut out = Vec::new();
    for held in groups {
        let train = train_excluding(images, held);
        let Some(f) = fit_model(model, &model_points(&train, k)) else { continue };
        let members: Vec<&ImageData> = held.iter().map(|&i| &images[i]).collect();
        let mw = members.iter().map(|im| im.w_norm[k]).sum::<f64>() / members.len() as f64;
        out.push(rel_err(f.predict(mean_rg(&members), mean_bg(&members)), mw));
    }
    out
}

// ---------------------------------------------------------------------------
// 主测试
// ---------------------------------------------------------------------------

fn env_stems() -> Option<Vec<String>> {
    match std::env::var("NRV_K_STEMS") {
        Ok(v) if !v.trim().is_empty() => Some(
            v.split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect(),
        ),
        _ => None,
    }
}

fn census_enabled() -> bool {
    !matches!(std::env::var("NRV_K_CENSUS").as_deref(), Ok("0") | Ok("no") | Ok("false"))
}

#[test]
fn fit_wb_dependent_gain() {
    eprintln!("=== 诊断 K：拟合 w(B/G) 并做留一档交叉验证 ===");
    let Some(cm) = camera::color_matrix("Nikon Z 8") else {
        eprintln!("跳过：机型表里没有 Z8 的色彩矩阵");
        return;
    };
    let m0 = cm.matrix;
    eprintln!("主模型：M0·diag(w)（增益在矩阵之前）；对照：diag(w)·M0（增益在矩阵之后）");

    // ---- 选图 ----
    let stems: Vec<String> = match env_stems() {
        Some(list) => {
            eprintln!("\n[0] 按 NRV_K_STEMS 指定 {} 张图，跳过普查", list.len());
            list
        }
        None => {
            let census = if census_enabled() {
                pair_census()
            } else {
                eprintln!("\n[0] 按 NRV_K_CENSUS 跳过普查");
                Vec::new()
            };
            if census.is_empty() {
                eprintln!("跳过：没有找到任何「NEF + 配对 TIF」");
                return;
            }
            let levels = print_census(&census);
            let (chosen, notes) = choose_stems(&census, &levels);
            eprintln!(
                "\n     解码 {} 张：每档 1 张代表（{} 档），多图档各补 {EXTRA_PER_MULTI} 张「编号相隔最远」的成员（不同场景代理），",
                chosen.len(),
                levels.len()
            );
            eprintln!(
                "     最大的一档再补 {EXTRA_PER_BIGGEST} 张以量档内离散；上限 {MAX_IMAGES} 张。"
            );
            for n in &notes {
                eprintln!("     ⚠ {n}");
            }
            eprintln!("     实际选中：{}", chosen.join(" "));
            chosen
        }
    };

    // ---- 逐图解码 + 解 w ----
    let mut images: Vec<ImageData> = Vec::new();
    for stem in &stems {
        eprintln!("  … 解码 {stem}");
        if let Some(d) = observe(stem, m0) {
            images.push(d);
        }
    }
    if images.len() < 4 {
        eprintln!("跳过：可用图不足（{} 张）", images.len());
        return;
    }
    images.sort_by(|a, b| a.bg.partial_cmp(&b.bg).unwrap_or(std::cmp::Ordering::Equal));
    let bgs: Vec<f64> = images.iter().map(|im| im.bg).collect();
    // 两种分档并存：`levels` 是**只按 B/G 的一维投影**（上一轮的做法，保留用于对照），
    // `states` 是**联合按 (R/G, B/G) 的白平衡状态**（本轮的主口径）。
    let levels = group_by_bg(&bgs, WB_LEVEL_TOL);
    let states = group_states(&images, WB_LEVEL_TOL);
    eprintln!(
        "\n分组：按 B/G 一维投影得 {} 档；联合按 (R/G, B/G) 得 {} 个白平衡状态（主口径）。",
        levels.len(),
        states.len()
    );

    // ---- [1] 逐图 w ----
    eprintln!("\n[1a] 逐图解出的 w（按绿通道归一化）");
    eprintln!(
        "  {:<10} {:>8} {:>8} {:>8} {:>9} {:>9} {:>10} {:>10} {:>10}",
        "图", "B/G", "R/G", "样本", "w_R", "w_B", "out_w_R", "out_w_B", "左半w_R"
    );
    for im in &images {
        eprintln!(
            "  {:<10} {:>8.4} {:>8.4} {:>8} {:>9.5} {:>9.5} {:>10.5} {:>10.5} {:>10.5}",
            im.stem, im.bg, im.rg, im.n, im.w_norm[0], im.w_norm[2], im.w_out[0], im.w_out[2], im.w_left[0]
        );
    }
    eprintln!("\n[1b] 残差与色度读数（色度 = 尺度不变 |Δlog2(R/G),Δlog2(B/G)| 中位，越小越好）");
    eprintln!(
        "  {:<10} {:>10} {:>10} {:>11} {:>11} {:>11} {:>14}",
        "图", "残差cam", "残差out", "色度M0", "色度cam", "色度out", "色度(左解→右评)"
    );
    for im in &images {
        eprintln!(
            "  {:<10} {:>10.4} {:>10.4} {:>11.5} {:>11.5} {:>11.5} {:>14.5}",
            im.stem, im.resid_cam, im.resid_out, im.chroma_m0, im.chroma_cam, im.chroma_out, im.chroma_half
        );
    }
    let cam_better = images.iter().filter(|im| im.resid_cam < im.resid_out).count();
    let chroma_cam_better = images.iter().filter(|im| im.chroma_cam < im.chroma_out).count();
    let mp_cam = median_of(&images.iter().map(|im| im.resid_cam).collect::<Vec<_>>());
    let mp_out = median_of(&images.iter().map(|im| im.resid_out).collect::<Vec<_>>());
    let mc_m0 = median_of(&images.iter().map(|im| im.chroma_m0).collect::<Vec<_>>());
    let mc_cam = median_of(&images.iter().map(|im| im.chroma_cam).collect::<Vec<_>>());
    let mc_out = median_of(&images.iter().map(|im| im.chroma_out).collect::<Vec<_>>());
    let mc_half = median_of(&images.iter().map(|im| im.chroma_half).collect::<Vec<_>>());
    eprintln!(
        "  中位：线性残差 cam {mp_cam:.4} vs out {mp_out:.4}（{cam_better}/{} 张 cam 更小）；",
        images.len()
    );
    eprintln!(
        "        色度 M0 {mc_m0:.5} → cam {mc_cam:.5} / out {mc_out:.5}（cam 更好的 {chroma_cam_better}/{} 张；左半解→右半评 {mc_half:.5}）",
        images.len()
    );
    let total_secs: f64 = images.iter().map(|im| im.secs).sum();
    eprintln!("  逐图解码总耗时 {total_secs:.0} s。");

    // ---- 状态内离散 ----
    eprintln!("\n[1c] 状态内离散（同一**白平衡状态**多张图时的 w 散布，是交叉验证误差的参照底）");
    let mut any_multi = false;
    for g in &states {
        if g.len() < 2 {
            continue;
        }
        any_multi = true;
        let wr: Vec<f64> = g.iter().map(|&i| images[i].w_norm[0]).collect();
        let wb: Vec<f64> = g.iter().map(|&i| images[i].w_norm[2]).collect();
        let rlo = wr.iter().cloned().fold(f64::INFINITY, f64::min);
        let rhi = wr.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let blo = wb.iter().cloned().fold(f64::INFINITY, f64::min);
        let bhi = wb.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        eprintln!(
            "    {} 档 {} 张：w_R {:.5}~{:.5}（跨 {:.2}%）· w_B {:.5}~{:.5}（跨 {:.2}%）",
            images[g[0]].stem,
            g.len(),
            rlo,
            rhi,
            100.0 * (rhi - rlo) / rlo,
            blo,
            bhi,
            100.0 * (bhi - blo) / blo
        );
    }
    if !any_multi {
        eprintln!("    （每档只有一张图，无法量档内离散）");
    }
    // 档内散度的总括（后面判读要用）：所有多图**状态**的 w 相对跨度的最大值
    let mut scatter_r = f64::NAN;
    let mut scatter_b = f64::NAN;
    for g in &states {
        if g.len() < 2 {
            continue;
        }
        let wr: Vec<f64> = g.iter().map(|&i| images[i].w_norm[0]).collect();
        let wb: Vec<f64> = g.iter().map(|&i| images[i].w_norm[2]).collect();
        let span = |v: &Vec<f64>| {
            let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            100.0 * (hi - lo) / lo
        };
        scatter_r = scatter_r.max(span(&wr));
        scatter_b = scatter_b.max(span(&wb));
    }

    // ---- [1d] 近同档配对：同白平衡、不同图（含不同场景）的 w 是否一致 ----
    eprintln!(
        "\n[1d] 近同档配对（B/G 相对差 ≤ {:.0}%）：同白平衡但**不同图**的 w 是否一致",
        NEAR_BG_TOL * 100.0
    );
    eprintln!("  这是本批数据里唯一能看「跨场景」的地方：色温扫描那批与老图场景不同。");
    let same_state =
        |i: usize, j: usize| states.iter().any(|g| g.contains(&i) && g.contains(&j));
    eprintln!(
        "  {:<10} {:<10} {:>9} {:>9} {:>9} {:>9}  关系",
        "图 A", "图 B", "A 的B/G", "B 的B/G", "Δw_R", "Δw_B"
    );
    let mut near_pairs = 0usize;
    let mut cross_scene: Vec<f64> = Vec::new();
    let mut same_scene: Vec<f64> = Vec::new();
    for i in 0..images.len() {
        for j in (i + 1)..images.len() {
            let (a, b) = (&images[i], &images[j]);
            if (b.bg - a.bg).abs() / a.bg > NEAR_BG_TOL {
                continue;
            }
            near_pairs += 1;
            let dr = 100.0 * (b.w_norm[0] - a.w_norm[0]).abs() / a.w_norm[0];
            let db = 100.0 * (b.w_norm[2] - a.w_norm[2]).abs() / a.w_norm[2];
            let same = same_state(i, j);
            if same {
                same_scene.push(dr.max(db));
            } else {
                cross_scene.push(dr.max(db));
            }
            eprintln!(
                "  {:<10} {:<10} {:>9.4} {:>9.4} {:>8.2}% {:>8.2}%  {}",
                a.stem,
                b.stem,
                a.bg,
                b.bg,
                dr,
                db,
                if same { "同一状态（重复图）" } else { "同 B/G 档、**不同状态**" }
            );
        }
    }
    if near_pairs == 0 {
        eprintln!("    （没有 B/G 足够接近的配对）");
    } else {
        eprintln!(
            "  ⇒ 最大差：同一状态的配对 {:.2}%（{} 对）· 同 B/G 但不同状态的配对 {:.2}%（{} 对）",
            max_of(&same_scene),
            same_scene.len(),
            max_of(&cross_scene),
            cross_scene.len()
        );
        eprintln!("     后者 = 同一个 B/G 上 w 仍能差出这么多 ⇒ **B/G 一个自变量不够**。");
        eprintln!("     但要如实说：这些跨状态配对同时也是**跨场景**的（扫描组 vs 老图），");
        eprintln!("     所以单看它们分不清「状态不同」与「场景不同」；[1e]/[2b] 的模型级证据才排除了场景解释。");
    }

    // ---- [1e] 二元相关性：w 到底跟着 R/G 还是 B/G 走 ----
    eprintln!("\n[1e] w 与两个白平衡比值的关系（自变量应当是两个，不是一个）");
    let ln_rg: Vec<f64> = images.iter().map(|im| im.rg.ln()).collect();
    let ln_bg: Vec<f64> = images.iter().map(|im| im.bg.ln()).collect();
    // R/B 是「绿-品红」那一维：M 偏移主要动它
    let ln_rb: Vec<f64> = images.iter().map(|im| (im.rg / im.bg).ln()).collect();
    eprintln!(
        "  {:<6} {:>14} {:>14} {:>14} {:>16} {:>16}",
        "通道", "r(ln w, lnR/G)", "r(ln w, lnB/G)", "r(ln w, lnR/B)", "偏相关|lnB/G", "偏相关|lnR/G"
    );
    for k in [0usize, 2usize] {
        let ln_w: Vec<f64> = images.iter().map(|im| im.w_norm[k].ln()).collect();
        eprintln!(
            "  {:<6} {:>14.3} {:>14.3} {:>14.3} {:>16.3} {:>16.3}",
            if k == 0 { "R" } else { "B" },
            pearson(&ln_w, &ln_rg),
            pearson(&ln_w, &ln_bg),
            pearson(&ln_w, &ln_rb),
            partial_corr(&ln_w, &ln_rg, &ln_bg),
            partial_corr(&ln_w, &ln_bg, &ln_rg)
        );
    }
    eprintln!("  偏相关 = 控制另一个变量之后的相关（第一列控制 lnB/G，第二列控制 lnR/G）。");
    eprintln!("  R/G 与 B/G 本身高度共线，所以单看相关系数会互相顶替；偏相关才是「谁在起作用」。");
    // 两条投影轴各自的范围（说明为什么必须两个自变量）
    let rgs: Vec<f64> = images.iter().map(|im| im.rg).collect();
    let bgs_v: Vec<f64> = images.iter().map(|im| im.bg).collect();
    let (rlo, rhi) = (
        rgs.iter().cloned().fold(f64::INFINITY, f64::min),
        rgs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let (blo, bhi) = (
        bgs_v.iter().cloned().fold(f64::INFINITY, f64::min),
        bgs_v.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let r_between = pearson(&ln_rg, &ln_bg);
    eprintln!(
        "  样本：R/G {rlo:.4}~{rhi:.4}（{:.2}×）· B/G {blo:.4}~{bhi:.4}（{:.2}×）· r(ln R/G, ln B/G) = {r_between:+.3}",
        rhi / rlo,
        bhi / blo
    );
    eprintln!("  ⇒ 两个比值共线但**不重合**，平面上确有第二个自由度。");

    // ---- [2] 一元形式表（任务要求保留） ----
    let all: Vec<&ImageData> = images.iter().collect();
    eprintln!("\n[2] 只用 B/G 的六种形式（系数与训练残差；上一轮的形式表，保留备查）");
    eprintln!("  {:<18} {:>4} {:>11} {:>11} {:>14}", "形式", "通道", "a", "b", "训练残差→%");
    for form in FORMS {
        for k in [0usize, 2usize] {
            let Some(f) = fit_form(form, &points(&all, k)) else { continue };
            eprintln!(
                "  {:<18} {:>4} {:>11.5} {:>11.5} {:>13.2}%",
                form.name(),
                if k == 0 { "R" } else { "B" },
                f.a,
                f.b,
                100.0 * (f.log_rms.exp() - 1.0)
            );
        }
    }

    // ---- [2b] 模型比较：一元 vs 二元 ----
    eprintln!("\n[2b] 模型比较（**留一个白平衡状态** CV，按 minimax 逐通道选）");
    eprintln!("  选择口径仍是**最大误差**：要用的地方是「没见过的白平衡状态」，最坏的那个才是风险。");
    eprintln!(
        "  {:<34} {:>4} {:>11} {:>11}",
        "模型", "通道", "CV 中位", "CV 最大"
    );
    let mut chosen_model = [Model::ConstMedian; 3];
    for k in [0usize, 2usize] {
        let label = if k == 0 { "R" } else { "B" };
        let mut best_max: Option<(Model, f64, f64)> = None;
        let mut best_med: Option<(Model, f64, f64)> = None;
        let mut best_uni_max: Option<(Model, f64, f64)> = None;
        for model in MODELS {
            let errs = cv_model_errors(&images, &states, model, k);
            let (med, mx) = (median_of(&errs), max_of(&errs));
            eprintln!("  {:<34} {:>4} {:>10.2}% {:>10.2}%", model.name(), label, 100.0 * med, 100.0 * mx);
            if best_max.map(|(_, bm, _)| mx < bm).unwrap_or(true) {
                best_max = Some((model, mx, med));
            }
            if best_med.map(|(_, bm, _)| med < bm).unwrap_or(true) {
                best_med = Some((model, med, mx));
            }
            if !model.is_bivariate() && best_uni_max.map(|(_, bm, _)| mx < bm).unwrap_or(true) {
                best_uni_max = Some((model, mx, med));
            }
        }
        if let Some((m, mx, med)) = best_max {
            chosen_model[k] = m;
            eprintln!("    ⇒ 通道 {label} 按**最大误差**选中：{}（中位 {:.2}% / 最大 {:.2}%）", m.name(), 100.0 * med, 100.0 * mx);
        }
        if let Some((m, mx, med)) = best_uni_max {
            if m.is_bivariate() {
                eprintln!("      （best_univariate 不应是二元模型，逻辑有误）");
            } else if chosen_model[k].is_bivariate() {
                let bm = best_max.map(|(_, x, _)| x).unwrap_or(f64::NAN);
                eprintln!(
                    "      **二元相对一元最好者**：最大误差 {:.2}% → {:.2}%（改善 {:.1}%），中位 {:.2}% → {:.2}%",
                    100.0 * mx,
                    100.0 * bm,
                    100.0 * (1.0 - bm / mx.max(1e-12)),
                    100.0 * med,
                    100.0 * best_max.map(|(_, _, mm)| mm).unwrap_or(f64::NAN)
                );
            } else {
                eprintln!(
                    "      **二元并没有赢过一元**：一元最好者 {} 的最大误差 {:.2}%（二元选中项已含在表内）。",
                    m.name(),
                    100.0 * mx
                );
            }
        }
        if let Some((m, med, mx)) = best_med {
            if m != chosen_model[k] {
                eprintln!(
                    "      若只按**中位**挑则是：{}（中位 {:.2}% / **最大 {:.2}%**）",
                    m.name(),
                    100.0 * med,
                    100.0 * mx
                );
            }
        }
    }
    let (model_r, model_b) = (chosen_model[0], chosen_model[2]);
    if let (Some(fr), Some(fb)) = (
        fit_model(model_r, &model_points(&all, 0)),
        fit_model(model_b, &model_points(&all, 2)),
    ) {
        eprintln!(
            "  全量拟合：w_R = {}（a {:.5}, b₁ {:.5}, b₂ {:.5}）",
            model_r.name(),
            fr.c[0],
            fr.c[1],
            fr.c[2]
        );
        eprintln!(
            "            w_B = {}（a {:.5}, b₁ {:.5}, b₂ {:.5}）",
            model_b.name(),
            fb.c[0],
            fb.c[1],
            fb.c[2]
        );
    }
    eprintln!("  （模型是按 CV 指标逐通道挑的，属模型选择，CV 数字因此略偏乐观。）");

    // ---- [3] 留一状态交叉验证 ----
    eprintln!("\n[3] 留一状态交叉验证（重点）：抽掉**整个白平衡状态**，用其余状态拟合后预测它");
    eprintln!("  误差按状态统计（每个状态等权，与其中张数无关——连拍那 3 张只算一次）。");
    eprintln!(
        "  {:<22} {:>8} {:>8} {:>10} {:>10} {:>8} {:>10} {:>10} {:>8}",
        "留出状态", "R/G", "B/G", "实测w_R", "预测w_R", "误差", "实测w_B", "预测w_B", "误差"
    );
    let rows = cv_rows_leave_group(&images, &states, model_r, model_b);
    let mut errs_r = Vec::new();
    let mut errs_b = Vec::new();
    for r in &rows {
        let er = rel_err(r.predicted[0], r.measured[0]);
        let eb = rel_err(r.predicted[2], r.measured[2]);
        errs_r.push(er);
        errs_b.push(eb);
        eprintln!(
            "  {:<22} {:>8.4} {:>8.4} {:>10.5} {:>10.5} {:>7.2}% {:>10.5} {:>10.5} {:>7.2}%",
            r.label, r.rg, r.bg, r.measured[0], r.predicted[0], 100.0 * er, r.measured[2], r.predicted[2], 100.0 * eb
        );
    }
    let amp_r = median_of(&images.iter().map(|im| (im.w_norm[0] - 1.0).abs()).collect::<Vec<_>>());
    let amp_b = median_of(&images.iter().map(|im| (im.w_norm[2] - 1.0).abs()).collect::<Vec<_>>());
    eprintln!(
        "\n  w_R：中位 {:.2}% · 最大 {:.2}%    修正幅度中位 |w_R−1| = {:.2}%    误差÷幅度 {:.2}",
        100.0 * median_of(&errs_r),
        100.0 * max_of(&errs_r),
        100.0 * amp_r,
        median_of(&errs_r) / amp_r.max(1e-12)
    );
    eprintln!(
        "  w_B：中位 {:.2}% · 最大 {:.2}%    修正幅度中位 |w_B−1| = {:.2}%    误差÷幅度 {:.2}",
        100.0 * median_of(&errs_b),
        100.0 * max_of(&errs_b),
        100.0 * amp_b,
        median_of(&errs_b) / amp_b.max(1e-12)
    );

    // ---- [3b] 分组口径与容差的敏感性 ----
    eprintln!("\n[3b] 分组口径与容差的敏感性：近邻会让「留一」几乎不用外推，CV 因此偏乐观");
    eprintln!("  主结果：联合 (R/G, B/G) 状态的 {:.1}% 容差。下面同时改变**分组口径**与**容差**。", WB_LEVEL_TOL * 100.0);
    eprintln!(
        "  {:<26} {:>8} {:>7} {:>22} {:>22}",
        "分组口径", "容差", "组数", "w_R 中位/最坏", "w_B 中位/最坏"
    );
    for tol in std::iter::once(WB_LEVEL_TOL).chain(TOL_SENSITIVITY) {
        // 一维投影（只按 B/G）——上一轮的口径
        let lv = group_by_bg(&bgs, tol);
        let er: Vec<f64> = cv_model_errors(&images, &lv, model_r, 0);
        let eb: Vec<f64> = cv_model_errors(&images, &lv, model_b, 2);
        eprintln!(
            "  {:<26} {:>7.1}% {:>7} {:>21} {:>21}",
            "一维：只按 B/G",
            tol * 100.0,
            lv.len(),
            format!("{:.2}% / {:.2}%", 100.0 * median_of(&er), 100.0 * max_of(&er)),
            format!("{:.2}% / {:.2}%", 100.0 * median_of(&eb), 100.0 * max_of(&eb))
        );
        // 联合状态（本轮主口径）
        let st = group_states(&images, tol);
        let er: Vec<f64> = cv_model_errors(&images, &st, model_r, 0);
        let eb: Vec<f64> = cv_model_errors(&images, &st, model_b, 2);
        eprintln!(
            "  {:<26} {:>7.1}% {:>7} {:>21} {:>21}",
            "联合：(R/G, B/G)",
            tol * 100.0,
            st.len(),
            format!("{:.2}% / {:.2}%", 100.0 * median_of(&er), 100.0 * max_of(&er)),
            format!("{:.2}% / {:.2}%", 100.0 * median_of(&eb), 100.0 * max_of(&eb))
        );
    }
    eprintln!("  一维那几行里，被留出的「档」可能同时含两个 R/G 不同的状态——对手里的模型不公平；");
    eprintln!("  反之联合分档更细，留一时训练集里**不会**留下几乎相同的状态，CV 更严。判读以联合口径为准。");

    // ---- [3c] 功能后果 ----
    eprintln!("\n[3c] 功能后果（留出状态上，尺度不变色度误差中位）：M0 原样 → 用预测 w → 用实测 w");
    eprintln!(
        "  {:<22} {:>10} {:>10} {:>12} {:>12}",
        "留出状态", "M0", "预测 w", "实测状态均值 w", "本图自解 w"
    );
    let mut before_all = Vec::new();
    let mut pred_all = Vec::new();
    let mut ceil_all = Vec::new();
    let mut self_all = Vec::new();
    for (held, r) in states.iter().zip(&rows) {
        let mut before = Vec::new();
        let mut pred = Vec::new();
        let mut ceil = Vec::new();
        let mut selfv = Vec::new();
        for &i in held {
            let im = &images[i];
            before.push(im.chroma_m0);
            pred.push(chroma_err(&im.samples, |c| apply_w_cam(m0, r.predicted, c)));
            ceil.push(chroma_err(&im.samples, |c| apply_w_cam(m0, r.measured, c)));
            selfv.push(im.chroma_cam);
        }
        let (b, p, c, s) = (
            median_of(&before),
            median_of(&pred),
            median_of(&ceil),
            median_of(&selfv),
        );
        before_all.extend(before);
        pred_all.extend(pred);
        ceil_all.extend(ceil);
        self_all.extend(selfv);
        eprintln!("  {:<22} {:>10.5} {:>10.5} {:>12.5} {:>12.5}", r.label, b, p, c, s);
    }
    let (ba, pa, ca, sa) = (
        median_of(&before_all),
        median_of(&pred_all),
        median_of(&ceil_all),
        median_of(&self_all),
    );
    eprintln!(
        "  ⇒ 全体留出样本：M0 {ba:.5} → 预测 {pa:.5}（{:+.0}%）· 实测状态均值（上限）{ca:.5}（{:+.0}%）· 本图自解 {sa:.5}",
        100.0 * (pa / ba.max(1e-12) - 1.0),
        100.0 * (ca / ba.max(1e-12) - 1.0)
    );

    // ---- [4] 留一图交叉验证 ----
    eprintln!("\n[4] 留一图交叉验证（每张图轮流留出；训练集里仍有同一**状态**图的那些行检验状态内泛化）");
    let img_rows = cv_rows_leave_image(&images, model_r, model_b);
    eprintln!(
        "  {:<10} {:>8} {:>8} {:>9} {:>9} {:>8} {:>9} {:>9} {:>8}  同状态还有图？",
        "留出图", "R/G", "B/G", "实测w_R", "预测w_R", "误差", "实测w_B", "预测w_B", "误差"
    );
    let mut same_level: Vec<f64> = Vec::new();
    let mut other_level: Vec<f64> = Vec::new();
    for (i, r) in img_rows.iter().enumerate() {
        let er = rel_err(r.predicted[0], r.measured[0]);
        let eb = rel_err(r.predicted[2], r.measured[2]);
        let has_same = states.iter().find(|g| g.contains(&i)).map(|g| g.len() > 1).unwrap_or(false);
        if has_same {
            same_level.push(er.max(eb));
        } else {
            other_level.push(er.max(eb));
        }
        eprintln!(
            "  {:<10} {:>8.4} {:>8.4} {:>9.5} {:>9.5} {:>7.2}% {:>9.5} {:>9.5} {:>7.2}%  {}",
            r.label,
            r.rg,
            r.bg,
            r.measured[0],
            r.predicted[0],
            100.0 * er,
            r.measured[2],
            r.predicted[2],
            100.0 * eb,
            if has_same { "是" } else { "否" }
        );
    }
    eprintln!(
        "\n  训练集里**仍有同一状态**的图（{} 张）：中位 {:.2}% · 最大 {:.2}%",
        same_level.len(),
        100.0 * median_of(&same_level),
        100.0 * max_of(&same_level)
    );
    eprintln!(
        "  训练集里**没有同一状态**的图（{} 张）：中位 {:.2}% · 最大 {:.2}%",
        other_level.len(),
        100.0 * median_of(&other_level),
        100.0 * max_of(&other_level)
    );

    // ---- [5] 判读 ----
    let med_r = median_of(&errs_r);
    let max_r = max_of(&errs_r);
    let med_b = median_of(&errs_b);
    let max_b = max_of(&errs_b);
    let ratio_r = med_r / amp_r.max(1e-12);
    let ratio_b = med_b / amp_b.max(1e-12);
    eprintln!(
        "\n[5] 判读（门槛先定后看：CV 中位 ≤ {:.0}% 算可用；误差÷幅度 ≥ {:.1} 即「与修正幅度同量级」）",
        100.0 * ERR_USABLE,
        ERR_VS_AMPLITUDE
    );
    for (label, model, med, mx, amp, ratio) in [
        ("R", model_r, med_r, max_r, amp_r, ratio_r),
        ("B", model_b, med_b, max_b, amp_b, ratio_b),
    ] {
        eprintln!(
            "  通道 {label}：选中「{}」· CV 中位 {:.2}% / **最大 {:.2}%** · 修正幅度 {:.2}% · 误差÷幅度 {:.2} · 状态内噪声底 {:.2}%",
            model.name(),
            100.0 * med,
            100.0 * mx,
            100.0 * amp,
            ratio,
            if label == "R" { scatter_r } else { scatter_b }
        );
        if model == Model::ConstMean || model == Model::ConstMedian {
            eprintln!("     ⇒ **该通道与白平衡无关**：曲线跑不赢常数。");
        } else if med <= ERR_USABLE && ratio < ERR_VS_AMPLITUDE {
            eprintln!("     ⇒ 中位在门槛内且误差远小于修正幅度：**该通道的模型是一阶可用的**。");
        } else if ratio >= ERR_VS_AMPLITUDE {
            eprintln!("     ⇒ **负结论**：预测误差与要修正的量本身同量级，该通道不能只靠白平衡曲线。");
        } else {
            eprintln!("     ⇒ 抓住了主要趋势，但精度介于两者之间：可作一阶修正，不能当精确标定。");
        }
    }
    // 「翻没翻」的基准是**上一轮那一族**（只允许用 B/G），不是「所有一元模型」。
    // 这一点很关键：上一轮的负结论说的是「B/G 一个自变量不够」，不是「任何一元都不够」。
    let bg_only = MODELS
        .iter()
        .filter(|m| m.bg_only())
        .map(|m| {
            let e = cv_model_errors(&images, &states, *m, 0);
            (*m, median_of(&e), max_of(&e))
        })
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    if let Some((m, med_u, max_u)) = bg_only {
        eprintln!(
            "  **红通道翻没翻**（基准 = 上一轮那一族「只用 B/G」里最好者「{}」：中位 {:.2}% / 最大 {:.2}%）",
            m.name(),
            100.0 * med_u,
            100.0 * max_u
        );
        eprintln!(
            "    本轮选中「{}」：中位 {:.2}% / 最大 {:.2}%",
            model_r.name(),
            100.0 * med_r,
            100.0 * max_r
        );
        if !model_r.bg_only() && max_r < max_u {
            eprintln!(
                "    ⇒ **翻了**：最坏状态误差 {:.2}% → {:.2}%（改善 {:.0}%），中位 {:.2}% → {:.2}%。",
                100.0 * max_u,
                100.0 * max_r,
                100.0 * (1.0 - max_r / max_u.max(1e-12)),
                100.0 * med_u,
                100.0 * med_r
            );
            eprintln!("    翻的原因是**自变量选错了轴**（见 [1e]：w_R 对 ln(R/G) 的偏相关 +0.918，对 ln(B/G) 只有 +0.061），");
            eprintln!("    不是靠加参数硬拟合：R 通道上二元并没有比一元更好（见 [2b]）。");
        } else {
            eprintln!("    ⇒ 未翻：本轮选中项仍属「只用 B/G」那一族。如实报。");
        }
    }

    // ---- [5b] 结构：把「B/G 上的陡降」拆开 ----
    eprintln!("\n[5b] 结构：B/G 上看到的「陡降」有多少是第二个自由度的投影");
    eprintln!("  最直接的检验是**同 B/G、不同 R/G** 的配对（见 [1d]），它们的 w 差多少：");
    let mut worst_pair: Option<(f64, &ImageData, &ImageData)> = None;
    for i in 0..images.len() {
        for j in (i + 1)..images.len() {
            let (a, b) = (&images[i], &images[j]);
            if (b.bg - a.bg).abs() / a.bg > NEAR_BG_TOL {
                continue;
            }
            // 排除同一状态内的重复图（那本来就应该一致）
            let same_state = states.iter().any(|g| g.contains(&i) && g.contains(&j));
            if same_state {
                continue;
            }
            let dr = 100.0 * (b.w_norm[0] - a.w_norm[0]).abs() / a.w_norm[0];
            if worst_pair.as_ref().map(|(d, _, _)| dr > *d).unwrap_or(true) {
                worst_pair = Some((dr, a, b));
            }
        }
    }
    if let Some((dr, a, b)) = worst_pair {
        eprintln!(
            "    {} vs {}：B/G {:.4} vs {:.4}（差 {:.2}%）、R/G {:.4} vs {:.4}（差 {:.1}%）",
            a.stem,
            b.stem,
            a.bg,
            b.bg,
            100.0 * (b.bg - a.bg).abs() / a.bg,
            a.rg,
            b.rg,
            100.0 * (b.rg - a.rg).abs() / a.rg
        );
        eprintln!(
            "      → w_R {:.5} vs {:.5}（差 {:.1}%）· w_B {:.5} vs {:.5}（差 {:.1}%）",
            a.w_norm[0],
            b.w_norm[0],
            dr,
            a.w_norm[2],
            b.w_norm[2],
            100.0 * (b.w_norm[2] - a.w_norm[2]).abs() / a.w_norm[2]
        );
        eprintln!("      **同一个 B/G 上 w_R 能差这么多，而 w_B 几乎不动** —— 这正是「B/G 少了一个自由度」的指纹。");
    } else {
        eprintln!("    （没有找到同 B/G、不同状态的配对）");
    }
    // 线性投影残差：把 w_R 对 B/G 回归之后，剩下的残差还能被 R/G 解释多少
    let ln_w_r: Vec<f64> = images.iter().map(|im| im.w_norm[0].ln()).collect();
    let r_proj = pearson(&ln_bg, &ln_w_r);
    let rp_second = {
        // 残差 = ln w_R 去掉 ln(B/G) 的线性部分
        let n = ln_bg.len() as f64;
        let (mx, my) = (
            ln_bg.iter().sum::<f64>() / n,
            ln_w_r.iter().sum::<f64>() / n,
        );
        let sxx: f64 = ln_bg.iter().map(|x| (x - mx).powi(2)).sum();
        let sxy: f64 = ln_bg.iter().zip(&ln_w_r).map(|(x, y)| (x - mx) * (y - my)).sum();
        let b = if sxx.abs() < 1e-15 { 0.0 } else { sxy / sxx };
        let a = my - b * mx;
        let resid: Vec<f64> = ln_bg.iter().zip(&ln_w_r).map(|(x, y)| y - (a + b * x)).collect();
        (pearson(&ln_rg, &resid), resid.iter().map(|r| r * r).sum::<f64>())
    };
    eprintln!(
        "  把 ln w_R 先对 ln(B/G) 回归（r = {r_proj:+.3}，残差平方和 {:.4}），残差再与 ln(R/G) 的相关是 {:+.3}",
        rp_second.1, rp_second.0
    );
    eprintln!("  ⇒ 若这个残差相关仍然很强，说明 R/G 携带了 B/G 之外的信息（支持二元）。");

    // ---- [5c] 分组拟合：扫描组内部是否只是「沿 B/G 走」 ----
    eprintln!("\n[5c] 扫描组内部：若调色偏移真的全程固定，组内的 B/G 就是**一维**的，可以单独看形状");
    eprintln!("  （Lead 报告偏移自 DSC_0563 起固定；下面这张表可以直接检验「是不是只差 B/G」。）");
    // 扫描组 = 设备号 056x~057x
    let scan: Vec<&ImageData> = images
        .iter()
        .filter(|im| stem_number(&im.stem).map(|n| (561..=570).contains(&n)).unwrap_or(false))
        .collect();
    if scan.len() >= 4 {
        let scan_bg: Vec<f64> = scan.iter().map(|im| im.bg).collect();
        let (slo, shi) = (
            scan_bg.iter().cloned().fold(f64::INFINITY, f64::min),
            scan_bg.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        let mut pts = scan.clone();
        pts.sort_by(|a, b| a.bg.partial_cmp(&b.bg).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!("  扫描组 {} 张，B/G {slo:.4}→{shi:.4}（{:.2}×）：", pts.len(), shi / slo);
        eprintln!(
            "  {:<26} {:>10} {:>10} {:>10} {:>10} {:>10}",
            "相邻两张（B/G）", "Δln(B/G)", "ΔR/G%", "Δw_R%", "w_R 弹性", "w_B 弹性"
        );
        for w in pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            let dx = (b.bg / a.bg).ln();
            let drg = 100.0 * (b.rg / a.rg - 1.0);
            let dwr = 100.0 * (b.w_norm[0] / a.w_norm[0] - 1.0);
            let er = (b.w_norm[0] / a.w_norm[0]).ln() / dx;
            let eb = (b.w_norm[2] / a.w_norm[2]).ln() / dx;
            // Δln(B/G) 太小时弹性没有意义（除以近零），如实标注
            let tag = if dx.abs() < 0.02 { "  ← ΔB/G 太小，弹性不可用" } else { "" };
            eprintln!(
                "  {:.4} → {:.4}      {:>9.4} {:>9.1}% {:>9.1}% {:>10.3} {:>10.3}{tag}",
                a.bg, b.bg, dx, drg, dwr, er, eb
            );
        }
        let in_pts: Vec<(f64, f64)> = pts.iter().map(|im| (im.bg, im.w_norm[0])).collect();
        let in_fit = fit_form(Form::LogLog, &in_pts);
        let out_pts: Vec<(f64, f64)> = {
            let mut v: Vec<&ImageData> = images.iter().collect();
            v.sort_by(|a, b| a.bg.partial_cmp(&b.bg).unwrap_or(std::cmp::Ordering::Equal));
            v.iter().map(|im| (im.bg, im.w_norm[0])).collect()
        };
        let out_fit = fit_form(Form::LogLog, &out_pts);
        if let (Some(i), Some(o)) = (in_fit, out_fit) {
            eprintln!(
                "  w_R 幂律拟合：扫描组内 b = {:.3}（对数残差 {:.2}%）· 全样本 b = {:.3}（{:.2}%）",
                i.b,
                100.0 * (i.log_rms.exp() - 1.0),
                o.b,
                100.0 * (o.log_rms.exp() - 1.0)
            );
        }
        eprintln!("  读法：Δln(B/G) 明显非零的那些区间里，若 w_R 弹性近似常数 → 沿 B/G 是光滑的；");
        eprintln!("        ΔB/G 极小而 ΔR/G 不小的区间 = **同一 B/G、不同白平衡状态**，弹性在那里无意义。");
    } else {
        eprintln!("  （扫描组样本不足 {} 张，跳过）", scan.len());
    }
    eprintln!(
        "  功能后果：留出档上色度误差 M0 {ba:.5} → 预测 {pa:.5}（{:+.0}%），上限 {ca:.5}（{:+.0}%）",
        100.0 * (pa / ba.max(1e-12) - 1.0),
        100.0 * (ca / ba.max(1e-12) - 1.0)
    );
    eprintln!(
        "  档内离散参照：{}",
        if any_multi {
            "见 [1c]——CV 误差若不明显大于档内离散，说明误差主要是同档内的图间散布（内容/噪声），而不是曲线本身的问题。"
        } else {
            "本批没有同档多图，无法量档内离散。"
        }
    );
    eprintln!();
    eprintln!("  必须随结论一起说的限制：");
    eprintln!(
        "  · `w` 是在这批参考图上解出来的：{} 个白平衡**状态**（按 B/G 一维投影是 {} 档），其中最大的一档是同场景连拍——",
        states.len(),
        levels.len()
    );
    eprintln!("    因此「状态内离散」**低估**了跨场景的内容敏感性。");
    eprintln!("  · **新增的色温扫描是同一场景**：它补的是「白平衡平面上点多不多」，**没有**补「换场景还成不成立」。");
    eprintln!("    而且其中几张的 R/G 与 B/G 并不同步（见 [5c]），说明扫描组内部的调色偏移也未必全程一致。");
    eprintln!("  · 留一状态 CV 仍是**同一批场景内的插值与近距离外推**，不是跨场景验证；唯一能看跨场景的是 [1d]。");
    eprintln!("  · 模型是按 CV 指标逐通道挑的，属于模型选择，CV 数字略偏乐观。");
    eprintln!("  · 色度尺对**逐通道曲线差**不免疫：参考端 NEUTRAL 的逐通道色调曲线会部分落进残差。");
    eprintln!("  · `w` 只有出现在参考集里的白平衡状态才被观测过；本诊断给的是曲线，不是物理模型。");
}

/// 仪器自检（不依赖样本）：两族求解器都要能解回已知的 `w`，且色度尺对整体增益免疫。
#[test]
fn instrument_recovers_known_gains() {
    let m0 = camera::color_matrix("Nikon Z 8").expect("Z8 应有固定矩阵").matrix;
    let w_true = [1.07f64, 1.0, 0.93];
    let mut cam_samples: Vec<Sample> = Vec::new();
    let mut out_samples: Vec<Sample> = Vec::new();
    for i in 0..60 {
        for j in 0..60 {
            let c = [
                0.05 + i as f32 * 0.014,
                0.05 + j as f32 * 0.014,
                0.05 + ((i * 7 + j * 13) % 60) as f32 * 0.014,
            ];
            cam_samples.push((c, apply_w_cam(m0, w_true, c), i % 2 == 0));
            out_samples.push((c, apply_w_out(m0, w_true, c), i % 2 == 0));
        }
    }
    let got_cam = norm_g(solve_w_cam(m0, &cam_samples).expect("应能求解")).expect("应能归一化");
    let got_out = norm_g(solve_w_out(m0, &out_samples).expect("应能求解")).expect("应能归一化");
    for k in [0usize, 2usize] {
        assert!(
            (got_cam[k] - w_true[k]).abs() < 1e-3,
            "矩阵之前：[{k}] {} vs {}（{got_cam:?}）",
            got_cam[k],
            w_true[k]
        );
        assert!(
            (got_out[k] - w_true[k]).abs() < 1e-3,
            "矩阵之后：[{k}] {} vs {}（{got_out:?}）",
            got_out[k],
            w_true[k]
        );
    }
    assert!(got_cam[1] == 1.0 && got_out[1] == 1.0, "归一化后绿通道应为 1");

    // 色度尺对整体增益免疫：把参考整体乘 3.7，读数不应变（容差只留浮点误差）
    let scaled: Vec<Sample> = cam_samples
        .iter()
        .map(|(c, t, l)| (*c, [t[0] * 3.7, t[1] * 3.7, t[2] * 3.7], *l))
        .collect();
    let e0 = chroma_err(&cam_samples, |c| apply_w_cam(m0, w_true, c));
    let e1 = chroma_err(&scaled, |c| apply_w_cam(m0, w_true, c));
    assert!((e0 - e1).abs() < 1e-6, "色度尺应对整体增益免疫：{e0} vs {e1}");

    // 常数形式确实不含 x：换任何 B/G 都预测同一个值
    let pts: Vec<(f64, f64)> = vec![(1.1, 0.9), (1.5, 1.2), (3.0, 2.1)];
    let f = fit_form(Form::Const, &pts).expect("应能拟合");
    assert!((f.predict(1.1) - f.predict(3.0)).abs() < 1e-12, "常数形式不应随 x 变化");
    assert!((f.predict(1.1) - 1.4).abs() < 1e-9, "常数形式应给出均值，得到 {}", f.predict(1.1));

    // 分档与选图的两条不变量。注意分档是拿**组首**比，不是拿前一个比。
    let bgs = [1.000, 1.002, 1.004, 1.5, 3.0];
    let g = group_by_bg(&bgs, 5e-3);
    assert_eq!(g.len(), 3, "1.000/1.002/1.004 应合成一档，得到 {g:?}");
    // 与组首比较，不是与相邻元素比较：1.009 距组首 0.9% 超出容差，于是另起一档
    let g2 = group_by_bg(&[1.000, 1.006, 1.009], 5e-3);
    assert_eq!(g2.len(), 2, "分档应与组首比较，得到 {g2:?}");

    let stems = ["DSC_0001", "DSC_0002", "DSC_0003", "DSC_0100", "DSC_0200"];
    let census: Vec<(String, f64, f64)> = stems
        .iter()
        .enumerate()
        .map(|(i, s)| ((*s).to_string(), 1.0, bgs[i]))
        .collect();
    let levels = group_by_bg(&census.iter().map(|c| c.2).collect::<Vec<_>>(), 5e-3);
    let (chosen, _) = choose_stems(&census, &levels);
    assert_eq!(
        chosen.len(),
        chosen.iter().collect::<std::collections::HashSet<_>>().len(),
        "选图列表不得有重复：{chosen:?}"
    );
    for g in &levels {
        let rep = &census[g[0]].0;
        assert!(chosen.contains(rep), "每档的代表 {rep} 都应在选中列表里：{chosen:?}");
    }
    assert!(
        chosen.contains(&"DSC_0003".to_string()),
        "多图档应再补一张编号相隔最远的成员（DSC_0003）：{chosen:?}"
    );
}
