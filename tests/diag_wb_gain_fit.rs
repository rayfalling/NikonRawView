//! 诊断 K：拟合 `w(B/G)` 并做**留一档交叉验证**。
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
    extras.retain(|s| !reps.contains(s));
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
    /// 留出对象的实测值（单张图，或该档多张图的均值）。
    measured: [f64; 3],
    /// 用其余档拟合出来的预测值。
    predicted: [f64; 3],
}

/// 在训练集上为 R、B 两个通道各拟合一个形式。
fn fit_channels(train: &[&ImageData], form_r: Form, form_b: Form) -> Option<(Fit, Fit)> {
    let fr = fit_form(form_r, &points(train, 0))?;
    let fb = fit_form(form_b, &points(train, 2))?;
    Some((fr, fb))
}

fn predict_vec(fr: &Fit, fb: &Fit, bg: f64) -> [f64; 3] {
    [fr.predict(bg), 1.0, fb.predict(bg)]
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

fn mean_bg(images: &[&ImageData]) -> f64 {
    if images.is_empty() {
        return f64::NAN;
    }
    images.iter().map(|im| im.bg).sum::<f64>() / images.len() as f64
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

/// 留一档：抽掉整整一档（多张图时全抽掉），用其余档拟合。
fn cv_rows_leave_level(
    images: &[ImageData],
    levels: &[Vec<usize>],
    form_r: Form,
    form_b: Form,
) -> Vec<CvRow> {
    let mut rows = Vec::new();
    for held in levels {
        let train = train_excluding(images, held);
        let Some((fr, fb)) = fit_channels(&train, form_r, form_b) else { continue };
        let members: Vec<&ImageData> = held.iter().map(|&i| &images[i]).collect();
        let bg = mean_bg(&members);
        rows.push(CvRow {
            label: if members.len() == 1 {
                members[0].stem.clone()
            } else {
                format!("{} 档（{} 张）", members[0].stem, members.len())
            },
            bg,
            measured: mean_w(&members),
            predicted: predict_vec(&fr, &fb, bg),
        });
    }
    rows
}

/// 留一图：每张图轮流留出。
fn cv_rows_leave_image(images: &[ImageData], form_r: Form, form_b: Form) -> Vec<CvRow> {
    let mut rows = Vec::new();
    for (i, held) in images.iter().enumerate() {
        let train = train_excluding(images, &[i]);
        let Some((fr, fb)) = fit_channels(&train, form_r, form_b) else { continue };
        rows.push(CvRow {
            label: held.stem.clone(),
            bg: held.bg,
            measured: held.w_norm,
            predicted: predict_vec(&fr, &fb, held.bg),
        });
    }
    rows
}

/// 单通道的留一档 CV 相对误差（逐通道选形式用）。
fn cv_channel_errors(
    images: &[ImageData],
    levels: &[Vec<usize>],
    form: Form,
    k: usize,
) -> Vec<f64> {
    let mut out = Vec::new();
    for held in levels {
        let train = train_excluding(images, held);
        let Some(f) = fit_form(form, &points(&train, k)) else { continue };
        let members: Vec<&ImageData> = held.iter().map(|&i| &images[i]).collect();
        let mw = members.iter().map(|im| im.w_norm[k]).sum::<f64>() / members.len() as f64;
        out.push(rel_err(f.predict(mean_bg(&members)), mw));
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
    let levels = group_by_bg(&bgs, WB_LEVEL_TOL);

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

    // ---- 档内离散 ----
    eprintln!("\n[1c] 档内离散（同档多张图时的 w 散布，是交叉验证误差的参照底）");
    let mut any_multi = false;
    for g in &levels {
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
    // 档内散度的总括（后面判读要用）：所有多图档的 w 相对跨度的最大值
    let mut scatter_r = f64::NAN;
    let mut scatter_b = f64::NAN;
    for g in &levels {
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
    let biggest = levels.iter().max_by_key(|g| g.len()).cloned().unwrap_or_default();
    eprintln!(
        "  {:<10} {:<10} {:>9} {:>9} {:>9} {:>9}  场景关系",
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
            let both_in_biggest = biggest.contains(&i) && biggest.contains(&j);
            if both_in_biggest {
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
                if both_in_biggest { "同档同场景（连拍）" } else { "同档，不同图/场景" }
            );
        }
    }
    if near_pairs == 0 {
        eprintln!("    （没有 B/G 足够接近的配对）");
    } else {
        eprintln!(
            "  ⇒ 最大差：同场景配对 {:.2}%（{} 对）· 跨图配对 {:.2}%（{} 对）",
            max_of(&same_scene),
            same_scene.len(),
            max_of(&cross_scene),
            cross_scene.len()
        );
        eprintln!("     跨图配对若与同场景配对同量级，说明 w 是白平衡档的属性、与场景无关。");
    }

    // ---- [2] 拟合 ----
    let all: Vec<&ImageData> = images.iter().collect();
    eprintln!("\n[2] 拟合 w(B/G)：{} 个形式 × 两个通道（系数与训练残差）", FORMS.len());
    eprintln!(
        "  {:<18} {:>4} {:>11} {:>11} {:>14}",
        "形式", "通道", "a", "b", "训练残差→%"
    );
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

    // ---- [2b] 逐通道选形式 ----
    eprintln!("\n[2b] 逐通道的留一档 CV 相对误差（中位 / 最大）——**每个通道各选各的形式**");
    eprintln!("  选择口径：**按最大误差（minimax）**挑。理由：要用的地方是「没见过的白平衡」，");
    eprintln!("  最坏的那一档才是风险；只看中位数会让「高原档很准、最暖档全错」的形式蒙混过关。");
    let mut chosen_form = [Form::Const; 3];
    for k in [0usize, 2usize] {
        let label = if k == 0 { "R" } else { "B" };
        let mut best_max: Option<(Form, f64, f64)> = None;
        let mut best_med: Option<(Form, f64, f64)> = None;
        for form in FORMS {
            let errs = cv_channel_errors(&images, &levels, form, k);
            let (med, mx) = (median_of(&errs), max_of(&errs));
            eprintln!(
                "  {:<22} {:>2}   中位 {:>7.2}%   最大 {:>7.2}%",
                form.name(),
                label,
                100.0 * med,
                100.0 * mx
            );
            if best_max.map(|(_, bm, _)| mx < bm).unwrap_or(true) {
                best_max = Some((form, mx, med));
            }
            if best_med.map(|(_, bm, _)| med < bm).unwrap_or(true) {
                best_med = Some((form, med, mx));
            }
        }
        if let Some((f, mx, med)) = best_max {
            chosen_form[k] = f;
            eprintln!("    ⇒ 通道 {label} 按**最大误差**选中：{}（中位 {:.2}% / 最大 {:.2}%）", f.name(), 100.0 * med, 100.0 * mx);
        }
        if let Some((f, med, mx)) = best_med {
            if f != chosen_form[k] {
                eprintln!(
                    "      若只按**中位**挑则是：{}（中位 {:.2}% / **最大 {:.2}%**）——中位与最坏档给出不同的形式，",
                    f.name(),
                    100.0 * med,
                    100.0 * mx
                );
                eprintln!("      这本身就是「中位数会掩盖最坏档」的实证。");
            }
        }
    }
    let (form_r, form_b) = (chosen_form[0], chosen_form[2]);
    if let (Some(fr), Some(fb)) =
        (fit_form(form_r, &points(&all, 0)), fit_form(form_b, &points(&all, 2)))
    {
        eprintln!(
            "  全量拟合：w_R = {}（a {:.5}, b {:.5}）· w_B = {}（a {:.5}, b {:.5}）",
            form_r.name(), fr.a, fr.b, form_b.name(), fb.a, fb.b
        );
    }
    eprintln!("  （形式是按 CV 指标选的，属于模型选择；下面的 CV 数字因此略偏乐观。）");
    if form_r == Form::Const || form_b == Form::Const {
        eprintln!("  注意：有通道选中了**常数**——那就等于说「该通道与 B/G 无关」，不是硬套曲线。");
    }

    // ---- [3] 留一档交叉验证 ----
    eprintln!("\n[3] 留一档交叉验证（重点）：抽掉整整一档，用其余档拟合后预测它");
    eprintln!("  误差按**档**统计（每档等权，与档内张数无关——最大的一档只算一次）。");
    eprintln!(
        "  {:<22} {:>8} {:>10} {:>10} {:>8} {:>10} {:>10} {:>8}",
        "留出档", "B/G", "实测w_R", "预测w_R", "误差", "实测w_B", "预测w_B", "误差"
    );
    let rows = cv_rows_leave_level(&images, &levels, form_r, form_b);
    let mut errs_r = Vec::new();
    let mut errs_b = Vec::new();
    for r in &rows {
        let er = rel_err(r.predicted[0], r.measured[0]);
        let eb = rel_err(r.predicted[2], r.measured[2]);
        errs_r.push(er);
        errs_b.push(eb);
        eprintln!(
            "  {:<22} {:>8.4} {:>10.5} {:>10.5} {:>7.2}% {:>10.5} {:>10.5} {:>7.2}%",
            r.label, r.bg, r.measured[0], r.predicted[0], 100.0 * er, r.measured[2], r.predicted[2], 100.0 * eb
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

    // ---- [3b] 档容差敏感性 ----
    eprintln!("\n[3b] 档容差敏感性：近邻档会让「留一档」变得几乎不用外推，CV 因此偏乐观");
    eprintln!("  主结果用 {:.1}% 容差；下面把容差放宽重算，看 CV 对「档怎么切」有多敏感。", WB_LEVEL_TOL * 100.0);
    eprintln!("  {:<8} {:>5} {:>24} {:>24}", "容差", "档数", "w_R 中位/最坏", "w_B 中位/最坏");
    for tol in std::iter::once(WB_LEVEL_TOL).chain(TOL_SENSITIVITY) {
        let lv = group_by_bg(&bgs, tol);
        let rows = cv_rows_leave_level(&images, &lv, form_r, form_b);
        let er: Vec<f64> = rows.iter().map(|r| rel_err(r.predicted[0], r.measured[0])).collect();
        let eb: Vec<f64> = rows.iter().map(|r| rel_err(r.predicted[2], r.measured[2])).collect();
        eprintln!(
            "  {:<8} {:>5} {:>23} {:>23}",
            format!("{:.1}%", tol * 100.0),
            lv.len(),
            format!("{:.2}% / {:.2}%", 100.0 * median_of(&er), 100.0 * max_of(&er)),
            format!("{:.2}% / {:.2}%", 100.0 * median_of(&eb), 100.0 * max_of(&eb))
        );
    }

    // ---- [3c] 功能后果 ----
    eprintln!("\n[3c] 功能后果（留出档上，尺度不变色度误差中位）：M0 原样 → 用预测 w → 用实测 w");
    eprintln!(
        "  {:<22} {:>10} {:>10} {:>12} {:>12}",
        "留出档", "M0", "预测 w", "实测档均值 w", "本图自解 w"
    );
    let mut before_all = Vec::new();
    let mut pred_all = Vec::new();
    let mut ceil_all = Vec::new();
    let mut self_all = Vec::new();
    for (held, r) in levels.iter().zip(&rows) {
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
        "  ⇒ 全体留出样本：M0 {ba:.5} → 预测 {pa:.5}（{:+.0}%）· 实测档均值（上限）{ca:.5}（{:+.0}%）· 本图自解 {sa:.5}",
        100.0 * (pa / ba.max(1e-12) - 1.0),
        100.0 * (ca / ba.max(1e-12) - 1.0)
    );

    // ---- [4] 留一图交叉验证 ----
    eprintln!("\n[4] 留一图交叉验证（每张图轮流留出；训练集里仍有同档图的那些行检验**档内泛化**）");
    let img_rows = cv_rows_leave_image(&images, form_r, form_b);
    eprintln!(
        "  {:<10} {:>8} {:>9} {:>9} {:>8} {:>9} {:>9} {:>8}  同档还有图？",
        "留出图", "B/G", "实测w_R", "预测w_R", "误差", "实测w_B", "预测w_B", "误差"
    );
    let mut same_level: Vec<f64> = Vec::new();
    let mut other_level: Vec<f64> = Vec::new();
    for (i, r) in img_rows.iter().enumerate() {
        let er = rel_err(r.predicted[0], r.measured[0]);
        let eb = rel_err(r.predicted[2], r.measured[2]);
        let has_same = levels.iter().find(|g| g.contains(&i)).map(|g| g.len() > 1).unwrap_or(false);
        if has_same {
            same_level.push(er.max(eb));
        } else {
            other_level.push(er.max(eb));
        }
        eprintln!(
            "  {:<10} {:>8.4} {:>9.5} {:>9.5} {:>7.2}% {:>9.5} {:>9.5} {:>7.2}%  {}",
            r.label,
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
        "\n  训练集里**仍有同档图**（{} 张）：中位 {:.2}% · 最大 {:.2}%",
        same_level.len(),
        100.0 * median_of(&same_level),
        100.0 * max_of(&same_level)
    );
    eprintln!(
        "  训练集里**没有同档图**（{} 张）：中位 {:.2}% · 最大 {:.2}%",
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
    for (label, form, med, mx, amp, ratio) in [
        ("R", form_r, med_r, max_r, amp_r, ratio_r),
        ("B", form_b, med_b, max_b, amp_b, ratio_b),
    ] {
        eprintln!(
            "  通道 {label}：选中「{}」· CV 中位 {:.2}% / **最大 {:.2}%** · 修正幅度 {:.2}% · 误差÷幅度 {:.2} · 档内噪声底 {:.2}%",
            form.name(),
            100.0 * med,
            100.0 * mx,
            100.0 * amp,
            ratio,
            if label == "R" { scatter_r } else { scatter_b }
        );
        if form == Form::Const || form == Form::Median {
            eprintln!("     ⇒ **该通道与 B/G 无关**：曲线跑不赢常数。");
        } else if med <= ERR_USABLE && ratio < ERR_VS_AMPLITUDE {
            eprintln!("     ⇒ 中位在门槛内且误差远小于修正幅度：**该通道的 w(B/G) 是一阶可用的模型**。");
        } else if ratio >= ERR_VS_AMPLITUDE {
            eprintln!("     ⇒ **负结论**：预测误差与要修正的量本身同量级，该通道不能只靠 B/G 曲线。");
        } else {
            eprintln!("     ⇒ 抓住了主要趋势，但精度介于两者之间：可作一阶修正，不能当精确标定。");
        }
    }

    // ---- [5b] 结构观察 ----
    eprintln!("\n[5b] 结构观察（比系数更重要）");
    // w_R 的「高原」：从最低 B/G 起，w_R 一直保持在首张的 1% 以内的最长前缀
    let w0 = images[0].w_norm[0];
    let mut plateau = 0usize;
    for (i, im) in images.iter().enumerate().skip(1) {
        if (im.w_norm[0] - w0).abs() / w0 <= 0.01 {
            plateau = i;
        } else {
            break;
        }
    }
    let pw: Vec<f64> = images[..=plateau].iter().map(|im| im.w_norm[0]).collect();
    let (plo, phi) = (
        pw.iter().cloned().fold(f64::INFINITY, f64::min),
        pw.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    eprintln!(
        "  w_R 在 B/G ≤ {:.4} 的 {} 张上几乎不动：{:.5} ~ {:.5}（跨 {:.2}%），与 {:.2}% 的修正幅度相比等于「不用修」；",
        images[plateau].bg,
        plateau + 1,
        plo,
        phi,
        100.0 * (phi - plo) / plo,
        100.0 * amp_r
    );
    eprintln!("  此后逐档下降：");
    for im in images.iter().skip(plateau + 1) {
        eprintln!("      B/G {:>8.4}  →  w_R {:>9.5}（比高原低 {:>5.1}%）", im.bg, im.w_norm[0], 100.0 * (1.0 - im.w_norm[0] / phi));
    }
    let exp_b = fit_form(form_b, &points(&all, 2)).map(|f| f.b).unwrap_or(f64::NAN);
    let (wb_lo, wb_hi) = (images[0].w_norm[2], images[images.len() - 1].w_norm[2]);
    eprintln!(
        "  w_B 则随 B/G 上升（{wb_lo:.5} → {wb_hi:.5}，×{:.2}），全量幂律拟合的指数 b = {exp_b:.3}；",
        wb_hi / wb_lo
    );
    eprintln!("  以上是原样读数；「陡降还是光滑过渡」由 [5c] 的局部弹性表判定，不在本节下结论。");

    // ---- [5c] 局部弹性：判定「陡降」还是「光滑过渡」 ----
    eprintln!("\n[5c] 局部弹性 d ln w / d ln(B/G)（幂律的局部指数）——直接回答「陡降还是光滑过渡」");
    eprintln!("  按**档均值**算相邻档之间（同档多张图先平均，避免除以近零的 Δln x）。");
    eprintln!(
        "  {:<22} {:>10} {:>10} {:>10} {:>10}",
        "相邻档（B/G）", "Δln(B/G)", "R 弹性", "B 弹性", "区间跨度"
    );
    let mut level_pts: Vec<(f64, f64, f64)> = Vec::new(); // (bg, w_R, w_B) 档均值
    for g in &levels {
        let n = g.len() as f64;
        let bg = g.iter().map(|&i| images[i].bg).sum::<f64>() / n;
        let wr = g.iter().map(|&i| images[i].w_norm[0]).sum::<f64>() / n;
        let wb = g.iter().map(|&i| images[i].w_norm[2]).sum::<f64>() / n;
        level_pts.push((bg, wr, wb));
    }
    let mut er_all: Vec<f64> = Vec::new();
    let mut eb_all: Vec<f64> = Vec::new();
    for w in level_pts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let dx = (b.0 / a.0).ln();
        if dx.abs() < 1e-9 {
            continue;
        }
        let er = (b.1 / a.1).ln() / dx;
        let eb = (b.2 / a.2).ln() / dx;
        er_all.push(er);
        eb_all.push(eb);
        eprintln!(
            "  {:.4} → {:.4}  {:>10.4} {:>10.3} {:>10.3} {:>9.0}%",
            a.0,
            b.0,
            dx,
            er,
            eb,
            100.0 * (b.0 / a.0 - 1.0)
        );
    }
    let (rmin, rmax) = (
        er_all.iter().cloned().fold(f64::INFINITY, f64::min),
        er_all.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let (bmin, bmax) = (
        eb_all.iter().cloned().fold(f64::INFINITY, f64::min),
        eb_all.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    eprintln!(
        "\n  R 弹性范围 {rmin:.3} ~ {rmax:.3}（跨度 {:.3}）· B 弹性范围 {bmin:.3} ~ {bmax:.3}（跨度 {:.3}）",
        rmax - rmin,
        bmax - bmin
    );
    eprintln!("  若某通道的弹性近似常数 → 单一幂律成立（光滑过渡）；若弹性在个别区间远大于其他区间 → 那里是「陡降」。");
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
        "  · 参考图场景各异，`w` 是在它们身上解出来的；白平衡有 {} 档，且最大的一档是同场景连拍——",
        levels.len()
    );
    eprintln!("    因此「档内离散」**低估**了跨场景的内容敏感性。");
    eprintln!("  · **新增的色温扫描是同一场景**：它补的是「B/G 轴上点多不多」，**没有**补「换场景还成不成立」。");
    eprintln!("    留一档 CV 仍是**相邻档之间的插值**，不是跨场景验证；唯一能看跨场景的是 [1d]。");
    eprintln!("  · 档内还可能有近邻档（B/G 只差零点几个百分点），那会让留一档 CV 偏乐观——见 [3b]。");
    eprintln!("  · 形式是按 CV 指标逐通道挑的，属于模型选择，CV 数字略偏乐观。");
    eprintln!("  · 色度尺对**逐通道曲线差**不免疫：参考端 NEUTRAL 的逐通道色调曲线会部分落进残差。");
    eprintln!("  · `w` 只有出现在参考集里的白平衡才被观测过；本诊断给的是曲线，不是物理模型。");
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
    for g in &levels {
        let rep = &census[g[0]].0;
        assert!(chosen.contains(rep), "每档的代表 {rep} 都应在选中列表里：{chosen:?}");
    }
    assert!(
        chosen.contains(&"DSC_0003".to_string()),
        "多图档应再补一张编号相隔最远的成员（DSC_0003）：{chosen:?}"
    );
}
