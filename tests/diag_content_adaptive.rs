//! 诊断 P：**修好一张图所需的修正，是否与"画面内容统计量"相关？**
//!
//! # 背景：为什么问这个
//!
//! Lead 实测：把"神谕逐通道曲线"换成用**标定集**拟合，跨图几乎不兑现
//! （ΔE00 均值 6.607 → 6.525，−1.2%）。也就是说，**修好一张图所需的修正，
//! 无法从别的图学出来**。只剩两种解释：
//!
//! - **A. 尼康的渲染有自适应成分** —— 那么"需要的修正"应当与画面的**内容**相关
//!   （亮度分布、动态范围、色度分布……）——本诊断的可检验形式；
//! - **B. 我们的解码有逐图误差** —— 那么"需要的修正"与任何内容统计量都**无关**
//!   （随机项或尚未识别的隐藏参数）。B 是零假设。
//!
//! # 协议：修正必须定义在"跨图基座"之上，而且要先量清"哪种参数化才有 headroom"
//!
//! 第一版把"需要的修正"定义成"从我们解码直接到参考"的曲线，结果各臂不可比：
//! 那个口径下"修正"要独自承担整条色调映射，一个 6 参数模型当然打不过标定集
//! 拟合出来的 129 点曲线。**所以本版把口径改成残差，并对三种参数化并列测量**：
//!
//! ```text
//! 跨图基座 T_i ：只在**标定集**（其余图的左半）上拟合的「共享曲线 + 17³ LUT」
//! 需要的修正 y ：把 T_i 的输出映射到参考所需的**残差修正**，三种表示各测一rip
//!                · 6 参数：编码域逐通道曲线在 0.5 处的偏离（mid）+（0.75−0.25）的斜率差
//!                · 3 参数：逐通道对角增益 log2 g
//!                · 完整：129 点逐通道曲线（只作神谕上界，用来量参数化的天花板）
//! ```
//!
//! 所有臂共用**同一个** `T_i`，差别只在给它加什么残差修正——这才是可控对照。
//! `T_i` 是用标定集拟合的，所以它对应 Lead 那份"跨图不兑现"的现状；"无修正"臂的
//! ΔE00 就是那个现状的水平。
//!
//! # 做法（四步，第 4 步是唯一裁判）
//!
//! 1. **把"需要的修正"压成标量。** 对每张图 `j`、每个留出图 `i`：在 `j` 的**左半**
//!    上把 `T_i(o)` 映射到参考，用 `fit_curves_per_channel` 拟 3 条逐通道曲线，再取
//!    编码域 `x = 0.5` 的偏离与 `0.75 − 0.25` 的斜率差（6 参数）与残差对角增益（3 个
//!    `log2 g`）。理由：诊断 N 把跨图惩罚的 96.6% 归到「输出端、逐通道的色调/尺度差」，
//!    `mid/slope` 是它的水平与斜率、增益是最省的自由度。
//!    **每张图报的 `y_i` 取 `T_i`（没见过这张图的基座）下的值。**
//! 2. **量内容统计量。** 17 个，全部**只从我们自己的解码算**（运行时可获得）——
//!    参考图上的统计量在生产里拿不到，用它做预测等于作弊。
//! 3. **相关矩阵。** 13 个修正标量 × 17 个统计量 = 221 个候选对。把已有的两个
//!    自变量 `cam_mul` 的 `R/G`、`B/G`（以及色温 K）放进同一张表作**对照**：
//!    只有某个内容统计量**明显强于**它们，才算找到第三个自变量。
//! 4. **留一图验证（十条臂）。** 关键是顺序：**先问这条路的 headroom 有多大**
//!    （各"神谕"臂 = 用留出图自己左半解出的修正），**再问内容预测吃到了多少**
//!    （内容臂 vs 同表示的"平均"臂）。若神谕都改善不了，那"内容能不能预测"这个问题
//!    本身就是空的——不能拿一个无用的量去要求预测。
//!
//!    | 臂 | 修正怎么来 | 它回答什么 |
//!    |---|---|---|
//!    | 无修正 | 不加 | 现状（跨图不兑现的水平） |
//!    | 平均(6参/增益) | 标定集各图 `y_j` 的均值 | "用标定集平均"这个基线 |
//!    | 内容(嵌套) | 每个参数在**标定集内**选 \|r\| 最大的统计量 + 一元 OLS | 诚实的可达增益 |
//!    | 内容(偷看) | 统计量用**全部图**（含留出图）选 | 乐观上界，量化选择偏差 |
//!    | 神谕(6参/增益/完整曲线) | 留出图自己的修正（用它左半解出） | 各参数化的 headroom 上界 |
//!
//! # 纪律
//!
//! **相关性 ≠ 可达增益。** 本项目已因"结论说得比证据满"返工三次（最近一次就是
//! 诊断 N 的神谕 96.6% vs 实际 −1.2%）。所以：报告候选对数量、给出 Bonferroni
//! 校正后的显著性门槛、并**以第 4 步的留一图结果为准**下判定。只做到第 3 步的结论
//! 必须标注"这是相关性，未验证可达性"。
//!
//! 另有两条不依赖相关系数的证据：
//!
//! - **同内容对照组**：数据里有两组"同内容"图——`DSC_0561~0570`（同一场景的色温扫描，
//!   白平衡从 2560K 到 10000K）与 `DSC_0001~0015`（42 秒内连拍）。它们的内容几乎相同，
//!   若"需要的修正"在组内就差很多，那它不可能由内容决定。
//! - **内容签名近邻**：用亮度分位 + 12×12 粗格归一化亮度算图与图之间的内容距离，
//!   不依赖事先知道哪几张是同一场景。
//!
//! `NRV_N_STEMS=DSC_0001,DSC_0141` 可只跑指定图（迭代用）。

mod common;

use common::samples_dir;
use nikonrawview::camera;
use nikonrawview::color;
use nikonrawview::deltae;
use nikonrawview::fit::{self, Sample};
use nikonrawview::icc;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::makernote;
use nikonrawview::mat3::{self, Mat3};
use nikonrawview::transform::BaseTransform;
use std::time::Instant;

// ---------------------------------------------------------------------------
// 口径（先定后看）
// ---------------------------------------------------------------------------

/// 默认图集：7 张不同场景 + 色温扫描组的 3 张（`DSC_0563/0567/0570`）。
///
/// 与诊断 N 的后续实验（`diag_per_channel_curve`）尽量重合，便于和"−1.2%"对齐；
/// 额外补 `DSC_0563` 是为了让**同内容组**有 3 个成员。
const DEFAULT_STEMS: [&str; 10] = [
    "DSC_0001", "DSC_0010", "DSC_0141", "DSC_0314", "DSC_0507", "DSC_0540", "DSC_0563",
    "DSC_0567", "DSC_0570", "DSC_8562",
];

/// 采样步长（与诊断 N 一致）。
const STEP: usize = 17;
/// 色调曲线箱数（与生产一致）。
const CURVE_BINS: usize = 128;
/// 生产用的 LUT 边长。
const LUT_EDGE: usize = 17;
/// 池化标定集与解残差修正时的抽稀步长（不影响评估）。
const POOL_STRIDE: usize = 2;

/// 取曲线偏离的编码域采样点：中点与两端。
const PROBE: [f32; 3] = [0.25, 0.5, 0.75];

/// 样本门槛：任一侧饱和/近黑都丢掉（与既有诊断一致）。
const CAM_LO: f32 = 1e-4;
const CAM_HI: f32 = 0.999;
const REF_LO: f32 = 1e-5;
const REF_HI: f32 = 0.999;
/// 单张图可用的最少样本数。
const MIN_SAMPLES: usize = 2000;

/// 内容签名里的粗格边长与分位数个数。
const GRID: usize = 12;
const NQ: usize = 17;

/// 内容统计量的个数与名字（顺序即 `ImageData::stats` 的顺序）。
const NSTAT: usize = 17;
const STAT_NAMES: [&str; NSTAT] = [
    "亮度均值(线性)",
    "亮度P50(编码)",
    "亮度P1(编码)",
    "亮度P5(编码)",
    "亮度P95(编码)",
    "亮度P99(编码)",
    "动态范围P99-P1",
    "暗部占比<0.05",
    "高光占比>0.9",
    "C*中位",
    "C*P95",
    "C*宽度P95-P50",
    "R/G中位",
    "B/G中位",
    "cam R/G",
    "cam B/G",
    "色温K",
];

/// 修正模型的 6 个参数名。
const PARAM_NAMES: [&str; 6] = ["中R", "中G", "中B", "斜R", "斜G", "斜B"];
/// 由 6 参数派生的 4 个物理量：整体电平、整体对比、红绿偏、蓝绿偏。
const DERIVED_NAMES: [&str; 4] = ["电平", "对比", "红绿", "蓝绿"];
/// 3 个残差对角增益（log2）。
const GAIN_NAMES: [&str; 3] = ["增益R", "增益G", "增益B"];
/// 修正标量总列数。
const NCOL: usize = 6 + 4 + 3;

/// "统计量与修正强相关"的门槛（任务给的 |r| > 0.8）。
const R_STRONG: f64 = 0.8;
/// "还值得看一眼"的门槛。
const R_WEAK: f64 = 0.5;
/// 内容预测要算"比平均基线更好"的最小改善（ΔE00 中位，绝对量）。
const MIN_GAIN: f64 = 0.05;

// ---------------------------------------------------------------------------
// 基础工具
// ---------------------------------------------------------------------------

fn mean_of(v: &[f64]) -> f64 {
    if v.is_empty() {
        f64::NAN
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

fn median_of(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[v.len() / 2]
}

/// 最近秩分位（与 `common::quantiles` 同一口径）。
fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let rank = (q * sorted.len() as f64).ceil().max(1.0) as usize;
    sorted[rank.saturating_sub(1).min(sorted.len() - 1)]
}

/// 取双方都有限的配对，其余丢掉（色温之类的量可能有缺项）。
fn finite_pairs(xs: &[f64], ys: &[f64]) -> Vec<(f64, f64)> {
    xs.iter()
        .zip(ys.iter())
        .filter(|(a, b)| a.is_finite() && b.is_finite())
        .map(|(a, b)| (*a, *b))
        .collect()
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let pairs = finite_pairs(xs, ys);
    let n = pairs.len();
    if n < 3 {
        return f64::NAN;
    }
    let mx = mean_of(&pairs.iter().map(|p| p.0).collect::<Vec<_>>());
    let my = mean_of(&pairs.iter().map(|p| p.1).collect::<Vec<_>>());
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (x, y) in &pairs {
        let dx = x - mx;
        let dy = y - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return f64::NAN;
    }
    sxy / (sxx.sqrt() * syy.sqrt())
}

/// 一元 OLS 预测：`ŷ = ȳ + r·(s_y/s_x)·(x − x̄)`。
///
/// 用这个写法而不是 `a + b·x`，是因为 `|r|` 小时它会**自动收缩到均值**。
/// 但收缩不足以挡住外推失控，所以外层用 [`ols_predict_clamped`]。
fn ols_predict(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let pairs = finite_pairs(xs, ys);
    if pairs.len() < 3 {
        return mean_of(&pairs.iter().map(|p| p.1).collect::<Vec<_>>());
    }
    let mx = mean_of(&pairs.iter().map(|p| p.0).collect::<Vec<_>>());
    let my = mean_of(&pairs.iter().map(|p| p.1).collect::<Vec<_>>());
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (a, b) in &pairs {
        let dx = a - mx;
        sxy += dx * (b - my);
        sxx += dx * dx;
    }
    if sxx <= 0.0 {
        return my;
    }
    my + (sxy / sxx) * (x - mx)
}

/// 一元 OLS 预测，**结果夹在标定集 `y` 的取值范围内**。
///
/// # 为什么必须夹
///
/// 不夹时实测会失控：`高光占比>0.9` 这类几乎恒为 0（跨度 ~1e-4）的统计量一旦被选中，
/// 斜率会随 x 的微小差异爆炸，`DSC_8562` 的 ΔE00 从 9.7 直接跳到 65——那不是"内容
/// 预测不行"，是数值失控，会把整条臂的成绩毁掉、也可能把结论带偏。夹住之后每个参数
/// 都落在标定集见过的范围内，比较才是公平的（这就是"不外推"这条常规纪律）。
fn ols_predict_clamped(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let pairs = finite_pairs(xs, ys);
    if pairs.len() < 3 {
        return mean_of(&pairs.iter().map(|p| p.1).collect::<Vec<_>>());
    }
    let lo = pairs.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
    let hi = pairs.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    ols_predict(xs, ys, x).clamp(lo, hi)
}

/// Pearson `|r|` 的双侧临界值表（`df = n − 2`，下标 0 对应 `df = 1`）。
///
/// 由 `r = t/√(t²+df)` 换算，`t` 取学生 t 分布 0.975 / 0.995 分位。
/// 写死在表里而不是现算：本文件不需要通用统计库，而把换算过程写错比查表更容易。
const R_CRIT_05: [f64; 30] = [
    0.997, 0.950, 0.878, 0.811, 0.755, 0.707, 0.666, 0.632, 0.602, 0.576, 0.553, 0.532, 0.514,
    0.497, 0.482, 0.468, 0.456, 0.444, 0.433, 0.423, 0.413, 0.404, 0.396, 0.388, 0.381, 0.374,
    0.367, 0.361, 0.355, 0.349,
];
const R_CRIT_01: [f64; 30] = [
    1.000, 0.990, 0.959, 0.917, 0.875, 0.834, 0.798, 0.765, 0.735, 0.708, 0.684, 0.661, 0.641,
    0.623, 0.606, 0.590, 0.575, 0.561, 0.549, 0.537, 0.526, 0.515, 0.505, 0.496, 0.487, 0.478,
    0.470, 0.463, 0.456, 0.449,
];

/// `df = n − 2` 时 `|r|` 的 5% / 1% 临界值；超出表范围时返回 `NaN`。
fn r_critical(df: usize) -> (f64, f64) {
    if df == 0 || df > 30 {
        return (f64::NAN, f64::NAN);
    }
    (R_CRIT_05[df - 1], R_CRIT_01[df - 1])
}

// ---------------------------------------------------------------------------
// 修正的表示：6 参数 + 3 增益
// ---------------------------------------------------------------------------

/// 一张图"需要的残差修正"：编码域逐通道曲线相对恒等的偏离。
///
/// `mid_c` 是 `d_c(0.5)`，`slope_c` 是 `d_c(0.75) − d_c(0.25)`。
/// 施加时展开成定义在 `x = 0 / 0.25 / 0.5 / 0.75 / 1` 上的 5 点曲线
/// （两端固定为恒等），与神谕曲线在同一套编码域里。
#[derive(Debug, Clone, Copy, PartialEq)]
struct Corr {
    mid: [f64; 3],
    slope: [f64; 3],
}

impl Corr {
    fn params(&self) -> [f64; 6] {
        [
            self.mid[0], self.mid[1], self.mid[2], self.slope[0], self.slope[1], self.slope[2],
        ]
    }

    fn from_params(p: [f64; 6]) -> Self {
        Self { mid: [p[0], p[1], p[2]], slope: [p[3], p[4], p[5]] }
    }

    /// 派生的 4 个物理量：整体电平、整体对比、红绿偏、蓝绿偏。
    fn derived(&self) -> [f64; 4] {
        let level = (self.mid[0] + self.mid[1] + self.mid[2]) / 3.0;
        let contrast = (self.slope[0] + self.slope[1] + self.slope[2]) / 3.0;
        [level, contrast, self.mid[0] - self.mid[1], self.mid[2] - self.mid[1]]
    }

    /// 展开成 3 条 5 点曲线。
    fn curves(&self) -> [Vec<f32>; 3] {
        let mut out: [Vec<f32>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for (c, slot) in out.iter_mut().enumerate() {
            let d = |x: f32| self.mid[c] as f32 + self.slope[c] as f32 * (2.0 * x - 1.0);
            *slot = vec![0.0, 0.25 + d(PROBE[0]), 0.5 + d(PROBE[1]), 0.75 + d(PROBE[2]), 1.0];
        }
        out
    }

    /// 与另一个修正的逐参数平均绝对差。
    fn dist(&self, other: &Corr) -> f64 {
        let a = self.params();
        let b = other.params();
        (0..6).map(|k| (a[k] - b[k]).abs()).sum::<f64>() / 6.0
    }

    /// 逐参数均值。
    fn mean_of(list: &[Corr]) -> Corr {
        let mut acc = [0.0f64; 6];
        for c in list {
            let p = c.params();
            for k in 0..6 {
                acc[k] += p[k];
            }
        }
        for v in &mut acc {
            *v /= list.len().max(1) as f64;
        }
        Corr::from_params(acc)
    }
}

/// 由「曲线前的值 → 参考」的配对解出 6 参数修正。
fn extract_corr(samples: &[Sample]) -> Corr {
    let curves = fit::fit_curves_per_channel(samples, CURVE_BINS);
    let mut mid = [0.0f64; 3];
    let mut slope = [0.0f64; 3];
    for c in 0..3 {
        let d = |x: f32| fit::apply_curve(&curves[c], x) - x;
        mid[c] = d(PROBE[1]) as f64;
        slope[c] = (d(PROBE[2]) - d(PROBE[0])) as f64;
    }
    Corr { mid, slope }
}

/// 对角增益 `g_c = Σ(o·t)/Σ(o²)` 的 `log2`（正数表示我们的值比参考暗）。
fn residual_gain(samples: &[Sample]) -> [f64; 3] {
    let mut num = [0.0f64; 3];
    let mut den = [0.0f64; 3];
    for s in samples {
        for c in 0..3 {
            let o = s.ours[c] as f64;
            num[c] += o * s.theirs[c] as f64;
            den[c] += o * o;
        }
    }
    let mut out = [0.0f64; 3];
    for c in 0..3 {
        out[c] = if den[c] > 0.0 { (num[c] / den[c]).log2() } else { 0.0 };
    }
    out
}

// ---------------------------------------------------------------------------
// 逐图数据
// ---------------------------------------------------------------------------

/// 一个样本：当前管线的输出（`M0·diag(w)·c`，w 已施加）与参考值。
#[derive(Clone, Copy)]
struct Pair {
    o: [f32; 3],
    t: [f32; 3],
}

struct ImageData {
    stem: String,
    /// 左半（解修正用）与右半（评估用）。
    left: Vec<Pair>,
    right: Vec<Pair>,
    /// 17 个内容统计量。
    stats: [f64; NSTAT],
    /// 内容签名（分位 + 粗格），用于找"同内容"的图。
    sig: Vec<f64>,
    secs: f64,
}

/// 文档记载的"同内容"组（来自诊断历史，不是本次拟合出来的）。
fn documented_group(stem: &str) -> Option<&'static str> {
    let num = stem.strip_prefix("DSC_")?.parse::<u32>().ok()?;
    if (561..=570).contains(&num) {
        Some("色温扫描")
    } else if (1..=15).contains(&num) {
        Some("42秒连拍")
    } else {
        None
    }
}

/// 解码一张图并收集样本、内容统计量与内容签名。
fn collect(stem: &str, m0: Mat3) -> Option<ImageData> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        common::skip("诊断 P", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }
    let t0 = Instant::now();
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
    // 相机空间解码（一次），再自己施加 M0 与 w——与现行管线的位置一致。
    let cam = match libraw::decode_with_output_for_test(&nef, &opts, OutputColor::Camera) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("跳过 {stem}：相机空间解码失败（{e}）");
            return None;
        }
    };
    let model = match camera::read_model(&nef) {
        Some(m) => m,
        None => {
            eprintln!("跳过 {stem}：读不出机型");
            return None;
        }
    };
    let entry = match camera::lookup(&model) {
        Some(e) => e,
        None => {
            eprintln!("跳过 {stem}：机型 {model} 不在边距表中");
            return None;
        }
    };
    let mg = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = match mg.effective(cam.width, cam.height) {
        Some(v) => v,
        None => {
            eprintln!("跳过 {stem}：边距与解码尺寸不自洽");
            return None;
        }
    };
    if (cw, ch) != (img.width, img.height) {
        eprintln!("跳过 {stem}：裁切后 {cw}×{ch} 与参考 {}×{} 不符", img.width, img.height);
        return None;
    }

    // 白平衡状态（内容统计量与 w 都要用它）
    let wb = libraw::read_wb(&nef).ok()?;
    let cam_rg = (wb[0] / wb[1]) as f64;
    let cam_bg = (wb[2] / wb[1]) as f64;
    let w = color::wb_adaptation(cam_rg, cam_bg);

    let mut left: Vec<Pair> = Vec::new();
    let mut right: Vec<Pair> = Vec::new();
    let mut lums: Vec<f64> = Vec::new();
    let mut chromas: Vec<f64> = Vec::new();
    let mut rgs: Vec<f64> = Vec::new();
    let mut bgs: Vec<f64> = Vec::new();
    let mut lum_sum = 0.0f64;
    let mut grid = [0.0f64; GRID * GRID];
    let mut gcnt = [0u32; GRID * GRID];

    for y in (0..ch).step_by(STEP) {
        for x in (0..cw).step_by(STEP) {
            let Some(px) = cam.at(x + mg.left, y + mg.top) else { continue };
            let Some(rp) = img.at(x, y) else { continue };
            let c = [
                px[0] as f32 / 65535.0,
                px[1] as f32 / 65535.0,
                px[2] as f32 / 65535.0,
            ];
            let inv = 1.0 / 65535.0;
            let t = plan.to_working_space([
                rp[0] as f32 * inv,
                rp[1] as f32 * inv,
                rp[2] as f32 * inv,
            ]);
            if c.iter().any(|v| *v <= CAM_LO || *v >= CAM_HI)
                || t.iter().any(|v| *v <= REF_LO || *v >= REF_HI)
            {
                continue;
            }
            let o = mat3::mul_vec(
                m0,
                [c[0] * w[0] as f32, c[1] * w[1] as f32, c[2] * w[2] as f32],
            );
            let row = Pair { o, t };
            if x < cw / 2 {
                left.push(row);
            } else {
                right.push(row);
            }

            // ---- 内容统计量（只从我们自己的解码算）----
            let lum = (0.2126 * o[0] + 0.7152 * o[1] + 0.0722 * o[2]).max(0.0);
            lum_sum += lum as f64;
            lums.push(fit::encode_srgb(lum) as f64);
            let lab = deltae::lab_from_prophoto_linear(o);
            chromas.push(((lab[1] as f64).powi(2) + (lab[2] as f64).powi(2)).sqrt());
            if o.iter().all(|v| *v > 1e-4) {
                rgs.push((o[0] / o[1]) as f64);
                bgs.push((o[2] / o[1]) as f64);
            }
            let gx = (x * GRID / cw).min(GRID - 1);
            let gy = (y * GRID / ch).min(GRID - 1);
            grid[gy * GRID + gx] += lum as f64;
            gcnt[gy * GRID + gx] += 1;
        }
    }
    if left.len() < MIN_SAMPLES || right.len() < MIN_SAMPLES {
        eprintln!("跳过 {stem}：左右半样本不足（{} / {}）", left.len(), right.len());
        return None;
    }

    // ---- 统计量 ----
    lums.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    chromas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    rgs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    bgs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n_lum = lums.len().max(1) as f64;
    let p50 = quantile(&lums, 0.5);
    let p1 = quantile(&lums, 0.01);
    let p99 = quantile(&lums, 0.99);
    let lum_mean_lin = lum_sum / n_lum;
    let shadow = lums.iter().filter(|v| **v < 0.05).count() as f64 / n_lum;
    let highlight = lums.iter().filter(|v| **v > 0.9).count() as f64 / n_lum;
    let c_med = quantile(&chromas, 0.5);
    let c_p95 = quantile(&chromas, 0.95);
    let temp_k = makernote::color_temperature(&std::fs::read(&nef).unwrap_or_default())
        .ok()
        .flatten()
        .map(|k| k as f64)
        .unwrap_or(f64::NAN);
    let stats = [
        lum_mean_lin,
        p50,
        p1,
        quantile(&lums, 0.05),
        quantile(&lums, 0.95),
        p99,
        p99 - p1,
        shadow,
        highlight,
        c_med,
        c_p95,
        c_p95 - c_med,
        quantile(&rgs, 0.5),
        quantile(&bgs, 0.5),
        cam_rg,
        cam_bg,
        temp_k,
    ];

    // ---- 内容签名：17 个亮度分位（除以中位）+ 12×12 粗格（除以均值）----
    let mut sig: Vec<f64> = Vec::with_capacity(NQ + GRID * GRID);
    let denom = if p50 > 1e-9 { p50 } else { 1.0 };
    for i in 0..NQ {
        let q = (i as f64 + 1.0) / (NQ as f64 + 1.0);
        sig.push(quantile(&lums, q) / denom);
    }
    let gmean = if lum_mean_lin > 1e-9 { lum_mean_lin } else { 1.0 };
    for i in 0..GRID * GRID {
        let cell = if gcnt[i] > 0 { grid[i] / gcnt[i] as f64 } else { 0.0 };
        sig.push(cell / gmean);
    }

    let im = ImageData {
        stem: stem.to_string(),
        left,
        right,
        stats,
        sig,
        secs: t0.elapsed().as_secs_f64(),
    };
    eprintln!(
        "  采集 {:<10} 左半 {:>6} / 右半 {:>6} 样本，{:.1} s",
        im.stem,
        im.left.len(),
        im.right.len(),
        im.secs
    );
    Some(im)
}

/// 两张图的内容距离：签名的均方根差（分位已除以中位、粗格已除以均值，所以
/// 对整体曝光不敏感，反映的是**画面结构**）。
fn content_distance(a: &ImageData, b: &ImageData) -> f64 {
    let n = a.sig.len().min(b.sig.len());
    if n == 0 {
        return f64::NAN;
    }
    let s: f64 = (0..n).map(|i| (a.sig[i] - b.sig[i]).powi(2)).sum();
    (s / n as f64).sqrt()
}

// ---------------------------------------------------------------------------
// 跨图基座与残差修正
// ---------------------------------------------------------------------------

/// 只用标定集拟合的"跨图基座"：共享曲线 + 17³ LUT（现行做法的形状）。
fn build_base(pooled: &[Sample]) -> BaseTransform {
    let curve = fit::fit_curve(pooled, CURVE_BINS);
    let lut = fit::fit_lut(pooled, &curve, LUT_EDGE);
    BaseTransform::from_fit(0x0000, "跨图基座", Vec::new(), curve, LUT_EDGE, lut)
}

/// 用 `left` 半样本解出「把 `base` 的输出映射到参考」所需的残差修正与残差增益。
fn residual_of(left: &[Pair], base: &BaseTransform, stride: usize) -> (Corr, [f64; 3]) {
    let samples: Vec<Sample> = left
        .iter()
        .step_by(stride)
        .map(|r| Sample { ours: base.apply(r.o), theirs: r.t })
        .collect();
    (extract_corr(&samples), residual_gain(&samples))
}

/// 同样口径下的**完整逐通道曲线**（129 点），用来量"参数化到底吃掉多少"。
fn residual_full_curves(left: &[Pair], base: &BaseTransform, stride: usize) -> [Vec<f32>; 3] {
    let samples: Vec<Sample> = left
        .iter()
        .step_by(stride)
        .map(|r| Sample { ours: base.apply(r.o), theirs: r.t })
        .collect();
    fit::fit_curves_per_channel(&samples, CURVE_BINS)
}

/// 增益预测与神谕增益的平均绝对差（单位 log2）。
fn gain_err(a: &[f64], b: &[f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).abs()).sum::<f64>() / 3.0
}

/// 在 `sel` 上为 `ys` 挑 |r| 最大的统计量，返回 `(统计量下标, |r|)`。
fn pick_stat(ys: &[f64], stats: &[[f64; NSTAT]], sel: &[usize]) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64)> = None;
    let yy: Vec<f64> = sel.iter().map(|&i| ys[i]).collect();
    for (k, _name) in STAT_NAMES.iter().enumerate() {
        let xs: Vec<f64> = sel.iter().map(|&i| stats[i][k]).collect();
        let r = pearson(&xs, &yy).abs();
        if r.is_finite() && best.map(|b| r > b.1).unwrap_or(true) {
            best = Some((k, r));
        }
    }
    best
}

/// 内容预测（通用维度）：对每个参数各选一个统计量（`select_set` 为选择集）。
fn predict_vec(
    images: &[ImageData],
    ys: &[Vec<f64>],
    fit_set: &[usize],
    select_set: &[usize],
    target: usize,
) -> (Vec<f64>, Vec<usize>) {
    let stats: Vec<[f64; NSTAT]> = images.iter().map(|im| im.stats).collect();
    let dim = ys.first().map(|v| v.len()).unwrap_or(0);
    let mut out = vec![0.0f64; dim];
    let mut picks = Vec::with_capacity(dim);
    for k in 0..dim {
        let yk: Vec<f64> = (0..images.len()).map(|i| ys[i][k]).collect();
        let Some((stat, _r)) = pick_stat(&yk, &stats, select_set) else {
            out[k] = mean_of(&fit_set.iter().map(|&i| ys[i][k]).collect::<Vec<_>>());
            picks.push(0);
            continue;
        };
        let xs: Vec<f64> = fit_set.iter().map(|&i| images[i].stats[stat]).collect();
        let yy: Vec<f64> = fit_set.iter().map(|&i| ys[i][k]).collect();
        out[k] = ols_predict_clamped(&xs, &yy, images[target].stats[stat]);
        picks.push(stat);
    }
    (out, picks)
}

// ---------------------------------------------------------------------------
// 打印工具
// ---------------------------------------------------------------------------

fn num(v: f64, prec: usize) -> String {
    if v.is_finite() {
        format!("{v:.prec$}")
    } else {
        String::from("-")
    }
}

/// 按**显示宽度**补空格：中日韩字符占两列。
fn is_wide(c: char) -> bool {
    matches!(c as u32,
        0x1100..=0x115F
            | 0x2460..=0x24FF
            | 0x2E80..=0x303E
            | 0x3041..=0x33FF
            | 0x3400..=0x4DBF
            | 0x4E00..=0x9FFF
            | 0xA000..=0xA4CF
            | 0xAC00..=0xD7A3
            | 0xF900..=0xFAFF
            | 0xFE30..=0xFE6F
            | 0xFF00..=0xFF60
            | 0xFFE0..=0xFFE6)
}

fn pad(s: &str, width: usize) -> String {
    let w: usize = s.chars().map(|c| if is_wide(c) { 2 } else { 1 }).sum();
    format!("{s}{}", " ".repeat(width.saturating_sub(w)))
}

// ---------------------------------------------------------------------------
// 主测试
// ---------------------------------------------------------------------------

#[test]
fn needed_correction_versus_image_content() {
    let Some(entry) = camera::color_matrix("Nikon Z 8") else {
        eprintln!("跳过：机型表里没有 Z 8 的矩阵");
        return;
    };
    let m0 = entry.matrix;

    let stems: Vec<String> = match std::env::var("NRV_N_STEMS") {
        Ok(s) if !s.trim().is_empty() => s.split(',').map(|x| x.trim().to_string()).collect(),
        _ => DEFAULT_STEMS.iter().map(|s| s.to_string()).collect(),
    };

    eprintln!("=== 诊断 P：需要的修正 vs 画面内容（{} 张图）===", stems.len());
    let mut images: Vec<ImageData> = Vec::new();
    for stem in &stems {
        if let Some(im) = collect(stem, m0) {
            images.push(im);
        }
    }
    let n = images.len();
    if n < 6 {
        eprintln!("跳过：可用图不足（{n} 张，至少 6 张）");
        return;
    }

    // ---------- 留一图基座 + 每张图"需要的修正" ----------
    // 对每个留出图 i：基座 T_i 只在其余图上拟合；corr[i][j] 是图 j 相对 T_i 的残差修正。
    // 报出去的 y_i 取 corr[i][i]——基座没见过这张图，这才是"它需要的修正"。
    eprintln!("\n拟合 {} 个跨图基座（每个留出一张图）…", n);
    let tb = Instant::now();
    let mut bases: Vec<BaseTransform> = Vec::with_capacity(n);
    let mut corr: Vec<Vec<Corr>> = Vec::with_capacity(n);
    let mut gain: Vec<Vec<[f64; 3]>> = Vec::with_capacity(n);
    // 每个留出图自己的**完整**残差曲线（只用它左半解出，用来量参数化的上限）。
    let mut full: Vec<[Vec<f32>; 3]> = Vec::with_capacity(n);
    for hole in 0..n {
        let mut pooled: Vec<Sample> = Vec::new();
        for (j, im) in images.iter().enumerate() {
            if j == hole {
                continue;
            }
            for r in im.left.iter().step_by(POOL_STRIDE) {
                pooled.push(Sample { ours: r.o, theirs: r.t });
            }
        }
        let base = build_base(&pooled);
        let mut cs = Vec::with_capacity(n);
        let mut gs = Vec::with_capacity(n);
        for im in &images {
            let (c, g) = residual_of(&im.left, &base, POOL_STRIDE);
            cs.push(c);
            gs.push(g);
        }
        full.push(residual_full_curves(&images[hole].left, &base, POOL_STRIDE));
        bases.push(base);
        corr.push(cs);
        gain.push(gs);
    }
    eprintln!("基座与残差修正用了 {:.1} s", tb.elapsed().as_secs_f64());

    // 每张图"需要的修正"：基座没见过它的那一个
    let y: Vec<Corr> = (0..n).map(|i| corr[i][i]).collect();
    let ygain: Vec<[f64; 3]> = (0..n).map(|i| gain[i][i]).collect();
    let ys: Vec<[f64; 6]> = y.iter().map(|c| c.params()).collect();
    let derived: Vec<[f64; 4]> = y.iter().map(|c| c.derived()).collect();

    // ---------- 一、逐图：需要的修正 ----------
    eprintln!("\n【一、每张图「需要的修正」（残差口径：基座由其余图拟合，没见过本图）】");
    let head = [
        pad("图", 11),
        pad("camR/G", 8),
        pad("camB/G", 8),
        pad("色温K", 7),
        pad("gR", 7),
        pad("gG", 7),
        pad("gB", 7),
        pad("中R", 8),
        pad("中G", 8),
        pad("中B", 8),
        pad("斜R", 8),
        pad("斜G", 8),
        pad("斜B", 8),
        pad("组", 10),
    ];
    eprintln!("{}", head.join(" "));
    for (i, im) in images.iter().enumerate() {
        let cells = [
            pad(&im.stem, 11),
            pad(&num(im.stats[14], 3), 8),
            pad(&num(im.stats[15], 3), 8),
            pad(&num(im.stats[16], 0), 7),
            pad(&num(ygain[i][0], 3), 7),
            pad(&num(ygain[i][1], 3), 7),
            pad(&num(ygain[i][2], 3), 7),
            pad(&num(y[i].mid[0], 4), 8),
            pad(&num(y[i].mid[1], 4), 8),
            pad(&num(y[i].mid[2], 4), 8),
            pad(&num(y[i].slope[0], 4), 8),
            pad(&num(y[i].slope[1], 4), 8),
            pad(&num(y[i].slope[2], 4), 8),
            pad(documented_group(&im.stem).unwrap_or("—"), 10),
        ];
        eprintln!("{}", cells.join(" "));
    }
    eprintln!("  口径：中X = 该通道修正曲线在编码域 0.5 处相对恒等的偏离；斜X = 0.75 处减 0.25 处；");
    eprintln!("        gX = 残差对角增益的 log2（正数表示基座输出比参考暗）。组 = 文档记载的同内容组。");

    eprintln!("\n【内容统计量（只从我们自己的解码算，运行时可获得）】");
    let mut head2 = vec![pad("统计量", 17)];
    for im in &images {
        head2.push(pad(im.stem.trim_start_matches("DSC_"), 8));
    }
    eprintln!("{}", head2.join(" "));
    for (k, name) in STAT_NAMES.iter().enumerate() {
        let mut cells = vec![pad(name, 17)];
        for im in &images {
            cells.push(pad(&num(im.stats[k], 4), 8));
        }
        eprintln!("{}", cells.join(" "));
    }

    // ---------- 二、相关矩阵 ----------
    // 行 = 17 个统计量，列 = 13 个修正标量（6 参数 + 4 派生量 + 3 增益）。
    struct Column {
        name: &'static str,
        values: Vec<f64>,
    }
    let mut columns: Vec<Column> = Vec::new();
    for (k, name) in PARAM_NAMES.iter().enumerate() {
        columns.push(Column { name, values: ys.iter().map(|p| p[k]).collect() });
    }
    for (k, name) in DERIVED_NAMES.iter().enumerate() {
        columns.push(Column { name, values: derived.iter().map(|p| p[k]).collect() });
    }
    for (k, name) in GAIN_NAMES.iter().enumerate() {
        columns.push(Column { name, values: ygain.iter().map(|p| p[k]).collect() });
    }
    debug_assert_eq!(columns.len(), NCOL);

    eprintln!("\n【二、相关矩阵（Pearson r；* = |r|≥0.8，+ = |r|≥0.5）】");
    let mut head3 = vec![pad("统计量", 17)];
    for c in &columns {
        head3.push(pad(c.name, 7));
    }
    eprintln!("{}", head3.join(" "));

    let mut all_pairs: Vec<(String, &'static str, f64)> = Vec::new();
    for (k, stat_name) in STAT_NAMES.iter().enumerate() {
        let xs: Vec<f64> = images.iter().map(|im| im.stats[k]).collect();
        let mut cells = vec![pad(stat_name, 17)];
        for c in &columns {
            let r = pearson(&xs, &c.values);
            all_pairs.push((stat_name.to_string(), c.name, r));
            let mark = if !r.is_finite() {
                " "
            } else if r.abs() >= R_STRONG {
                "*"
            } else if r.abs() >= R_WEAK {
                "+"
            } else {
                " "
            };
            cells.push(pad(&format!("{}{mark}", num(r, 2)), 7));
        }
        eprintln!("{}", cells.join(" "));
    }

    let df = n - 2;
    let (rc05, rc01) = r_critical(df);
    let n_pairs = NSTAT * NCOL;
    let bonf_alpha = 0.05 / n_pairs as f64;
    // t(1 − α/2, df) ≈ 5.6（df=8，α≈2.3e-4）；r = t/√(t²+df)。只用于给出量级。
    let rc_bonf = 5.6 / (5.6f64 * 5.6 + df as f64).sqrt();
    eprintln!("\n  候选对数量：{NCOL} 个修正标量 × {NSTAT} 个统计量 = **{n_pairs} 对**。");
    eprintln!("  n={n} ⇒ df={df}：|r| ≥ {rc05:.3} 为 p<0.05，≥ {rc01:.3} 为 p<0.01。");
    eprintln!("  但 {n_pairs} 对一起看必须做多重比较校正：Bonferroni 后单对 α={bonf_alpha:.2e}，");
    eprintln!("  对应 |r| ≳ {rc_bonf:.2}（估算）。**按 0.8 门槛挑出来的'强相关'过不了这条线。**");

    let best = all_pairs
        .iter()
        .filter(|p| p.2.is_finite())
        .max_by(|a, b| a.2.abs().partial_cmp(&b.2.abs()).unwrap_or(std::cmp::Ordering::Equal));
    let controls = ["cam R/G", "cam B/G", "色温K"];
    let best_control = all_pairs
        .iter()
        .filter(|p| controls.contains(&p.0.as_str()) && p.2.is_finite())
        .max_by(|a, b| a.2.abs().partial_cmp(&b.2.abs()).unwrap_or(std::cmp::Ordering::Equal));
    let best_content = all_pairs
        .iter()
        .filter(|p| !controls.contains(&p.0.as_str()) && p.2.is_finite())
        .max_by(|a, b| a.2.abs().partial_cmp(&b.2.abs()).unwrap_or(std::cmp::Ordering::Equal));
    if let Some(b) = best {
        eprintln!("  全表最强（含对照量）：{} × {} ⇒ r = {:+.3}", b.0, b.1, b.2);
    }
    if let Some(b) = best_content {
        eprintln!("  内容量里最强：{} × {} ⇒ r = {:+.3}", b.0, b.1, b.2);
    }
    if let Some(c) = best_control {
        eprintln!("  对照（已有的两个自变量与色温）最强：{} × {} ⇒ r = {:+.3}", c.0, c.1, c.2);
    }
    eprintln!(
        "  候选量共 **{NSTAT} 个统计量**（{} 个内容量 + 3 个对照量）× {NCOL} 个修正标量。",
        NSTAT - 3
    );

    // ---------- 三、同内容对照 ----------
    eprintln!("\n【三、同内容对照：内容几乎相同，需要的修正是否也相同？】");
    let mut groups: Vec<(&'static str, Vec<usize>)> = Vec::new();
    for (i, im) in images.iter().enumerate() {
        if let Some(g) = documented_group(&im.stem) {
            match groups.iter_mut().find(|(name, _)| *name == g) {
                Some((_, v)) => v.push(i),
                None => groups.push((g, vec![i])),
            }
        }
    }
    let sd_of = |idx: &[usize]| -> Vec<f64> {
        (0..6)
            .map(|k| {
                let v: Vec<f64> = idx.iter().map(|&i| ys[i][k]).collect();
                let m = mean_of(&v);
                (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
            })
            .collect()
    };
    let all_idx: Vec<usize> = (0..n).collect();
    let overall_sd = sd_of(&all_idx);
    eprintln!(
        "  全体（n={n}）6 参数标准差：{}",
        (0..6)
            .map(|k| format!("{} {:.4}", PARAM_NAMES[k], overall_sd[k]))
            .collect::<Vec<_>>()
            .join("  ")
    );
    for (name, idx) in &groups {
        if idx.len() < 2 {
            continue;
        }
        let sd = sd_of(idx);
        let mut dists: Vec<f64> = Vec::new();
        for a in 0..idx.len() {
            for b in (a + 1)..idx.len() {
                dists.push(content_distance(&images[idx[a]], &images[idx[b]]));
            }
        }
        let d_med = median_of(&mut dists);
        eprintln!(
            "  组「{name}」（{} 张：{}）内容距离中位 {:.4}",
            idx.len(),
            idx.iter().map(|&i| images[i].stem.as_str()).collect::<Vec<_>>().join(", "),
            d_med
        );
        let ratio: Vec<String> = (0..6)
            .map(|k| {
                if overall_sd[k] > 1e-9 {
                    format!("{} {:.0}%", PARAM_NAMES[k], 100.0 * sd[k] / overall_sd[k])
                } else {
                    format!("{} -", PARAM_NAMES[k])
                }
            })
            .collect();
        eprintln!(
            "    组内 6 参数标准差 {}",
            (0..6)
                .map(|k| format!("{} {:.4}", PARAM_NAMES[k], sd[k]))
                .collect::<Vec<_>>()
                .join("  ")
        );
        eprintln!("    组内标准差 / 全体标准差：{}", ratio.join("  "));
    }
    eprintln!("  内容签名最近邻（不依赖事先知道哪几张同场景）：");
    for i in 0..n {
        let mut nn = (f64::MAX, 0usize);
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = content_distance(&images[i], &images[j]);
            if d < nn.0 {
                nn = (d, j);
            }
        }
        eprintln!(
            "    {:<10} → {:<10} 内容距离 {:.4}   需要的修正差（6 参数平均）{:.4}",
            images[i].stem,
            images[nn.1].stem,
            nn.0,
            y[i].dist(&y[nn.1])
        );
    }

    // ---------- 四、留一图验证 ----------
    // 所有臂共用同一个基座 T_i，差别只在给它加什么残差修正。
    // 三个表示并列，是为了先量清"哪种参数化才有 headroom"——若神谕都改善不了，
    // 那么"内容能不能预测"这个问题本身就是空的（不能拿一个无用的量去要求预测）。
    eprintln!("\n【四、留一图验证（唯一裁判）：所有臂共用同一个基座 T_i】");
    const ARM_NAMES: [&str; 10] = [
        "无修正",
        "平均(6参)",
        "内容(嵌套,6参)",
        "内容(偷看,6参)",
        "神谕(6参)",
        "平均(增益)",
        "内容(嵌套,增益)",
        "内容(偷看,增益)",
        "神谕(增益)",
        "神谕(完整曲线)",
    ];
    const A_NONE: usize = 0;
    const A_MEAN6: usize = 1;
    const A_NEST6: usize = 2;
    const A_CHERRY6: usize = 3;
    const A_ORACLE6: usize = 4;
    const A_MEANG: usize = 5;
    const A_NESTG: usize = 6;
    const A_CHERRYG: usize = 7;
    const A_ORACLEG: usize = 8;
    const A_ORACLEFULL: usize = 9;
    let n_arm = ARM_NAMES.len();
    let mut med: Vec<Vec<f64>> = vec![Vec::new(); n_arm];
    let mut avg: Vec<Vec<f64>> = vec![Vec::new(); n_arm];
    let mut errs: Vec<Vec<f64>> = vec![Vec::new(); n_arm];
    let mut picks: Vec<Vec<String>> = vec![Vec::new(); n_arm];
    // 每张图的左半（训练侧）：`(无修正, 6参神谕, 完整曲线神谕)` 的 ΔE00 中位。
    let mut left_side: Vec<(f64, f64, f64)> = Vec::with_capacity(n);

    for hole in 0..n {
        let cal: Vec<usize> = (0..n).filter(|j| *j != hole).collect();
        let base = &bases[hole];
        let base_out: Vec<[f32; 3]> = images[hole].right.iter().map(|r| base.apply(r.o)).collect();
        let want: Vec<[f32; 3]> = images[hole].right.iter().map(|r| r.t).collect();
        let eval = |f: &dyn Fn([f32; 3]) -> [f32; 3]| -> (f64, f64) {
            let mut v: Vec<f64> = base_out
                .iter()
                .zip(want.iter())
                .map(|(x, t)| deltae::delta_e_prophoto(f(*x), *t))
                .collect();
            let mean = mean_of(&v);
            (median_of(&mut v), mean)
        };

        // 标定集各图相对本基座的修正（用于平均臂与内容臂）
        let params_all: Vec<Vec<f64>> =
            (0..n).map(|j| corr[hole][j].params().to_vec()).collect();
        let gains_all: Vec<Vec<f64>> = (0..n).map(|j| gain[hole][j].to_vec()).collect();
        let oracle6 = corr[hole][hole];
        let oracleg = gain[hole][hole];

        let mut push = |slot: usize, r: (f64, f64), e: f64| {
            med[slot].push(r.0);
            avg[slot].push(r.1);
            errs[slot].push(e);
        };

        // (0) 无修正
        push(A_NONE, eval(&|x| x), f64::NAN);
        // (1)(4) 平均修正 / 神谕（6 参数）
        let mean6 = Corr::mean_of(&cal.iter().map(|&j| corr[hole][j]).collect::<Vec<_>>());
        let mean6c = mean6.curves();
        let oracle6c = oracle6.curves();
        push(A_MEAN6, eval(&|x| fit::forward_curves(&mean6c, x)), mean6.dist(&oracle6));
        push(A_ORACLE6, eval(&|x| fit::forward_curves(&oracle6c, x)), 0.0);
        // (2)(3) 内容预测（6 参数）
        for (slot, select_all) in [(A_NEST6, false), (A_CHERRY6, true)] {
            let select_set: Vec<usize> = if select_all { (0..n).collect() } else { cal.clone() };
            let (pred, pk) = predict_vec(&images, &params_all, &cal, &select_set, hole);
            let c = Corr::from_params([pred[0], pred[1], pred[2], pred[3], pred[4], pred[5]]);
            let curves = c.curves();
            push(slot, eval(&|x| fit::forward_curves(&curves, x)), c.dist(&oracle6));
            picks[slot].push(pk.iter().map(|&k| STAT_NAMES[k]).collect::<Vec<_>>().join(" / "));
        }
        // (5)(8) 平均 / 神谕（3 参数对角增益，施加在基座输出的线性值上）
        let mean_gain: Vec<f64> = (0..3)
            .map(|k| mean_of(&cal.iter().map(|&j| gain[hole][j][k]).collect::<Vec<_>>()))
            .collect();
        // 注意：`gain` 存的是 **log2** 增益（表格里显示的就是它），施加时必须先取 2 的幂。
        let apply_gain = |g: &[f64]| {
            let lin = [
                2f32.powf(g[0] as f32),
                2f32.powf(g[1] as f32),
                2f32.powf(g[2] as f32),
            ];
            move |x: [f32; 3]| [x[0] * lin[0], x[1] * lin[1], x[2] * lin[2]]
        };
        push(A_MEANG, eval(&apply_gain(&mean_gain)), gain_err(&mean_gain, &oracleg));
        push(A_ORACLEG, eval(&apply_gain(&oracleg)), 0.0);
        // (6)(7) 内容预测（3 参数增益）
        for (slot, select_all) in [(A_NESTG, false), (A_CHERRYG, true)] {
            let select_set: Vec<usize> = if select_all { (0..n).collect() } else { cal.clone() };
            let (pred, pk) = predict_vec(&images, &gains_all, &cal, &select_set, hole);
            push(slot, eval(&apply_gain(&pred)), gain_err(&pred, &oracleg));
            picks[slot].push(pk.iter().map(|&k| STAT_NAMES[k]).collect::<Vec<_>>().join(" / "));
        }
        // (9) 神谕（完整 129 点逐通道曲线）
        let fc = &full[hole];
        push(A_ORACLEFULL, eval(&|x| fit::forward_curves(fc, x)), f64::NAN);

        // 同一张图的**左半（训练）**上，神谕修正兑现得怎样？
        // 这是"需要的修正到底是不是一个全局每图对象"的直接检验：若左半大幅改善而右半不改善，
        // 那它就不是每图一个的函数，而是**逐区域**的——那样任何"每图一个修正"的方案都注定失败。
        {
            let stride = 4;
            let lt: Vec<[f32; 3]> =
                images[hole].left.iter().step_by(stride).map(|r| r.t).collect();
            let lo: Vec<[f32; 3]> = images[hole]
                .left
                .iter()
                .step_by(stride)
                .map(|r| base.apply(r.o))
                .collect();
            let med_of = |f: &dyn Fn([f32; 3]) -> [f32; 3]| -> f64 {
                let mut v: Vec<f64> = lo
                    .iter()
                    .zip(lt.iter())
                    .map(|(x, t)| deltae::delta_e_prophoto(f(*x), *t))
                    .collect();
                median_of(&mut v)
            };
            let o6 = oracle6.curves();
            left_side.push((
                med_of(&|x| x),
                med_of(&|x| fit::forward_curves(&o6, x)),
                med_of(&|x| fit::forward_curves(fc, x)),
            ));
        }

        eprintln!(
            "  留出 {:<10} 无修正 {:.3} | 6参: 平均 {:.3} 嵌套 {:.3} 偷看 {:.3} 神谕 {:.3} | \
             增益: 平均 {:.3} 嵌套 {:.3} 偷看 {:.3} 神谕 {:.3} | 完整曲线神谕 {:.3}",
            images[hole].stem,
            med[A_NONE][hole],
            med[A_MEAN6][hole],
            med[A_NEST6][hole],
            med[A_CHERRY6][hole],
            med[A_ORACLE6][hole],
            med[A_MEANG][hole],
            med[A_NESTG][hole],
            med[A_CHERRYG][hole],
            med[A_ORACLEG][hole],
            med[A_ORACLEFULL][hole],
        );
        if hole == 0 {
            eprintln!("      （各臂选中的统计量见汇总；这里只列第一张，避免刷屏）");
        }
    }

    eprintln!("\n=== 各臂汇总（{n} 张留出图，指标为 ΔE00 中位）===");
    eprintln!(
        "  {:<16} {:>10} {:>10} {:>12} {:>12}",
        "臂", "中位均值", "均值", "相对无修正", "修正参数偏差"
    );
    let baseline = mean_of(&med[A_NONE]);
    for (slot, name) in ARM_NAMES.iter().enumerate() {
        let m = mean_of(&med[slot]);
        let a = mean_of(&avg[slot]);
        let es: Vec<f64> = errs[slot].iter().copied().filter(|v| v.is_finite()).collect();
        eprintln!(
            "  {:<16} {:>10} {:>10} {:>12} {:>12}",
            name,
            num(m, 3),
            num(a, 3),
            format!("{:+.3}", baseline - m),
            if es.is_empty() { String::from("-") } else { num(mean_of(&es), 4) }
        );
    }
    eprintln!(
        "  「相对无修正」为正 = 比不加修正更好；headroom 就是各「神谕」臂的那一列。"
    );

    eprintln!("\n  内容臂各留出图选中的统计量（每张图 6 个参数各选一次，共 {} 次）：", 6 * n);
    for slot in [A_NEST6, A_CHERRY6, A_NESTG, A_CHERRYG] {
        let mut counts: Vec<(String, usize)> = Vec::new();
        for line in &picks[slot] {
            for name in line.split(" / ") {
                match counts.iter_mut().find(|(k, _)| k == name) {
                    Some((_, v)) => *v += 1,
                    None => counts.push((name.to_string(), 1)),
                }
            }
        }
        counts.sort_by_key(|a| std::cmp::Reverse(a.1));
        eprintln!(
            "    {:<16} {}",
            ARM_NAMES[slot],
            counts.iter().map(|(k, v)| format!("{k}×{v}")).collect::<Vec<_>>().join("  ")
        );
    }

    // ---------- 五、神谕修正的左半(训练) vs 右半(留出) ----------
    eprintln!("\n【五、神谕修正能不能在同一张图的另一半上兑现？（ΔE00 中位）】");
    eprintln!("  {:<11} {:>26} {:>26}", "图", "左半（训练侧）", "右半（留出侧）");
    for hole in 0..n {
        let (lb, l6, lf) = left_side[hole];
        eprintln!(
            "  {:<11} {:>26} {:>26}",
            images[hole].stem,
            format!("{lb:.3} → {l6:.3} → {lf:.3}"),
            format!(
                "{:.3} → {:.3} → {:.3}",
                med[A_NONE][hole], med[A_ORACLE6][hole], med[A_ORACLEFULL][hole]
            )
        );
    }
    let train6 = mean_of(&left_side.iter().map(|(b, o, _)| b - o).collect::<Vec<_>>());
    let trainfull = mean_of(&left_side.iter().map(|(b, _, f)| b - f).collect::<Vec<_>>());
    eprintln!(
        "  均值改善：左半（训练）6 参 {:+.3} / 完整曲线 {:+.3}；\
         右半（留出）6 参 {:+.3} / 完整曲线 {:+.3}",
        train6,
        trainfull,
        baseline - mean_of(&med[A_ORACLE6]),
        baseline - mean_of(&med[A_ORACLEFULL])
    );

    // ---------- 判定 ----------
    // 三步，顺序不能换：先问这条路有没有 headroom，再问内容预测吃到了多少，最后看相关性的门槛。
    let best_content_r = best_content.map(|b| b.2.abs()).unwrap_or(f64::NAN);
    let room6 = baseline - mean_of(&med[A_ORACLE6]);
    let roomg = baseline - mean_of(&med[A_ORACLEG]);
    let roomfull = baseline - mean_of(&med[A_ORACLEFULL]);
    let gain6 = mean_of(&med[A_MEAN6]) - mean_of(&med[A_NEST6]);
    let gaing = mean_of(&med[A_MEANG]) - mean_of(&med[A_NESTG]);
    let wins6 = med[A_NEST6].iter().zip(&med[A_MEAN6]).filter(|(a, b)| **a < **b).count();
    let winsg = med[A_NESTG].iter().zip(&med[A_MEANG]).filter(|(a, b)| **a < **b).count();

    eprintln!("\n========== 判定 ==========");
    eprintln!(
        "  ① headroom（神谕相对无修正，ΔE00 中位均值）：6 参数 {:+.3}，3 参数增益 {:+.3}，\
         完整逐通道曲线 {:+.3}",
        room6, roomg, roomfull
    );
    eprintln!(
        "  ② 内容预测相对「标定集平均」吃到多少：6 参数 {:+.3}（赢 {wins6}/{n} 张），\
         3 参数增益 {:+.3}（赢 {winsg}/{n} 张）；门槛 {:+.2}",
        gain6, gaing, MIN_GAIN
    );
    eprintln!(
        "  ③ 内容量最强相关 |r| = {:.3}（{}）；Bonferroni 校正门槛 ≈ {:.2}",
        best_content_r,
        best_content.map(|b| format!("{} × {}", b.0, b.1)).unwrap_or_default(),
        rc_bonf
    );

    let passes_gate = best_content_r.is_finite() && best_content_r >= rc_bonf;
    let captures = (roomg > MIN_GAIN && gaing > MIN_GAIN) || (room6 > MIN_GAIN && gain6 > MIN_GAIN);
    if passes_gate && captures {
        eprintln!(
            "  ⇒ **A 成立**：存在过校正门槛的强相关统计量，且在有 headroom 的表示上兑现。"
        );
    } else if roomg > MIN_GAIN || room6 > MIN_GAIN || roomfull > MIN_GAIN {
        eprintln!(
            "  ⇒ **A 不成立**：这条路上**确实有 headroom**（神谕 {room6:+.3} / {roomg:+.3} / {roomfull:+.3}，\
             即「每张图各补一个修正」是真能改善的），但**内容预测吃不掉它**——相对标定集平均只有 \
             {gain6:+.3}（6 参）/ {gaing:+.3}（增益）；最强相关 {best_content_r:.3} 也\
             {}。⇒ 需要的修正与这 {NSTAT} 个统计量无关，应转向 B：找一个还没被识别的逐图参数。",
            if passes_gate { "过了校正门槛但没兑现" } else { "没过 Bonferroni 校正门槛" }
        );
    } else {
        eprintln!(
            "  ⇒ **A 不成立，而且这条路的 headroom 接近零**：三种参数化的神谕修正都改善不了\
             （{room6:+.3} / {roomg:+.3} / {roomfull:+.3}）。在现行基座之上，「每张图再补一个修正」\
             不是一个有用的自由度——既然需要的修正几乎为零，问「它是否与内容相关」就没有意义。"
        );
    }
    eprintln!(
        "  边界：6 参数（编码域中点 + 两端的斜率）**概括不掉**需要的修正（神谕 {room6:+.3}，\
         比不加还差）；有 headroom 的是完整逐通道曲线（{roomfull:+.3}）与 3 参数增益（{roomg:+.3}）。\
         本诊断直接测了增益表示、以及完整曲线的低维投影；**没有**直接测「内容能不能预测整条曲线」\
         （n={n} 不足以验证一个 129 点对象）。"
    );
    eprintln!(
        "\n限制：n={n}；同内容组 {} 组（色温扫描 3 张 + 42 秒连拍 2 张）；\
         内容统计量只覆盖亮度/色度分布这类一阶量，「画面结构」只在内容签名里用到，未进入相关分析。",
        groups.len()
    );
}
