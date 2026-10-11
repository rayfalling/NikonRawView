//! 诊断 N：**跨图迁移惩罚（transfer penalty）的形状归因**。
//!
//! # 要回答的问题
//!
//! 端到端复测发现：每张图自己拟合「曲线 + 17³ LUT」时 ΔE00 中位 1.6~4.8，而把**一套**
//! 拟合好的变换套到别的图上时是 7.9~9.6——**差 2~5 倍**。也就是说我们那套变换本质上
//! 是「每张图各自的经验拟合」，不是一个可迁移的相机变换。
//!
//! 本诊断在**留一图**协议下把多出来的那部分误差拆开：
//!
//! ```text
//! (a) 自拟合：只用留出图自己的样本（左半）拟合曲线+LUT → 在它的右半上评估
//! (b) 跨图　：用**其余全部图**的样本（左半）池化拟合同一套 → 在留出图的右半上评估
//! penalty = (b) − (a)
//! ```
//!
//! 然后**逐级给跨图那套加自由度**（每一级都只用留出图的**左半**拟合，在**右半**评估），
//! 看每一级吃掉多少 penalty：
//!
//! | 级别 | 加什么 | 若吃掉大部分说明 |
//! |---|---|---|
//! | L1 | 3 参数对角增益 `w`（加在**矩阵之前**，即白平衡的位置） | 形状是逐通道增益型 |
//! | L2 | 1 条**共享**色调曲线（施加在跨图变换的输出上） | 形状是共享色调映射型 |
//! | L3 | 3 条**逐通道**色调曲线 | 需要逐通道曲线 |
//! | L4 | 完整 3×3（施加在输出上，线性混合） | 色彩转换环节仍有系统性差异 |
//! | L5 | 池化样本上重拟 **33³** LUT | 网格太稀，17³ 装不下 |
//!
//! 另有一列**累积链**（L1 → L1+L3 → L1+L3+L4），用来回答「几级合起来能吃掉多少」。
//!
//! # 这些级别是**神谕**，不是可达的改善
//!
//! L1~L4 都要用**留出图自己的参考值**来拟合——生产上拿到一张新图时并没有参考。
//! 它们衡量的是**penalty 的形状**（「如果模型多这一种自由度，误差会掉多少」），
//! **不是**「我们因此能改好多少」。任何把这里的数字当成可达增益的读法都是错的。
//!
//! # 覆盖度量（与表达力问题分开量）
//!
//! 对每张留出图另算一个量：它的样本落在 LUT 网格上时，有多少**质量**落在
//! **标定集一个样本都没有**的节点里（`unsupported_mass`）。17³ LUT 只存节点均值、
//! 空节点靠邻域填充，所以落在无支撑节点上的样本拿到的值本质上是**外推**。
//!
//! - `unsupported_mass` 高 → penalty 是**覆盖问题**（标定集没覆盖到留出图用到的颜色）
//! - `unsupported_mass` 低而 penalty 仍大 → 是**表达力问题**
//!
//! 两者的修法完全不同，所以必须分开量。**注意**：LUT 做三线性插值，节点无样本不等于
//! 该处一定差（邻居够近就行），所以 `unsupported_mass` 是覆盖问题的**上界**指标。
//!
//! # 纪律与限制
//!
//! - 五级**都吃不掉** penalty 本身就是最重要的结论：那说明缺的是一个还没被识别的环节。
//!   如实报，不硬找解释。
//! - 9 张图、两类场景；留一图仍然**共享场景类别**，「跨场景」这一维依然没被验证。
//! - 模型选择/协议都写在输出里；`(a)`/`(b)` 都用**空间留出**（左半拟合、右半评估），
//!   所以数字与 Lead 端到端那份（可能是全图拟合+全图评估）不会逐位相同。
//! - `NRV_N_STEMS=DSC_0001,DSC_0141` 可只跑指定图（迭代用）。

mod common;

use common::samples_dir;
use nikonrawview::camera;
use nikonrawview::color;
use nikonrawview::deltae;
use nikonrawview::fit::{self, Sample};
use nikonrawview::icc;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::mat3::{self, Mat3};
use nikonrawview::transform::BaseTransform;

// ---------------------------------------------------------------------------
// 口径与门槛（先定后看）
// ---------------------------------------------------------------------------

/// 默认 9 张图：覆盖两类场景（老图 + 同场景色温扫描）与白平衡范围两端。
const DEFAULT_STEMS: [&str; 9] = [
    "DSC_0001", "DSC_0141", "DSC_0314", "DSC_0507", "DSC_0540", "DSC_0563", "DSC_0567", "DSC_7737",
    "DSC_8562",
];

/// 采样步长。
const STEP: usize = 17;
/// 色调曲线的箱数（与生产一致）。
const CURVE_BINS: usize = 128;
/// 生产用的 LUT 边长。
const LUT_EDGE: usize = 17;
/// 高分辨率对照。
const LUT_EDGE_FINE: usize = 33;

/// 样本门槛：任一侧饱和/近黑都丢掉（与既有诊断一致）。
const CAM_LO: f32 = 1e-4;
const CAM_HI: f32 = 0.999;
const REF_LO: f32 = 1e-5;
const REF_HI: f32 = 0.999;
/// 半幅可用的最少样本数。
const MIN_SAMPLES: usize = 2000;

/// L1 的坐标下降：每轮三分搜索的迭代次数与轮数。
const W_TERNARY_ITERS: usize = 20;
const W_PASSES: usize = 3;
/// L1 拟合时对训练样本的抽稀步长（降开销）。
const W_FIT_STRIDE: usize = 8;
/// L1 每个分量的搜索区间。**放宽到 0.25~4.0**：这是神谕探针，区间卡住会把它的价值低估
/// （第一版取 0.4~2.0，实测有两个解顶在 2.0 上）。触边的解数会在输出里报出来。
const W_LO: f64 = 0.25;
const W_HI: f64 = 4.0;

/// 「某一级吃掉了大部分 penalty」的门槛（先定后看）。
const EATEN_STRONG: f64 = 0.6;
/// 算「吃掉比例」时的最小 penalty 分母；低于它的图只报绝对值，不报比例。
const MIN_PENALTY_FOR_FRACTION: f64 = 1.0;

// ---------------------------------------------------------------------------
// 基础工具
// ---------------------------------------------------------------------------

/// 一个样本：相机空间值、工作空间值（`M0·c`，即变换的输入）、参考值。
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

/// ΔE00 的中位数与 P95。
///
/// 闭包直接收**整个三元组**，由调用方显式挑字段（`r.0` 相机空间 / `r.1` 工作空间）——
/// 这样「拿工作空间值去当相机空间值用」这类错就写不出来（本诊断第一版正是这么错的：
/// L1 把 `M0·diag(w)` 施在了 `M0·c` 上，于是 L1 全线变差、w 也乱飞）。
fn de_stats<F: Fn(&Triple) -> [f32; 3]>(rows: &[Triple], f: F) -> (f64, f64) {
    let mut des: Vec<f64> = rows.iter().map(|r| deltae::delta_e_prophoto(f(r), r.2)).collect();
    if des.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    des.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p95 = des[((0.95 * des.len() as f64).ceil() as usize).saturating_sub(1).min(des.len() - 1)];
    (median_of(&mut des), p95)
}

/// 按生产公式把工作空间值变成 LUT 节点序号。
///
/// 与 `fit::lut_coord`（私有）同一套映射：**编码域**上等距取整。仪器自检里会用
/// `apply_lut` 在节点坐标处的取值验证这个映射（防止这里与实现漂移）。
fn lut_index(v: [f32; 3], edge: usize) -> usize {
    let q = |x: f32| ((fit::encode_srgb(x) * (edge - 1) as f32).round() as usize).min(edge - 1);
    (q(v[0]) * edge + q(v[1])) * edge + q(v[2])
}

/// 施加「矩阵之前的 3 通道增益」后再过变换：`T(M0·diag(w)·c)`。
fn apply_through_w(t: &BaseTransform, m0: Mat3, w: [f64; 3], c: [f32; 3]) -> [f32; 3] {
    let u = mat3::mul_vec(
        m0,
        [c[0] * w[0] as f32, c[1] * w[1] as f32, c[2] * w[2] as f32],
    );
    t.apply(u)
}

/// 一次坐标下降（目标 = 编码域平方误差）。
fn descend_w(t: &BaseTransform, m0: Mat3, sub: &[&Triple], w0: [f64; 3]) -> [f64; 3] {
    let loss = |w: [f64; 3]| -> f64 {
        let mut s = 0f64;
        let mut n = 0usize;
        for r in sub {
            let p = apply_through_w(t, m0, w, r.0);
            for (pk, tk) in p.iter().zip(r.2.iter()) {
                let e = (fit::encode_srgb(*pk) - fit::encode_srgb(*tk)) as f64;
                s += e * e;
                n += 1;
            }
        }
        s / n.max(1) as f64
    };
    let mut w = w0;
    for _ in 0..W_PASSES {
        for k in 0..3 {
            let (mut lo, mut hi) = (W_LO, W_HI);
            for _ in 0..W_TERNARY_ITERS {
                let m1 = lo + (hi - lo) / 3.0;
                let m2 = hi - (hi - lo) / 3.0;
                let mut a = w;
                a[k] = m1;
                let mut b = w;
                b[k] = m2;
                if loss(a) < loss(b) {
                    hi = m2;
                } else {
                    lo = m1;
                }
            }
            w[k] = (lo + hi) / 2.0;
        }
    }
    w
}

/// 在 `train` 上解矩阵之前的 3 通道增益 `w`，使 `T(M0·diag(w)·c)` 接近参考。
///
/// 经过曲线+LUT 之后目标对 `w` 是**非线性**的，没有闭式解。
///
/// **两次防护**，都是被实测打出来的：
/// 1. 目标函数取编码域（`encode_srgb`）的平方误差，不是线性域——线性域的残差被高光主导，
///    优化它会把高光拉对、中间调拉坏，而评审用的是 ΔE00；
/// 2. **多起点 + 按拟合半幅的 ΔE00 选解**。三分搜索假设沿每个坐标单峰，这里并不成立：
///    实测单起点会收敛到 `w_G≈2.8` 这类极端解（编码域损失更小，但 ΔE00 反而差 2.4 倍）。
///    用真正的评估指标在**拟合半幅**上挑，至少保证报出来的是「这个自由度最好能做到多少」。
fn fit_w_through(t: &BaseTransform, m0: Mat3, train: &[Triple]) -> ([f64; 3], f64) {
    let sub: Vec<&Triple> = train.iter().step_by(W_FIT_STRIDE).collect();
    if sub.is_empty() {
        return ([1.0, 1.0, 1.0], f64::NAN);
    }
    let starts = [
        [1.0f64, 1.0, 1.0],
        [0.8, 1.0, 1.25],
        [1.25, 1.0, 0.8],
        [1.0, 1.0, 1.0],
    ];
    let mut best: Option<([f64; 3], f64)> = None;
    for s in starts {
        let w = descend_w(t, m0, &sub, s);
        let de = de_stats(train, |r| apply_through_w(t, m0, w, r.0)).0;
        if best.map(|(_, b)| de < b).unwrap_or(true) {
            best = Some((w, de));
        }
    }
    best.unwrap_or(([1.0, 1.0, 1.0], f64::NAN))
}

/// 逐通道 1-D 曲线：把第 `k` 个通道单独喂给 `fit_curve`（它按通道分箱取中位数），
/// 复用生产同一套实现，避免自己重写一遍单调化。
fn fit_channel_curve(rows: &[Triple], k: usize) -> Vec<f32> {
    let s: Vec<Sample> = rows
        .iter()
        .map(|r| Sample { ours: [r.1[k]; 3], theirs: [r.2[k]; 3] })
        .collect();
    fit::fit_curve(&s, CURVE_BINS)
}

/// 施加逐通道曲线。
fn apply_channels(curves: &[Vec<f32>; 3], v: [f32; 3]) -> [f32; 3] {
    [
        fit::forward_curve(&curves[0], [v[0]; 3])[0],
        fit::forward_curve(&curves[1], [v[1]; 3])[1],
        fit::forward_curve(&curves[2], [v[2]; 3])[2],
    ]
}

// ---------------------------------------------------------------------------
// 逐图样本
// ---------------------------------------------------------------------------

struct ImageData {
    stem: String,
    /// 画面左半（用于拟合）与右半（用于评估）。
    left: Vec<Triple>,
    right: Vec<Triple>,
    secs: f64,
}

/// 解码一张图并收集（相机空间、工作空间、参考）三元组。
fn collect(stem: &str, m0: Mat3) -> Option<ImageData> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        common::skip("诊断 N", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }
    let t0 = std::time::Instant::now();
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
    Some(ImageData { stem: stem.to_string(), left, right, secs: t0.elapsed().as_secs_f64() })
}

/// 池化一组图的左半。
fn pool(images: &[ImageData], skip: usize) -> Vec<Triple> {
    let mut out = Vec::new();
    for (i, im) in images.iter().enumerate() {
        if i == skip {
            continue;
        }
        out.extend_from_slice(&im.left);
    }
    out
}

fn build_transform(samples: &[Sample], edge: usize, name: &str) -> BaseTransform {
    let curve = fit::fit_curve(samples, CURVE_BINS);
    let lut = fit::fit_lut(samples, &curve, edge);
    BaseTransform::from_fit(0x0000, name, Vec::new(), curve, edge, lut)
}

// ---------------------------------------------------------------------------
// 覆盖度量
// ---------------------------------------------------------------------------

/// 一张留出图相对标定集的网格覆盖情况。
struct Coverage {
    /// 标定集占用的节点数（无样本的节点不在内）。
    cal_nodes: usize,
    /// 标定集的节点占用率（占 edge³）。
    cal_fill: f64,
    /// 留出图用到的节点数。
    held_nodes: usize,
    /// 留出图**样本质量**落在「标定集无样本」节点里的比例（覆盖问题的上界指标）。
    unsupported_mass: f64,
    /// 只按节点数算的重叠率：`|H ∩ cal| / |H|`。
    node_overlap: f64,
}

fn coverage(cal: &[Sample], held: &[Triple], curve: &[f32], edge: usize) -> Coverage {
    let n = edge * edge * edge;
    let mut cal_cnt = vec![0u32; n];
    for s in cal {
        let c = fit::forward_curve(curve, s.ours);
        cal_cnt[lut_index(c, edge)] += 1;
    }
    let mut held_cnt = vec![0u32; n];
    for r in held {
        let c = fit::forward_curve(curve, r.1);
        held_cnt[lut_index(c, edge)] += 1;
    }
    let cal_nodes = cal_cnt.iter().filter(|v| **v > 0).count();
    let held_nodes = held_cnt.iter().filter(|v| **v > 0).count();
    let unsupported: u64 = (0..n)
        .filter(|i| held_cnt[*i] > 0 && cal_cnt[*i] == 0)
        .map(|i| held_cnt[i] as u64)
        .sum();
    let total: u64 = held_cnt.iter().map(|v| *v as u64).sum();
    let both = (0..n).filter(|i| held_cnt[*i] > 0 && cal_cnt[*i] > 0).count();
    Coverage {
        cal_nodes,
        cal_fill: cal_nodes as f64 / n as f64,
        held_nodes,
        unsupported_mass: if total == 0 { f64::NAN } else { unsupported as f64 / total as f64 },
        node_overlap: if held_nodes == 0 { f64::NAN } else { both as f64 / held_nodes as f64 },
    }
}

// ---------------------------------------------------------------------------
// 主测试
// ---------------------------------------------------------------------------

fn env_stems() -> Vec<String> {
    match std::env::var("NRV_N_STEMS") {
        Ok(v) if !v.trim().is_empty() => v
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect(),
        _ => DEFAULT_STEMS.iter().map(|s| s.to_string()).collect(),
    }
}

#[test]
fn transfer_penalty_shape() {
    eprintln!("=== 诊断 N：跨图迁移惩罚的形状归因 ===");
    let Some(cm) = camera::color_matrix("Nikon Z 8") else {
        eprintln!("跳过：机型表里没有 Z8 的色彩矩阵");
        return;
    };
    let m0 = cm.matrix;
    eprintln!(
        "协议：留一图；标定集 = 其余图的**左半**池化；每张留出图的每一级都在它的**左半**拟合、**右半**评估。"
    );
    eprintln!(
        "门槛（先定后看）：某一级吃掉 ≥ {:.0}% 的 penalty 才算「解释了这个形状」。",
        100.0 * EATEN_STRONG
    );

    let stems = env_stems();
    eprintln!("\n[0] 将要解码的图（{} 张）：{}", stems.len(), stems.join(" "));
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
    let decode_secs: f64 = images.iter().map(|im| im.secs).sum();
    eprintln!("  解码总耗时 {decode_secs:.0} s");

    // 每张留出图的结果
    struct Row {
        stem: String,
        self_de: f64,
        self_train_de: f64,
        transfer_de: f64,
        penalty: f64,
        /// 七个级别在**留出半幅**（评估）上的 ΔE00 中位，序同 `LEVEL_NAMES`。
        eval_de: [f64; 7],
        /// 同一批级别在**拟合半幅**上的 ΔE00 中位——与上面之差 = 这一级自身过拟合的程度。
        fit_de: [f64; 7],
        w: [f64; 3],
        cov: Coverage,
        transfer_rms: f64,
    }
    const LEVEL_NAMES: [&str; 7] = [
        "L1 增益（矩阵前 3 参数）",
        "L2 共享色调曲线",
        "L3 逐通道色调曲线",
        "L4 完整 3×3",
        "L5 LUT 17³→33³",
        "累积 L1+L3",
        "累积 L1+L3+L4",
    ];

    let mut rows: Vec<Row> = Vec::new();
    eprintln!("\n[1] 逐张留出图：拟合、评估");
    for held in 0..images.len() {
        let cal_rows = pool(&images, held);
        let cal: Vec<Sample> = to_samples(&cal_rows);
        let h = &images[held];
        eprint!(
            "  {}：自拟合…",
            h.stem
        );

        // (a) 自拟合：只用留出图左半；同时记下**训练误差**（在左半上评估），用来显出自拟合
        //     自己有多容易过拟合——留出半幅上的数字才是可比的。
        let self_t = build_transform(&to_samples(&h.left), LUT_EDGE, "self");
        let (self_de, _) = de_stats(&h.right, |r| self_t.apply(r.1));
        let (self_train_de, _) = de_stats(&h.left, |r| self_t.apply(r.1));

        // (b) 跨图：池化左半
        let tr = build_transform(&cal, LUT_EDGE, "transfer");
        let (transfer_de, _) = de_stats(&h.right, |r| tr.apply(r.1));
        let penalty = transfer_de - self_de;
        eprint!(" 跨图…");

        // 跨图变换在留出图上的相对线性残差（诊断用）
        let mut ss = 0f64;
        let mut tt = 0f64;
        for r in &h.right {
            let p = tr.apply(r.1);
            for (pk, tk) in p.iter().zip(r.2.iter()) {
                ss += ((*pk - *tk) as f64).powi(2);
                tt += (*tk as f64).powi(2);
            }
        }
        let transfer_rms = (ss / tt.max(1e-12)).sqrt();

        // L1：矩阵之前的 3 通道增益（用留出图左半解；**用相机空间字段 r.0**）
        let (w, _) = fit_w_through(&tr, m0, &h.left);
        let l1_c = |r: &Triple| apply_through_w(&tr, m0, w, r.0);
        let (l1, l1_fit) = (de_stats(&h.right, l1_c).0, de_stats(&h.left, l1_c).0);
        eprint!(" L1…");

        // L2：共享色调曲线（施加在跨图输出上）
        let post: Vec<Sample> = h
            .left
            .iter()
            .map(|r| Sample { ours: tr.apply(r.1), theirs: r.2 })
            .collect();
        let curve_shared = fit::fit_curve(&post, CURVE_BINS);
        let l2_c = |r: &Triple| fit::forward_curve(&curve_shared, tr.apply(r.1));
        let (l2, l2_fit) = (de_stats(&h.right, l2_c).0, de_stats(&h.left, l2_c).0);
        eprint!(" L2…");

        // L3：3 条逐通道曲线
        let shared_out: Vec<Triple> = h.left.iter().map(|r| (r.0, tr.apply(r.1), r.2)).collect();
        let ch_curves = [
            fit_channel_curve(&shared_out, 0),
            fit_channel_curve(&shared_out, 1),
            fit_channel_curve(&shared_out, 2),
        ];
        let l3_c = |r: &Triple| apply_channels(&ch_curves, tr.apply(r.1));
        let (l3, l3_fit) = (de_stats(&h.right, l3_c).0, de_stats(&h.left, l3_c).0);
        eprint!(" L3…");

        // L4：完整 3×3（施加在跨图输出上，线性最小二乘）
        let pairs: Vec<([f32; 3], [f32; 3])> =
            h.left.iter().map(|r| (tr.apply(r.1), r.2)).collect();
        let (l4, l4_fit) = match color::derive_matrix(&pairs) {
            Some(d) => {
                let c = |r: &Triple| mat3::mul_vec(d.matrix, tr.apply(r.1));
                (de_stats(&h.right, c).0, de_stats(&h.left, c).0)
            }
            None => (f64::NAN, f64::NAN),
        };
        eprint!(" L4…");

        // L5：池化样本上重拟 33³
        let tr_fine = build_transform(&cal, LUT_EDGE_FINE, "transfer-33");
        let l5_c = |r: &Triple| tr_fine.apply(r.1);
        let (l5, l5_fit) = (de_stats(&h.right, l5_c).0, de_stats(&h.left, l5_c).0);
        eprint!(" L5…");

        // 累积链：L1 → +逐通道曲线 → +3×3（每一级都在上一级的输出上重解）
        let stage1: Vec<Triple> =
            h.left.iter().map(|r| (r.0, apply_through_w(&tr, m0, w, r.0), r.2)).collect();
        let c_curves = [
            fit_channel_curve(&stage1, 0),
            fit_channel_curve(&stage1, 1),
            fit_channel_curve(&stage1, 2),
        ];
        let c13_c = |r: &Triple| apply_channels(&c_curves, apply_through_w(&tr, m0, w, r.0));
        let (c13, c13_fit) = (de_stats(&h.right, c13_c).0, de_stats(&h.left, c13_c).0);
        let stage2: Vec<([f32; 3], [f32; 3])> = stage1
            .iter()
            .map(|r| (apply_channels(&c_curves, r.1), r.2))
            .collect();
        let (c134, c134_fit) = match color::derive_matrix(&stage2) {
            Some(d) => {
                let c = |r: &Triple| {
                    mat3::mul_vec(
                        d.matrix,
                        apply_channels(&c_curves, apply_through_w(&tr, m0, w, r.0)),
                    )
                };
                (de_stats(&h.right, c).0, de_stats(&h.left, c).0)
            }
            None => (f64::NAN, f64::NAN),
        };
        eprintln!(" 完成");

        rows.push(Row {
            stem: h.stem.clone(),
            self_de,
            self_train_de,
            transfer_de,
            penalty,
            eval_de: [l1, l2, l3, l4, l5, c13, c134],
            fit_de: [l1_fit, l2_fit, l3_fit, l4_fit, l5_fit, c13_fit, c134_fit],
            w,
            cov: coverage(&cal, &h.left, &tr.curve, LUT_EDGE),
            transfer_rms,
        });
    }

    // ---- [2] 主表 ----
    eprintln!("\n[2] 逐图：(a) 自拟合 · (b) 跨图 · penalty · 各级后的 ΔE00 中位");
    eprintln!(
        "  {:<10} {:>9} {:>9} {:>8} {:>8} {:>9} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "留出图", "(a)训练", "(a)留出", "(b)跨图", "penalty", "+L1增益", "+L2共享", "+L3逐通道", "+L4 3×3", "+L5 33³", "L1+L3", "L1+L3+L4"
    );
    for r in &rows {
        eprintln!(
            "  {:<10} {:>9.3} {:>9.3} {:>8.3} {:>8.3} {:>9.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3}",
            r.stem,
            r.self_train_de,
            r.self_de,
            r.transfer_de,
            r.penalty,
            r.eval_de[0],
            r.eval_de[1],
            r.eval_de[2],
            r.eval_de[3],
            r.eval_de[4],
            r.eval_de[5],
            r.eval_de[6]
        );
    }
    let med = |v: Vec<f64>| median_of(&mut v.clone());
    let m_a = med(rows.iter().map(|r| r.self_de).collect());
    let m_at = med(rows.iter().map(|r| r.self_train_de).collect());
    let m_b = med(rows.iter().map(|r| r.transfer_de).collect());
    eprintln!(
        "  中位：自拟合（训练）{m_at:.3} · 自拟合（留出）{m_a:.3} · 跨图 {m_b:.3} · penalty {:.3}（跨图是自拟合留出的 {:.1}×）",
        m_b - m_a,
        m_b / m_a.max(1e-9)
    );
    eprintln!(
        "  注：自拟合的**训练**误差 {m_at:.3} 明显低于它在留出半幅上的 {m_a:.3}——我们这套「曲线+17³ LUT」"
    );
    eprintln!("  自己就会过拟合（4913 个节点摊在约 8 万个样本上，每节点十来个样本）。所以与 Lead 那份");
    eprintln!("  「同图拟合+同图评估」的 (a) 1.6~4.8 相比，本表的 (a) 更保守；penalty 的**比值**才是同口径的。");

    // ---- [3] 各级吃掉多少 penalty ----
    eprintln!("\n[3] 每一级吃掉的 penalty（正 = 改善；比例 = 该级吃掉的 / penalty）");
    let levels: [(&str, usize); 7] = [
        ("L1 增益（矩阵前 3 参数）", 0),
        ("L2 共享色调曲线", 1),
        ("L3 逐通道色调曲线", 2),
        ("L4 完整 3×3", 3),
        ("L5 LUT 17³→33³", 4),
        ("累积 L1+L3", 5),
        ("累积 L1+L3+L4", 6),
    ];
    eprintln!(
        "  {:<24} {:>12} {:>14} {:>12} {:>10} {:>10}",
        "级别", "吃掉(中位)", "吃掉比例(中位)", "剩余(均值)", "变差图数", "到基线图数"
    );
    let mut summary: Vec<(String, f64, f64, f64, usize)> = Vec::new();
    for (name, idx) in levels {
        let get = |r: &Row| r.eval_de[idx];
        let eaten: Vec<f64> = rows.iter().map(|r| r.transfer_de - get(r)).collect();
        // 比例只在 penalty 足够大的图上算：有 2/9 张的 penalty 近似为 0（甚至为负），
        // 拿它当分母会把「比例」变成噪声放大（−309%、−522% 这种数字没有意义）。
        let usable: Vec<f64> = rows
            .iter()
            .filter(|r| r.penalty >= MIN_PENALTY_FOR_FRACTION)
            .map(|r| (r.transfer_de - get(r)) / r.penalty)
            .collect();
        let neg = eaten.iter().filter(|v| **v < 0.0).count();
        let reaches = rows.iter().filter(|r| get(r) <= r.self_de).count();
        let (me, mf) = (median_of(&mut eaten.clone()), median_of(&mut usable.clone()));
        eprintln!(
            "  {:<24} {:>12.3} {:>13.1}% {:>12.3} {:>10} {:>10}",
            name,
            me,
            100.0 * mf,
            mean_of(&rows.iter().map(get).collect::<Vec<_>>()),
            neg,
            reaches
        );
        summary.push((name.to_string(), me, mf, 0.0, neg));
    }
    let degenerate = rows.iter().filter(|r| r.penalty < MIN_PENALTY_FOR_FRACTION).count();
    eprintln!(
        "  注：「吃掉比例」只在 penalty ≥ {MIN_PENALTY_FOR_FRACTION:.1} 的 {} 张图上算（排除 {degenerate} 张 penalty≈0 的图）；",
        rows.len() - degenerate
    );
    eprintln!("      「到基线图数」= 该级之后 ΔE00 已 ≤ 自拟合留出基线的图数。");

    // ---- [3b] 每一级自己的过拟合：拟合半幅 vs 留出半幅 ----
    eprintln!("\n[3b] 每一级自身的过拟合：**拟合半幅**上的 ΔE00 中位 vs **留出半幅**上的（差越大越可疑）");
    eprintln!(
        "  {:<24} {:>14} {:>14} {:>10}",
        "级别", "拟合半幅(中位)", "留出半幅(中位)", "差"
    );
    for (idx, name) in LEVEL_NAMES.iter().enumerate() {
        let fit: Vec<f64> = rows.iter().map(|r| r.fit_de[idx]).collect();
        let ev: Vec<f64> = rows.iter().map(|r| r.eval_de[idx]).collect();
        let (mf, me) = (median_of(&mut fit.clone()), median_of(&mut ev.clone()));
        eprintln!("  {name:<24} {mf:>14.3} {me:>14.3} {:>10.3}", me - mf);
    }
    {
        let fit: Vec<f64> = rows.iter().map(|r| r.self_train_de).collect();
        let ev: Vec<f64> = rows.iter().map(|r| r.self_de).collect();
        let (mf, me) = (median_of(&mut fit.clone()), median_of(&mut ev.clone()));
        eprintln!("  {:<24} {:>14.3} {:>14.3} {:>10.3}", "（对照）(a) 自拟合", mf, me, me - mf);
    }

    // ---- [4] L1 解出的 w（可与诊断 M 的模型对照） ----
    eprintln!("\n[4] L1 解出的 w（矩阵之前的三通道增益，=1 表示不需要修正）");
    eprintln!(
        "  {:<10} {:>10} {:>10} {:>10}   {:>10} {:>10}",
        "留出图", "w_R", "w_G", "w_B", "w_R/w_G", "w_B/w_G"
    );
    for r in &rows {
        eprintln!(
            "  {:<10} {:>10.5} {:>10.5} {:>10.5}   {:>10.5} {:>10.5}",
            r.stem,
            r.w[0],
            r.w[1],
            r.w[2],
            r.w[0] / r.w[1],
            r.w[2] / r.w[1]
        );
    }
    let w_b: Vec<f64> = rows.iter().map(|r| r.w[2]).collect();
    let w_r: Vec<f64> = rows.iter().map(|r| r.w[0]).collect();
    eprintln!(
        "  范围：w_R {:.4}~{:.4}（{:.2}×）· w_B {:.4}~{:.4}（{:.2}×）",
        w_r.iter().cloned().fold(f64::INFINITY, f64::min),
        w_r.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        w_r.iter().cloned().fold(f64::NEG_INFINITY, f64::max) / w_r.iter().cloned().fold(f64::INFINITY, f64::min),
        w_b.iter().cloned().fold(f64::INFINITY, f64::min),
        w_b.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        w_b.iter().cloned().fold(f64::NEG_INFINITY, f64::max) / w_b.iter().cloned().fold(f64::INFINITY, f64::min)
    );
    let bound = 1e-3;
    let at_bound = rows
        .iter()
        .filter(|r| r.w.iter().any(|v| (*v - W_LO).abs() < bound || (*v - W_HI).abs() < bound))
        .count();
    eprintln!(
        "  触到搜索边界（{W_LO}~{W_HI}）的解：{at_bound}/{} 张{}",
        rows.len(),
        if at_bound > 0 { " ⇒ 这些图的 gain 需求超出区间，说明该级被区间限制住了，实际价值可能更高" } else { "" }
    );

    // ---- [5] 覆盖度量 ----
    eprintln!("\n[5] LUT 网格覆盖：留出图的样本有多少质量落在**标定集无样本**的节点里");
    eprintln!(
        "  {:<10} {:>12} {:>10} {:>14} {:>12} {:>10}",
        "留出图", "标定占用率", "留出节点数", "无支撑质量占比", "节点重叠率", "penalty"
    );
    for r in &rows {
        eprintln!(
            "  {:<10} {:>11.2}% {:>10} {:>13.2}% {:>11.2}% {:>10.3}",
            r.stem,
            100.0 * r.cov.cal_fill,
            r.cov.held_nodes,
            100.0 * r.cov.unsupported_mass,
            100.0 * r.cov.node_overlap,
            r.penalty
        );
    }
    let uns: Vec<f64> = rows.iter().map(|r| r.cov.unsupported_mass).collect();
    let pens: Vec<f64> = rows.iter().map(|r| r.penalty).collect();
    let pen_after_l1: Vec<f64> = rows.iter().map(|r| r.transfer_de - r.eval_de[0]).collect();
    let pen_after_all: Vec<f64> = rows.iter().map(|r| r.transfer_de - r.eval_de[6]).collect();
    eprintln!(
        "  中位：标定占用率 {:.2}%（{} 格中 {} 格有样本）· 无支撑质量 {:.2}%",
        100.0 * median_of(&mut rows.iter().map(|r| r.cov.cal_fill).collect::<Vec<_>>()),
        LUT_EDGE * LUT_EDGE * LUT_EDGE,
        rows[0].cov.cal_nodes,
        100.0 * median_of(&mut uns.clone())
    );
    eprintln!(
        "  相关系数 r(无支撑质量, penalty) = {:+.3} · r(无支撑质量, L1 之后的剩余) = {:+.3} · r(无支撑质量, 全部之后剩余) = {:+.3}",
        pearson(&uns, &pens),
        pearson(&uns, &pen_after_l1),
        pearson(&uns, &pen_after_all)
    );

    // ---- [6] 判定 ----
    let best = summary
        .iter()
        .filter(|(n, _, _, _, _)| !n.starts_with("累积"))
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    let best_cum = summary.iter().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!("\n[6] 判定");
    eprintln!(
        "  1) penalty 量级：中位 {:.3}（跨图 {:.3} vs 自拟合 {:.3}，{:.1}×）",
        m_b - m_a,
        m_b,
        m_a,
        m_b / m_a.max(1e-9)
    );
    if let Some((name, me, mf, _, neg)) = best {
        eprintln!(
            "  2) 单级最好者是「{name}」：吃掉 {me:.3}（{:.0}%），但有 {neg}/{} 张图反而变差",
            100.0 * mf,
            rows.len()
        );
        if *mf >= 1.0 {
            eprintln!(
                "     ⇒ 它不只是「吃掉 penalty」，而是**比自拟合基线还好**（>100%）。原因是自拟合基线本身在过拟合："
            );
            eprintln!(
                "     它把 4913 个节点拟合在约 8 万个样本上，而这一级只有几百个参数、又借了跨图变换的正则化。"
            );
            eprintln!("     所以 penalty 里有一部分是「自拟合的过拟合」，不全是「跨图失败」。");
        } else if *mf >= EATEN_STRONG {
            eprintln!("     ⇒ 该级**解释了主要形状**（≥{:.0}%）。", 100.0 * EATEN_STRONG);
        } else {
            eprintln!(
                "     ⇒ 没有哪一级吃掉 ≥{:.0}%：**单一自由度解释不了 penalty**。",
                100.0 * EATEN_STRONG
            );
        }
    }
    if let Some((name, me, mf, _, _)) = best_cum {
        eprintln!("  3) 累积最好者是「{name}」：吃掉 {me:.3}（{:.0}%）", 100.0 * mf);
    }
    let med_uns = median_of(&mut uns.clone());
    let r_uns = pearson(&uns, &pens);
    if med_uns > 0.5 && r_uns > 0.5 {
        eprintln!(
            "  4) 覆盖：无支撑质量中位 {:.0}% 且与 penalty 正相关（r {r_uns:+.2}）⇒ penalty 里**有覆盖问题的成分**。",
            100.0 * med_uns
        );
    } else if med_uns < 0.2 {
        eprintln!(
            "  4) 覆盖：无支撑质量中位只有 {:.0}% ⇒ **不是覆盖问题**（标定集基本覆盖了留出图用到的颜色），",
            100.0 * med_uns
        );
        eprintln!("     那么 penalty 属于**表达力/未知环节**。");
    } else {
        eprintln!(
            "  4) 覆盖：无支撑质量中位 {:.0}%、与 penalty 相关 {r_uns:+.2} ⇒ 覆盖与表达力**各占一部分**，不能只归给一边。",
            100.0 * med_uns
        );
    }
    eprintln!();
    eprintln!("  限制（必须随结论一起转述）：");
    eprintln!("  · L1~L4 都是**神谕**：用留出图自己的参考值拟合，生产上拿不到。它们量的是 penalty 的**形状**，不是可达增益。");
    eprintln!("  · {} 张图、两类场景（老图 + 一组同场景色温扫描）；留一图仍然**共享场景类别**，跨场景这一维未被验证。", rows.len());
    eprintln!("  · (a)/(b) 都用空间留出（左半拟合、右半评估），与 Lead 端到端那份（可能全图拟合+全图评估）不会逐位相同。");
    eprintln!("  · LUT 三线性插值：节点无样本不等于该处一定差，`无支撑质量`是覆盖问题的**上界**指标。");
    eprintln!("  · L5 不是嵌套级（它重拟了一整棵树），只作对照。");
    eprintln!();
    eprintln!("  附：跨图变换在留出图上的相对线性残差中位 {:.4}", median_of(&mut rows.iter().map(|r| r.transfer_rms).collect::<Vec<_>>()));
}

/// 仪器自检（不依赖样本）：
///
/// 1. `lut_index` 与 `apply_lut` 的节点映射必须一致（防止我这里与实现漂移）；
/// 2. `fit_w_through` 在「变换是恒等、目标是 w⊙c」时应解回 w；
/// 3. 逐通道曲线在恒等数据上应接近恒等。
#[test]
fn instrument_checks() {
    let m0 = camera::color_matrix("Nikon Z 8").expect("Z8 应有固定矩阵").matrix;

    // 1) 节点映射：在节点坐标处 apply_lut 必须精确返回该节点的值
    let edge = 5usize;
    let n = edge * edge * edge;
    let mut lut = Vec::with_capacity(n * 3);
    for i in 0..n {
        lut.push(i as f32 * 0.001);
        lut.push(i as f32 * 0.002);
        lut.push(i as f32 * 0.003);
    }
    for i in [0usize, 1, 7, 60, n - 1] {
        let coord = |axis: usize| -> f32 {
            let idx = match axis {
                0 => i / (edge * edge),
                1 => (i / edge) % edge,
                _ => i % edge,
            };
            fit::decode_srgb(idx as f32 / (edge - 1) as f32)
        };
        let v = [coord(0), coord(1), coord(2)];
        assert_eq!(lut_index(v, edge), i, "节点映射应与 apply_lut 一致（{v:?}）");
        let got = fit::apply_lut(&lut, edge, v);
        for k in 0..3 {
            let want = lut[i * 3 + k];
            assert!((got[k] - want).abs() < 1e-4, "节点 {i} 通道 {k}：{got:?} vs {want}");
        }
    }

    // 2) w 求解器：变换取恒等、目标取 w⊙c 时应收敛到 w
    let ident = BaseTransform::new(0, "identity");
    let w_true = [1.15f64, 1.0, 0.85];
    let mut train: Vec<Triple> = Vec::new();
    for i in 0..40 {
        let c = [0.05 + i as f32 * 0.02, 0.10 + i as f32 * 0.017, 0.07 + i as f32 * 0.013];
        // 恒等变换下，工作空间值 = M0·c；目标 = M0·(w⊙c)
        let ours = mat3::mul_vec(m0, c);
        let t = mat3::mul_vec(m0, [c[0] * w_true[0] as f32, c[1] * w_true[1] as f32, c[2] * w_true[2] as f32]);
        train.push((c, ours, t));
    }
    let (w, _) = fit_w_through(&ident, m0, &train);
    for k in 0..3 {
        assert!(
            (w[k] - w_true[k]).abs() < 0.02,
            "w 求解器应收敛到 {w_true:?}，得到 {w:?}"
        );
    }

    // 3) 逐通道曲线在恒等数据上应接近恒等
    let rows: Vec<Triple> = (0..40)
        .map(|i| {
            let v = [0.05 + i as f32 * 0.02; 3];
            (v, v, v)
        })
        .collect();
    let curve = fit_channel_curve(&rows, 0);
    // 只查样本覆盖到的范围（`rows` 的值域是 0.05~0.83）；两端之外曲线按端点延拓，不参与这条自检。
    for x in [0.1f32, 0.35, 0.6, 0.8] {
        let y = fit::forward_curve(&curve, [x; 3])[0];
        assert!((y - x).abs() < 0.02, "恒等数据上曲线应近似恒等：{x} → {y}");
    }

    // 4) 覆盖度量：拿自己的样本当标定集时，无支撑质量必须为 0
    let held: Vec<Triple> = rows.clone();
    let cal: Vec<Sample> = to_samples(&rows);
    let cov = coverage(&cal, &held, &curve, 4);
    assert!(
        cov.unsupported_mass < 1e-12,
        "自覆盖时无支撑质量应为 0，得到 {}",
        cov.unsupported_mass
    );

    // 5) 网格分辨率本身不是实现问题：在**解析已知且足够密**的样本上，33³ 必须优于 17³。
    //    这一步把「L5 变差」是数据的锅还是实现的锅分开——真实数据上节点稀疏、节点均值有噪声，
    //    这里没有噪声，只剩三线性插值误差。
    let synth: Vec<Sample> = {
        let mut v = Vec::with_capacity(64 * 64 * 64);
        for i in 0..64 {
            for j in 0..64 {
                for k in 0..64 {
                    let ours = [
                        i as f32 / 63.0,
                        j as f32 / 63.0,
                        k as f32 / 63.0,
                    ];
                    let theirs = [
                        ours[0].powf(1.6),
                        ours[1].powf(1.6),
                        ours[2].powf(1.6),
                    ];
                    v.push(Sample { ours, theirs });
                }
            }
        }
        v
    };
    let synth_curve = fit::fit_curve(&synth, CURVE_BINS);
    let test: Vec<Sample> = {
        let mut v = Vec::with_capacity(16 * 16 * 16);
        for i in 0..16 {
            for j in 0..16 {
                for k in 0..16 {
                    let ours = [
                        (i as f32 + 0.5) / 16.0,
                        (j as f32 + 0.5) / 16.0,
                        (k as f32 + 0.5) / 16.0,
                    ];
                    let theirs = [ours[0].powf(1.6), ours[1].powf(1.6), ours[2].powf(1.6)];
                    v.push(Sample { ours, theirs });
                }
            }
        }
        v
    };
    let err = |edge: usize| -> f64 {
        let lut = fit::fit_lut(&synth, &synth_curve, edge);
        let mut s = 0f64;
        let mut n = 0usize;
        for t in &test {
            let p = fit::apply_lut(&lut, edge, fit::forward_curve(&synth_curve, t.ours));
            for (pk, tk) in p.iter().zip(t.theirs.iter()) {
                s += ((*pk - *tk) as f64).abs();
                n += 1;
            }
        }
        s / n as f64
    };
    let (e17, e33) = (err(17), err(33));
    assert!(
        e33 < e17,
        "密集无噪声数据上 33³ 应优于 17³：{e33} vs {e17}（若不然，问题在实现而不在数据）"
    );
}
