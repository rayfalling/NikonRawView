//! 诊断 E：**色度差的分阶段归因**——它在管线的哪一级产生。
//!
//! # 要回答的问题
//!
//! 三张参考图的 ΔE00 中位数是 `DSC_0001` 1.698（竖构图近中性）、`DSC_0141` 3.543、
//! `DSC_8562` 4.767（后两张横构图彩色），都没到「中位数 ≤ 1.0」。后两张上「我们原始
//! 解码的 C* 中位」只有参考的 0.53~0.63 倍。本诊断把这个色度差**定位到管线的某一级**：
//!
//! ```text
//! ① 原始解码（相机空间）            libraw::decode_camera_linear
//! ② 施加相机矩阵后（ProPhoto 线性）  libraw::decode_working_space
//! ③ 仅曲线之后                      fit::forward_curve(transform.curve, v)
//! ④ 曲线 + LUT 之后                 transform.apply(v)
//! ```
//!
//! ②③④ 由脚手架 `fitted_fixture_for` 的 `DiagSample` 与 `BaseTransform` 直接给出；
//! ① 需要自己再解一次 NEF（`decode_camera_linear`，每张约 20 秒、内部实际解码两遍）。
//!
//! # 判定用哪把尺子
//!
//! 三把尺子按可信度递增，判定只用最后一把：
//!
//! 1. **C\* 比**（任务书要求报的那个）：受亮度污染——② 还没施加色调曲线，画面整体偏暗，
//!    C* 自然就低。实测只加曲线（③）就能把 `DSC_0141` 的全样本 C* 比从 0.535 抬到 0.978，
//!    可见 ② 的低 C* 大半是亮度造成的。**它不能单独用来判色度。**
//! 2. `Δlog2(R/G)`、`Δlog2(B/G)`：尺度不变，任何逐像素曝光缩放都不改变它。
//! 3. 同样两个量，但只取**明度已配对**（`|ΔL*| ≤ 2`）的样本：把"我们的亮度本来就和参考
//!    不一样"这一条也排除掉，剩下的差只能是色度本身。**判定用它。**
//!
//! # ① 行怎么读
//!
//! ① 是相机空间，「把相机 RGB 直接按 ProPhoto 基解读」得到的 Lab 没有跨基色的色度含义，
//! 那一行的绝对值只能当**基线**用（看矩阵这一步把什么改掉了）。真正有含义的「①表现如何」
//! 是两行【上界】：把①用**理想的线性色彩转换**映射到参考上，能做到多好——
//!
//! - `① + 最优 3×3`：最小二乘意义下最好的相机矩阵。它若能把色度差消掉，说明差在
//!   **矩阵这一步**（即 ①②之间，色彩转换环节选错了变换）；它若消不掉，说明差在
//!   **① 本身或其上游**（解码、白平衡、黑电平、以及参考端与解码端之间的非线性差异）。
//!   同集求解的那一行是乐观上界；另有一行**左半求解→右半评估**的空间留出，判定以它为准。
//! - `① + 最优对角增益`：只允许逐通道增益（白平衡型修正）。它若够用，说明差的形态是
//!   通道缩放而不是矩阵形状。
//!
//! 用法：默认跑 `DSC_0001` / `DSC_0141` / `DSC_8562` 三张；`NRV_DIAG_STEM=...` 只跑一张。

mod common;

use common::{fitted_fixture_for, quantiles, samples_dir};
use nikonrawview::{camera, color, deltae, fit, libraw, mat3};
use std::time::Instant;

// ---------------------------------------------------------------------------
// 口径与常量
// ---------------------------------------------------------------------------

/// ΔE00 的达标线（与 `deltae::Report` 一致）。
const TARGET_MEDIAN: f64 = 1.0;

/// C* 比的门槛：参考 C* 低于此值的样本接近中性，比值是小分母上的噪声。
const CSTAR_GATE: f64 = 5.0;

/// `R/G`、`B/G` 差的门槛：三通道都须为正且不太小，否则比值无意义。
const RATIO_GATE: f32 = 1e-4;

/// ①→② 线性拟合只采信**未饱和也未压黑**的样本：两端的裁切会把裁切行为误当成矩阵。
const LINEAR_GATE: f32 = 0.002;
const LINEAR_GATE_HI: f32 = 0.95;

/// 「拟合前就已存在明显色度差」的判定门槛，单位 **log2**：明度配对的有色样本上，
/// `|Δlog2(R/G)|` 或 `|Δlog2(B/G)|` 超过它才算明显（0.10 ≈ 比值差 7%）。
///
/// 刻意不用 C* 比做这个判定：C* 受亮度影响，同一个色度在更暗的样本上 C* 更小，
/// 于是"色调曲线还没施加"这一事实本身就会把 C* 比拉偏——那是明度差，不是色度差。
const PREEXISTING_GATE: f64 = 0.10;

/// 「明度已配对」的门槛：`|ΔL*|` 不超过它才认为亮度对上了（L* 差 2 上下即肉眼可辨的明度差）。
const LUM_MATCH: f32 = 2.0;

/// 八段色相的中心（度）——与 `render.rs` 的 `BAND_CENTERS` 同一套分段，
/// 这样分段表与 4.2 的色相混合器说的是同一件事。
const BAND_CENTERS: [f32; 8] = [0.0, 30.0, 60.0, 120.0, 180.0, 240.0, 270.0, 300.0];
const BAND_NAMES: [&str; 8] = ["红", "橙", "黄", "绿", "青", "蓝", "紫", "品红"];

/// 任务 3.3 记录的实测相机矩阵（相机 → ProPhoto），只用于对照打印，不参与任何判定。
const DOCUMENTED_MATRIX: [[f32; 3]; 3] = [
    [0.68890, 0.32861, -0.01758],
    [0.01733, 1.26110, -0.27849],
    [-0.01457, -0.10767, 1.12218],
];

/// 分级表里各行的下标。
const L_CAM: usize = 0;
const L_WORK: usize = 1;
const L_CURVE: usize = 2;
const L_FULL: usize = 3;
const L_ORACLE_MATRIX: usize = 4;
const L_ORACLE_GAIN: usize = 5;

// ---------------------------------------------------------------------------
// 基础换算
// ---------------------------------------------------------------------------

fn lab(v: [f32; 3]) -> [f32; 3] {
    deltae::lab_from_prophoto_linear(v)
}

/// Lab 彩度 C*。
fn chroma(l: [f32; 3]) -> f64 {
    ((l[1] as f64).powi(2) + (l[2] as f64).powi(2)).sqrt()
}

/// Lab 色相角（度，0..360）。
fn hue_deg(l: [f32; 3]) -> f32 {
    let h = (l[2] as f64).atan2(l[1] as f64).to_degrees();
    (if h < 0.0 { h + 360.0 } else { h }) as f32
}

/// 解码层的 16 位输出 → 0..1。与脚手架 `fitted_fixture_for` 的换算保持一致。
fn norm96(p: [u16; 3]) -> [f32; 3] {
    [p[0] as f32 / 65535.0, p[1] as f32 / 65535.0, p[2] as f32 / 65535.0]
}

/// 逐通道对数比 `[log2(R/G), log2(B/G)]`；任一通道不满足门槛时返回 `None`。
///
/// 取对数是因为比值的偏离是**乘性**的：`log2` 差 +1 就是「我们的 R/G 是参考的两倍」，
/// 在不同亮度上可直接比较。
fn log_ratios(v: [f32; 3]) -> Option<[f64; 2]> {
    if v.iter().any(|x| !x.is_finite() || *x <= RATIO_GATE) {
        return None;
    }
    Some([
        (v[0] as f64 / v[1] as f64).log2(),
        (v[2] as f64 / v[1] as f64).log2(),
    ])
}

fn median_of(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[v.len() / 2]
}

/// 逐通道增益的最小二乘最优解 `g_c = Σ(o·t) / Σ(o²)`（每个通道各自缩放）。
fn best_gain(pairs: &[([f32; 3], [f32; 3])]) -> [f64; 3] {
    let mut g = [1f64; 3];
    for (c, gc) in g.iter_mut().enumerate() {
        let (mut num, mut den) = (0f64, 0f64);
        for (a, b) in pairs {
            num += a[c] as f64 * b[c] as f64;
            den += (a[c] as f64) * (a[c] as f64);
        }
        if den > 0.0 {
            *gc = num / den;
        }
    }
    g
}

fn max_abs(a: [f32; 3], b: [f32; 3]) -> f32 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).fold(0f32, f32::max)
}

fn mat_max_diff(a: &[[f32; 3]; 3], b: &[[f32; 3]; 3]) -> f32 {
    let mut worst = 0f32;
    for (ra, rb) in a.iter().zip(b.iter()) {
        for (x, y) in ra.iter().zip(rb.iter()) {
            worst = worst.max((x - y).abs());
        }
    }
    worst
}

fn fmt_mat(m: &[[f32; 3]; 3]) -> String {
    m.iter()
        .map(|r| format!("[{:>9.5} {:>9.5} {:>9.5}]", r[0], r[1], r[2]))
        .collect::<Vec<_>>()
        .join(" ")
}

/// 某个色相落在哪一段（取最近的段中心，与 `render.rs` 的权重函数同一规则）。
fn band_index(hue: f32) -> usize {
    let h = hue.rem_euclid(360.0);
    let mut best = (f32::MAX, 0usize);
    for (i, c) in BAND_CENTERS.iter().enumerate() {
        let mut d = (h - c).abs();
        if d > 180.0 {
            d = 360.0 - d;
        }
        if d < best.0 {
            best = (d, i);
        }
    }
    best.1
}

// ---------------------------------------------------------------------------
// 表格的排版：中日韩字符在等宽终端里占两列，只按字符数补空格会错位
// ---------------------------------------------------------------------------

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

fn num(v: f64, prec: usize) -> String {
    if v.is_finite() {
        format!("{v:.prec$}")
    } else {
        String::from("-")
    }
}

fn signed(v: f64, prec: usize) -> String {
    if v.is_finite() {
        format!("{v:+.prec$}")
    } else {
        String::from("-")
    }
}

// ---------------------------------------------------------------------------
// 一行统计
// ---------------------------------------------------------------------------

/// 某一级与参考的逐样本比较结果。
struct Row {
    label: &'static str,
    de_median: f64,
    de_p95: f64,
    /// `中位 C*(我们) / 中位 C*(参考)`，**全部样本**——与「原始解码的 C* 中位只有参考的
    /// 0.53~0.63 倍」是同一个口径。
    cstar_ratio_all: f64,
    /// 逐样本 `C*_ours / C*_ref` 的中位数（只统计参考 C* ≥ 门槛的样本，避开小分母）。
    cstar_ratio: f64,
    cstar_n: usize,
    /// `log2((R/G)_ours / (R/G)_ref)` 的中位数（只统计三通道皆正的样本）。
    rg_dev: f64,
    bg_dev: f64,
    ratio_n: usize,
    /// 同样的两个量，但只统计**参考 C* ≥ 门槛**的有色样本。
    ///
    /// 全样本的中位数会被大量近中性像素稀释：那些像素的 R/G、B/G 本来就近 1，偏差再大也
    /// 不说明色度问题。判断色度必须看有色样本。
    rg_dev_c: f64,
    bg_dev_c: f64,
    ratio_c_n: usize,
    /// 全部误差量里明度项的**合计占比**（`deltae::summarize_prophoto` 的可加口径）。
    lightness_share: f64,
}

impl Row {
    fn chroma_share(&self) -> f64 {
        1.0 - self.lightness_share
    }
}

fn build_row(label: &'static str, ours: &[[f32; 3]], theirs: &[[f32; 3]]) -> Row {
    let n = ours.len().min(theirs.len());
    let mut de = Vec::with_capacity(n);
    let mut cstar = Vec::new();
    let mut c_all_ours = Vec::with_capacity(n);
    let mut c_all_theirs = Vec::with_capacity(n);
    let mut rg = Vec::new();
    let mut bg = Vec::new();
    let mut rg_c = Vec::new();
    let mut bg_c = Vec::new();
    let mut pairs = Vec::with_capacity(n);

    for (o, t) in ours.iter().zip(theirs.iter()) {
        de.push(deltae::delta_e_prophoto(*o, *t));
        let co = chroma(lab(*o));
        let ct = chroma(lab(*t));
        c_all_ours.push(co);
        c_all_theirs.push(ct);
        let colored = ct >= CSTAR_GATE;
        if colored {
            cstar.push(co / ct);
        }
        if let (Some(a), Some(b)) = (log_ratios(*o), log_ratios(*t)) {
            rg.push(a[0] - b[0]);
            bg.push(a[1] - b[1]);
            if colored {
                rg_c.push(a[0] - b[0]);
                bg_c.push(a[1] - b[1]);
            }
        }
        pairs.push((*o, *t));
    }

    let (de_median, de_p95, _max) = quantiles(&mut de);
    // 明度 / 色度的占比用合计口径：它是唯一可加、合计恰为 1 的那个（见 deltae 的文档）。
    let summary = deltae::summarize_prophoto(&pairs);
    let med_ours = median_of(&mut c_all_ours);
    let med_theirs = median_of(&mut c_all_theirs);

    Row {
        label,
        de_median,
        de_p95,
        cstar_ratio_all: if med_theirs > 0.0 { med_ours / med_theirs } else { f64::NAN },
        cstar_ratio: median_of(&mut cstar),
        cstar_n: cstar.len(),
        rg_dev: median_of(&mut rg),
        bg_dev: median_of(&mut bg),
        ratio_n: rg.len(),
        rg_dev_c: median_of(&mut rg_c),
        bg_dev_c: median_of(&mut bg_c),
        ratio_c_n: rg_c.len(),
        lightness_share: summary.total_share[0],
    }
}

/// **明度配对**子集上的统计：只取 `|ΔL*| ≤ 门槛` 的样本。
///
/// # 为什么还要这一步
///
/// 前面两列的比值差虽然对**曝光缩放**不变，却挡不住"我们的亮度本来就和参考不一样"：
/// 同一份色度在更暗的样本上算出的 C* 更小，而 R/G、B/G 又只在**同一亮度**上才等价于
/// 色度差。把亮度差排除之后，剩下的差就只能是色度本身——这是本诊断里最干净的一把尺子。
struct MatchedRow {
    label: &'static str,
    n: usize,
    cstar: f64,
    rg: f64,
    bg: f64,
}

impl MatchedRow {
    /// 两个比值差里较大的那个（绝对值，单位 log2）——「色度差有多大」就取它。
    fn worst_ratio(&self) -> f64 {
        self.rg.abs().max(self.bg.abs())
    }
}

fn matched_luminance(label: &'static str, ours: &[[f32; 3]], theirs: &[[f32; 3]]) -> MatchedRow {
    let mut n = 0usize;
    let mut cs = Vec::new();
    let mut rg = Vec::new();
    let mut bg = Vec::new();
    for (o, t) in ours.iter().zip(theirs.iter()) {
        let lo = lab(*o);
        let lt = lab(*t);
        if (lo[0] - lt[0]).abs() > LUM_MATCH {
            continue;
        }
        n += 1;
        let ct = chroma(lt);
        if ct >= CSTAR_GATE {
            cs.push(chroma(lo) / ct);
        }
        if let (Some(a), Some(b)) = (log_ratios(*o), log_ratios(*t)) {
            rg.push(a[0] - b[0]);
            bg.push(a[1] - b[1]);
        }
    }
    MatchedRow {
        label,
        n,
        cstar: median_of(&mut cs),
        rg: median_of(&mut rg),
        bg: median_of(&mut bg),
    }
}

// ---------------------------------------------------------------------------
// 按色相分段的细分
// ---------------------------------------------------------------------------

struct BandRow {
    name: &'static str,
    n: usize,
    /// ②的 C* 比、④的 C* 比。
    c2: f64,
    c4: f64,
    /// ②、④ 各自的 `Δlog2(R/G)` 与 `Δlog2(B/G)`。
    rg2: f64,
    rg4: f64,
    bg2: f64,
    bg4: f64,
    de2: f64,
    de4: f64,
}

fn build_bands(work: &[[f32; 3]], full: &[[f32; 3]], theirs: &[[f32; 3]]) -> Vec<BandRow> {
    let mut n = [0usize; 8];
    let mut c2: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut c4: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut rg2: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut rg4: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut bg2: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut bg4: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut de2: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());
    let mut de4: [Vec<f64>; 8] = std::array::from_fn(|_| Vec::new());

    for ((o, f), t) in work.iter().zip(full.iter()).zip(theirs.iter()) {
        let lt = lab(*t);
        let ct = chroma(lt);
        if ct < CSTAR_GATE {
            continue;
        }
        let b = band_index(hue_deg(lt));
        n[b] += 1;
        c2[b].push(chroma(lab(*o)) / ct);
        c4[b].push(chroma(lab(*f)) / ct);
        de2[b].push(deltae::delta_e_prophoto(*o, *t));
        de4[b].push(deltae::delta_e_prophoto(*f, *t));
        if let (Some(a), Some(rt)) = (log_ratios(*o), log_ratios(*t)) {
            rg2[b].push(a[0] - rt[0]);
            bg2[b].push(a[1] - rt[1]);
        }
        if let (Some(a), Some(rt)) = (log_ratios(*f), log_ratios(*t)) {
            rg4[b].push(a[0] - rt[0]);
            bg4[b].push(a[1] - rt[1]);
        }
    }

    (0..8)
        .map(|i| BandRow {
            name: BAND_NAMES[i],
            n: n[i],
            c2: median_of(&mut c2[i]),
            c4: median_of(&mut c4[i]),
            rg2: median_of(&mut rg2[i]),
            rg4: median_of(&mut rg4[i]),
            bg2: median_of(&mut bg2[i]),
            bg4: median_of(&mut bg4[i]),
            de2: median_of(&mut de2[i]),
            de4: median_of(&mut de4[i]),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 一张图的完整报告
// ---------------------------------------------------------------------------

struct ImageReport {
    stem: String,
    size: (usize, usize),
    rotated: bool,
    eval_count: usize,
    fit_count: usize,
    decode_secs: f64,
    wb: [f32; 3],
    curve_len: usize,
    lut_edge: usize,
    /// ①→② 的最小二乘矩阵与其拟合质量。
    m12: [[f32; 3]; 3],
    pair_n: usize,
    pair_rms: f64,
    pair_max: f64,
    m12_vs_documented: f32,
    /// ①→参考 的最小二乘矩阵（上界）与最优对角增益（上界）。
    m_oracle: [[f32; 3]; 3],
    gain: [f32; 3],
    /// 交叉验证（左半求解、右半评估）的两个上界，以及两半各自解出的 3×3。
    cross_matrix: Option<Row>,
    cross_gain: Option<Row>,
    cross_matched: Option<MatchedRow>,
    m_left: Option<[[f32; 3]; 3]>,
    m_right: Option<[[f32; 3]; 3]>,
    rows: Vec<Row>,
    bands: Vec<BandRow>,
    /// 明度配对子集上的各级统计。
    matched: Vec<MatchedRow>,
}

fn analyze(stem: &str) -> Option<ImageReport> {
    let fx = fitted_fixture_for(stem, true)?;
    assert!(fx.eval.len() > 100, "{stem}：评估样本只有 {} 个，统计无意义", fx.eval.len());

    let nef = samples_dir().join(format!("{stem}.NEF"));
    let opts = libraw::Options { demosaic: libraw::Demosaic::Dht, user_mul: None };
    let t0 = Instant::now();
    let cam = match libraw::decode_camera_linear(&nef, &opts) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("① 解码失败（{stem}）：{e}");
            return None;
        }
    };
    let decode_secs = t0.elapsed().as_secs_f64();

    let model = camera::read_model(&nef)?;
    let entry = camera::lookup(&model)?;
    let mg = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = mg.effective(cam.width, cam.height)?;
    // ① 与 ② 是同一个文件的两次解码，有效尺寸必须一致，否则逐级坐标会整体错位——
    // 那样得出的"分级表"比不对齐还危险：它会安静地给出错误结论。
    assert_eq!(
        (cw, ch),
        fx.size,
        "{stem}：① 与 ② 两次解码的有效尺寸不一致（{}×{} vs {}×{}）",
        cw,
        ch,
        fx.size.0,
        fx.size.1
    );

    let theirs: Vec<[f32; 3]> = fx.eval.iter().map(|s| s.theirs_linear).collect();
    let work: Vec<[f32; 3]> = fx.eval.iter().map(|s| s.ours_linear).collect();
    assert_eq!(theirs.len(), work.len());

    // ---------- ① 原始解码（相机空间） ----------
    let mut raw = Vec::with_capacity(work.len());
    let mut missing = 0usize;
    for s in &fx.eval {
        match cam.at(s.xy.0 + mg.left, s.xy.1 + mg.top) {
            Some(p) => raw.push(norm96(p)),
            None => {
                missing += 1;
                raw.push([f32::NAN; 3]);
            }
        }
    }
    assert_eq!(missing, 0, "{stem}：{missing} 个评估样本在①里取不到像素");

    // ---------- ③ 仅曲线、④ 曲线 + LUT ----------
    let t = &fx.transform;
    let curve: Vec<[f32; 3]> = work.iter().map(|v| fit::forward_curve(&t.curve, *v)).collect();
    let full: Vec<[f32; 3]> = work.iter().map(|v| t.apply(*v)).collect();
    let worst_full = full
        .iter()
        .zip(fx.eval.iter())
        .map(|(a, s)| max_abs(*a, s.ours_transformed))
        .fold(0f32, f32::max);
    assert!(
        worst_full < 1e-6,
        "{stem}：④ 的重算与脚手架的 ours_transformed 不符（最大差 {worst_full:e}）——\
         分级表里的④必须就是被评估的那一级，否则整张表都不可信"
    );

    // ---------- ①→②：是不是一个线性变换 ----------
    let usable = |v: &[f32; 3]| {
        v.iter()
            .all(|x| x.is_finite() && *x > LINEAR_GATE && *x < LINEAR_GATE_HI)
    };
    let p12: Vec<([f32; 3], [f32; 3])> = raw
        .iter()
        .zip(work.iter())
        .filter(|(a, b)| usable(a) && usable(b))
        .map(|(a, b)| (*a, *b))
        .collect();
    let d12 = color::derive_matrix(&p12).expect("①→② 的样本退化，无法求解矩阵");
    assert!(d12.samples > 1000, "{stem}：①→② 可用样本只有 {}", d12.samples);
    assert!(
        d12.rms < 1e-3,
        "{stem}：①→② 不是线性变换（rms {:e}）——那么「施加矩阵」与「解码」就分不成两级",
        d12.rms
    );

    // ---------- 上界一：理想的线性色彩转换（最优 3×3，同集求解） ----------
    let p1t: Vec<([f32; 3], [f32; 3])> = raw
        .iter()
        .zip(theirs.iter())
        .filter(|(a, b)| usable(a) && usable(b))
        .map(|(a, b)| (*a, *b))
        .collect();
    let oracle = color::derive_matrix(&p1t).expect("①→参考 的样本退化，无法求解矩阵");
    let oracle_matrix: Vec<[f32; 3]> =
        raw.iter().map(|v| mat3::mul_vec(oracle.matrix, *v)).collect();

    // ---------- 上界二：只允许逐通道增益（白平衡型修正，同集求解） ----------
    let gain = best_gain(&p1t);
    let gain_applied: Vec<[f32; 3]> = raw
        .iter()
        .map(|v| {
            [
                (v[0] as f64 * gain[0]) as f32,
                (v[1] as f64 * gain[1]) as f32,
                (v[2] as f64 * gain[2]) as f32,
            ]
        })
        .collect();

    // ---------- 交叉验证：左半求解、右半评估 ----------
    // 「同集求解」的上界有可能只是把这一张画面的色彩分布背了下来。因此再做一次空间留出：
    // 用画面左半的样本求变换、去评估右半，反之亦然。两半拍的是不同内容，交叉结果仍然好，
    // 才说明这个色度差能被一个**全局**色彩变换解释，而不是靠记住了这批像素。
    let (mut left, mut right): (Vec<_>, Vec<_>) = (Vec::new(), Vec::new());
    for ((a, b), s) in raw.iter().zip(theirs.iter()).zip(fx.eval.iter()) {
        if !usable(a) || !usable(b) {
            continue;
        }
        if s.xy.0 < cw / 2 {
            left.push((*a, *b));
        } else {
            right.push((*a, *b));
        }
    }
    let cross = if left.len() > 100 && right.len() > 100 {
        let m_left = color::derive_matrix(&left).expect("左半样本退化");
        let m_right = color::derive_matrix(&right).expect("右半样本退化");
        let g_left = best_gain(&left);
        let g_right = best_gain(&right);
        let cross_matrix: Vec<[f32; 3]> = raw
            .iter()
            .zip(fx.eval.iter())
            .map(|(v, s)| {
                let m = if s.xy.0 < cw / 2 { &m_right } else { &m_left };
                mat3::mul_vec(m.matrix, *v)
            })
            .collect();
        let cross_gain: Vec<[f32; 3]> = raw
            .iter()
            .zip(fx.eval.iter())
            .map(|(v, s)| {
                let g = if s.xy.0 < cw / 2 { &g_right } else { &g_left };
                [
                    (v[0] as f64 * g[0]) as f32,
                    (v[1] as f64 * g[1]) as f32,
                    (v[2] as f64 * g[2]) as f32,
                ]
            })
            .collect();
        Some((
            build_row("交叉：最优 3×3", &cross_matrix, &theirs),
            build_row("交叉：对角增益", &cross_gain, &theirs),
            m_left.matrix,
            m_right.matrix,
            matched_luminance("交叉：最优 3×3", &cross_matrix, &theirs),
        ))
    } else {
        None
    };

    let rows = vec![
        build_row("① 原始解码（相机空间）", &raw, &theirs),
        build_row("② 施加相机矩阵后", &work, &theirs),
        build_row("③ 仅曲线之后", &curve, &theirs),
        build_row("④ 曲线 + LUT 之后", &full, &theirs),
        build_row("⑤【上界】①+最优 3×3", &oracle_matrix, &theirs),
        build_row("⑥【上界】①+最优对角增益", &gain_applied, &theirs),
    ];
    debug_assert_eq!(rows.len(), 6);

    let (cross_matrix, cross_gain, m_left, m_right, cross_matched) = match cross {
        Some((a, b, c, d, e)) => (Some(a), Some(b), Some(c), Some(d), Some(e)),
        None => (None, None, None, None, None),
    };

    let matched = vec![
        matched_luminance("① 原始解码（相机空间）", &raw, &theirs),
        matched_luminance("② 施加相机矩阵后", &work, &theirs),
        matched_luminance("③ 仅曲线之后", &curve, &theirs),
        matched_luminance("④ 曲线 + LUT 之后", &full, &theirs),
        matched_luminance("⑤【上界】①+最优 3×3", &oracle_matrix, &theirs),
    ];

    Some(ImageReport {
        stem: stem.to_string(),
        size: fx.size,
        rotated: cam.rotated,
        eval_count: fx.eval.len(),
        fit_count: fx.fit_count,
        decode_secs,
        wb: cam.wb,
        curve_len: t.curve.len(),
        lut_edge: t.lut_edge,
        m12: d12.matrix,
        pair_n: d12.samples,
        pair_rms: d12.rms,
        pair_max: d12.max_abs,
        m12_vs_documented: mat_max_diff(&d12.matrix, &DOCUMENTED_MATRIX),
        m_oracle: oracle.matrix,
        gain: [gain[0] as f32, gain[1] as f32, gain[2] as f32],
        cross_matrix,
        cross_gain,
        cross_matched,
        m_left,
        m_right,
        rows,
        bands: build_bands(&work, &full, &theirs),
        matched,
    })
}

// ---------------------------------------------------------------------------
// 打印
// ---------------------------------------------------------------------------

fn print_table(rows: &[Row]) {
    let head = [
        pad("级别", 26),
        pad("ΔE00中位", 9),
        pad("P95", 8),
        pad("C*比全", 8),
        pad("C*比彩", 8),
        pad("彩样本", 8),
        pad("比样本", 8),
        pad("Δlog2R/G", 10),
        pad("Δlog2B/G", 10),
        pad("明度%", 7),
        pad("色度%", 7),
    ];
    eprintln!("\n{}", head.join(" "));
    for r in rows {
        let cells = [
            pad(r.label, 26),
            pad(&num(r.de_median, 3), 9),
            pad(&num(r.de_p95, 3), 8),
            pad(&num(r.cstar_ratio_all, 3), 8),
            pad(&num(r.cstar_ratio, 3), 8),
            pad(&r.cstar_n.to_string(), 8),
            pad(&r.ratio_n.to_string(), 8),
            pad(&signed(r.rg_dev, 3), 10),
            pad(&signed(r.bg_dev, 3), 10),
            pad(&num(100.0 * r.lightness_share, 1), 7),
            pad(&num(100.0 * r.chroma_share(), 1), 7),
        ];
        eprintln!("{}", cells.join(" "));
    }
    eprintln!(
        "  口径：C*比全 = 中位 C*(我们) / 中位 C*(参考)，全部样本；\
         C*比彩 = 逐样本比的中位，只统计参考 C* ≥ {CSTAR_GATE:.0} 的样本（避开小分母）。"
    );
    eprintln!(
        "  Δlog2(R/G)、Δlog2(B/G) = log2(我们的比值 / 参考的比值) 的中位，只统计三通道皆正的样本；\
         它是**尺度不变量**：任何逐像素曝光缩放都不改变它。"
    );
    eprintln!("  明度% / 色度% 为**全部误差量的合计占比**（可加、合计 100%）。");
}

fn print_bands(bands: &[BandRow]) {
    eprintln!("\n【按参考色相分段（八段中心与 render.rs 的 BAND_CENTERS 同一套）】");
    let head = [
        pad("段", 8),
        pad("样本", 8),
        pad("②C*比", 8),
        pad("④C*比", 8),
        pad("②Δlog2R/G", 11),
        pad("④Δlog2R/G", 11),
        pad("②Δlog2B/G", 11),
        pad("④Δlog2B/G", 11),
        pad("②ΔE00", 8),
        pad("④ΔE00", 8),
    ];
    eprintln!("{}", head.join(" "));
    for b in bands {
        if b.n == 0 {
            continue;
        }
        let cells = [
            pad(b.name, 8),
            pad(&b.n.to_string(), 8),
            pad(&num(b.c2, 3), 8),
            pad(&num(b.c4, 3), 8),
            pad(&signed(b.rg2, 3), 11),
            pad(&signed(b.rg4, 3), 11),
            pad(&signed(b.bg2, 3), 11),
            pad(&signed(b.bg4, 3), 11),
            pad(&num(b.de2, 3), 8),
            pad(&num(b.de4, 3), 8),
        ];
        eprintln!("{}", cells.join(" "));
    }
    eprintln!(
        "  分段依据是**参考**的色相角；C* 比受明度影响（同色度在更亮处 C* 更大），\
         所以只有 Δlog2(R/G)、Δlog2(B/G) 才是无亮度污染的色度指标。"
    );
}

fn print_matched(rows: &[MatchedRow], cross: Option<&MatchedRow>) {
    eprintln!(
        "\n【明度配对子集（|ΔL*| ≤ {LUM_MATCH:.0}）：亮度差被排除，剩下的只能是色度】"
    );
    let head = [
        pad("级别", 26),
        pad("样本", 10),
        pad("C*比", 8),
        pad("Δlog2R/G", 11),
        pad("Δlog2B/G", 11),
    ];
    eprintln!("{}", head.join(" "));
    for m in rows.iter().chain(cross.iter().copied()) {
        if m.n == 0 {
            continue;
        }
        let cells = [
            pad(m.label, 26),
            pad(&m.n.to_string(), 10),
            pad(&num(m.cstar, 3), 8),
            pad(&signed(m.rg, 3), 11),
            pad(&signed(m.bg, 3), 11),
        ];
        eprintln!("{}", cells.join(" "));
    }
    eprintln!(
        "  这一表里只有两个**比值**列是干净的尺子：C* 比即使在 |ΔL*| ≤ 2 的窗口内也仍有残余的亮度污染——\
         暗部（L* < 8）C* 近似正比于 Y，而那时 ±2 个 L* 对应的相对亮度差可以很大。判定一律用比值列。"
    );
}

fn print_verdict(r: &ImageReport) {
    let cam = &r.rows[L_CAM];
    let work = &r.rows[L_WORK];
    let curve = &r.rows[L_CURVE];
    let full = &r.rows[L_FULL];
    let om = &r.rows[L_ORACLE_MATRIX];
    let og = &r.rows[L_ORACLE_GAIN];
    let m_cam = &r.matched[0];
    let m_work = &r.matched[1];
    let m_curve = &r.matched[2];
    let m_full = &r.matched[3];

    // 判定一律用**明度配对**的有色样本上的比值差：它是尺度不变量，又不含亮度污染。
    let pre = m_work.worst_ratio();
    let post = m_full.worst_ratio();
    let chroma_word = if post <= pre * 0.3 {
        "把色度差基本修掉了"
    } else if post <= pre * 0.7 {
        "把色度差压小了（修，但没修完）"
    } else if post >= pre * 1.3 {
        "**把色度差放大了**"
    } else {
        "基本没改变色度差（既没引入也没修复）"
    };

    eprintln!("\n【判读（{}）】", r.stem);
    eprintln!(
        " · ΔE00 中位：①(相机空间) {:.3} → ②(矩阵后) {:.3} → ③(仅曲线) {:.3} → ④(曲线+LUT) {:.3}（目标 ≤ {:.1}）",
        cam.de_median, work.de_median, curve.de_median, full.de_median, TARGET_MEDIAN
    );
    eprintln!(
        " · 全样本 C* 比（中位之比，与「原始解码只有参考的 0.53~0.63 倍」同口径）：\
         ① {:.3} → ② {:.3} → ③ {:.3} → ④ {:.3}",
        cam.cstar_ratio_all, work.cstar_ratio_all, curve.cstar_ratio_all, full.cstar_ratio_all
    );
    eprintln!(
        "   注：① 行的 Lab 是把相机 RGB 直接按 ProPhoto 基解读的，跨基色的绝对值没有色度含义，\
         只用于看矩阵这一步改了什么。"
    );
    eprintln!(
        " · 有色样本（参考 C* ≥ {CSTAR_GATE:.0}，{} 个）上的比值差 Δlog2(R/G)：② {:+.3} → ③ {:+.3} → ④ {:+.3}；\
         Δlog2(B/G)：② {:+.3} → ③ {:+.3} → ④ {:+.3}",
        work.ratio_c_n,
        work.rg_dev_c,
        curve.rg_dev_c,
        full.rg_dev_c,
        work.bg_dev_c,
        curve.bg_dev_c,
        full.bg_dev_c
    );
    eprintln!(
        " · 明度配对子集上的比值差（最干净的尺子）Δlog2(R/G)：① {:+.3} → ② {:+.3} → ③ {:+.3} → ④ {:+.3}；\
         Δlog2(B/G)：① {:+.3} → ② {:+.3} → ③ {:+.3} → ④ {:+.3}",
        m_cam.rg, m_work.rg, m_curve.rg, m_full.rg, m_cam.bg, m_work.bg, m_curve.bg, m_full.bg
    );
    eprintln!(
        " · 拟合前（②）的色度差 {:.3}（log2），拟合后（④）{:.3} ⇒ ③④ {}",
        pre, post, chroma_word
    );
    eprintln!(
        " · 【乐观上界｜同集求解】①+最优 3×3：ΔE00 中位 {:.3}，C* 比 {:.3}；\
         ①+最优对角增益：ΔE00 中位 {:.3}，C* 比 {:.3}",
        om.de_median, om.cstar_ratio, og.de_median, og.cstar_ratio
    );
    match (&r.cross_matrix, &r.cross_gain, &r.cross_matched) {
        (Some(cm), Some(cg), Some(cx)) => eprintln!(
            " · 【乐观上界｜左半求解→右半评估】最优 3×3：ΔE00 中位 {:.3}，C* 比 {:.3}，\
             明度配对的 Δlog2(R/G) {:+.3}、Δlog2(B/G) {:+.3}；对角增益：ΔE00 中位 {:.3}，C* 比 {:.3}",
            cm.de_median, cm.cstar_ratio, cx.rg, cx.bg, cg.de_median, cg.cstar_ratio
        ),
        _ => eprintln!(" · 交叉验证样本不足，未做左右半留出"),
    }
    eprintln!(
        " · 误差量归因（合计口径）：② 明度 {:.1}% / 色度 {:.1}%；④ 明度 {:.1}% / 色度 {:.1}%",
        100.0 * work.lightness_share,
        100.0 * work.chroma_share(),
        100.0 * full.lightness_share,
        100.0 * full.chroma_share()
    );

    // ---------- 判定 ----------
    let preexisting = pre >= PREEXISTING_GATE;
    // 「线性色彩转换能不能解释这个色度差」以**留出**的交叉结果为准，同集上界只作参考。
    let cross_worst = r.cross_matched.as_ref().map(|c| c.worst_ratio()).unwrap_or(f64::NAN);
    let cross_de = r.cross_matrix.as_ref().map(|c| c.de_median).unwrap_or(f64::NAN);
    let matrix_can = cross_worst.is_finite()
        && cross_worst <= 0.5 * pre
        && cross_de <= 0.9 * work.de_median;
    let origin = if preexisting {
        if matrix_can {
            "色度差主要在①②之间（色彩转换环节）——一个全局线性色彩变换就能把色度差基本消掉，\
             说明①的色度信息是够的，差出在①→②所用的那个变换上"
        } else {
            "色度差主要在①（原始解码本身）或其上游——连留出的全局线性色彩变换都消不掉"
        }
    } else if post >= PREEXISTING_GATE {
        "色度差主要在③④（拟合环节把它引入或放大了）"
    } else {
        "这张图上没有可归因的色度差：两级都只是轻微偏离"
    };
    eprintln!(" 判定：{origin}");

    if post >= PREEXISTING_GATE || pre >= PREEXISTING_GATE {
        if pre >= PREEXISTING_GATE && post <= pre * 0.3 {
            eprintln!(
                " 补充：③④ 把色度差从 {pre:.3} 压到 {post:.3}（{:.0}% 被修掉）⇒ 拟合环节是**修复者**，\
                 不是引入者。",
                100.0 * (1.0 - post / pre)
            );
        } else if post >= pre * 1.3 {
            eprintln!(" 补充：③④ 不但没修，还把色度差从 {pre:.3} 放大到 {post:.3}。");
        } else {
            eprintln!(" 补充：③④ 对色度差几乎没有影响（{pre:.3} → {post:.3}）。");
        }
    } else {
        eprintln!(" 补充：前后都在门槛以下 ⇒ 这张图上没有可归因的色度问题。");
    }

    // ② 的 C* 偏低常被当成"色度不够"，但只加曲线就能抬回来 ⇒ 那是亮度造成的，不是色度缺失。
    if work.cstar_ratio_all < 0.75 && curve.cstar_ratio_all >= 1.25 * work.cstar_ratio_all {
        eprintln!(
            " 关键澄清：只加曲线（③）就把全样本 C* 比从 {:.3} 抬到 {:.3} ⇒ ② 的低 C* 主要是\
             「还没施加色调曲线、画面整体偏暗」（② 的明度项占 {:.1}%）造成的，**不是色度信息缺失**。\
             这正是不能拿 C* 判色度的原因。",
            work.cstar_ratio_all,
            curve.cstar_ratio_all,
            100.0 * work.lightness_share
        );
    }

    if full.lightness_share >= 0.5 {
        eprintln!(
            " 注意：④ 的误差量有 {:.1}% 是**明度**项 ⇒ 当前的 ΔE00 {:.3} 已经不是色度主导，\
             继续在色度上加力收益有限。",
            100.0 * full.lightness_share,
            full.de_median
        );
    } else {
        eprintln!(
            " 注意：④ 的误差量有 {:.1}% 在色度侧 ⇒ 残留误差仍是色度主导。",
            100.0 * full.chroma_share()
        );
    }
}

fn print_report(r: &ImageReport) {
    eprintln!("\n================ {} ================", r.stem);
    eprintln!(
        "有效尺寸 {}×{}（旋转 {}），拟合样本 {}，留出评估样本 {}，① 解码耗时 {:.1} s",
        r.size.0,
        r.size.1,
        if r.rotated { "是" } else { "否" },
        r.fit_count,
        r.eval_count,
        r.decode_secs
    );
    eprintln!(
        "白平衡系数 R/G/B = {:.4} / {:.4} / {:.4}；基准变换：曲线 {} 点，LUT 边长 {}",
        r.wb[0], r.wb[1], r.wb[2], r.curve_len, r.lut_edge
    );
    eprintln!(
        "①→② 最小二乘矩阵（{} 对未裁切样本，rms {:.2e}，最大残差 {:.2e}）：",
        r.pair_n, r.pair_rms, r.pair_max
    );
    eprintln!("   {}", fmt_mat(&r.m12));
    eprintln!(
        "   ①→② 是精确的线性变换，所以「施加矩阵」这一步本身不会凭空造出或丢掉色度；\
         该矩阵与任务 3.3 记录值的最大分量差 {:.2e}（同一机型级标定，差异来自样本群不同）。",
        r.m12_vs_documented
    );
    eprintln!(
        "①→参考 的最优 3×3（同集求解，乐观上界）：{}",
        fmt_mat(&r.m_oracle)
    );
    eprintln!(
        "①→参考 的最优对角增益（同集求解，乐观上界）：{:.4} / {:.4} / {:.4}",
        r.gain[0], r.gain[1], r.gain[2]
    );
    if let (Some(ml), Some(mr)) = (&r.m_left, &r.m_right) {
        eprintln!("交叉验证用的两半矩阵——左半：{}", fmt_mat(ml));
        eprintln!("                     右半：{}", fmt_mat(mr));
    }

    print_table(&r.rows);
    print_matched(&r.matched, r.cross_matched.as_ref());
    print_bands(&r.bands);
    print_verdict(r);
}

fn print_global(reports: &[ImageReport]) {
    eprintln!("\n================ 三图汇总 ================");
    let head = [
        pad("图", 12),
        pad("②色度差", 10),
        pad("④色度差", 10),
        pad("②C*比全", 10),
        pad("④C*比全", 10),
        pad("②ΔE00", 9),
        pad("④ΔE00", 9),
        pad("交叉3×3 ΔE00", 14),
        pad("交叉3×3 色度差", 16),
    ];
    eprintln!("{}", head.join(" "));
    eprintln!("（色度差 = 明度配对的有色样本上 |Δlog2(R/G)| 与 |Δlog2(B/G)| 的较大者，单位 log2）");
    for r in reports {
        let cm = r.cross_matrix.as_ref();
        let cx = r.cross_matched.as_ref();
        let cells = [
            pad(&r.stem, 12),
            pad(&num(r.matched[L_WORK].worst_ratio(), 3), 10),
            pad(&num(r.matched[L_FULL].worst_ratio(), 3), 10),
            pad(&num(r.rows[L_WORK].cstar_ratio_all, 3), 10),
            pad(&num(r.rows[L_FULL].cstar_ratio_all, 3), 10),
            pad(&num(r.rows[L_WORK].de_median, 3), 9),
            pad(&num(r.rows[L_FULL].de_median, 3), 9),
            pad(&cm.map(|c| num(c.de_median, 3)).unwrap_or_else(|| String::from("-")), 14),
            pad(&cx.map(|c| num(c.worst_ratio(), 3)).unwrap_or_else(|| String::from("-")), 16),
        ];
        eprintln!("{}", cells.join(" "));
    }

    if reports.len() >= 2 {
        let mut worst = 0f32;
        for w in reports.windows(2) {
            worst = worst.max(mat_max_diff(&w[0].m12, &w[1].m12));
        }
        eprintln!(
            "①→② 矩阵的跨图最大分量差 {worst:.2e}（对 1.0 而言约 {:.2}%）：远大于 16 位量化底噪 1.5e-05，\
             说明解码层的**有效**转换随图略变；但量级很小，不影响「矩阵按机型固定」这个近似。",
            100.0 * worst as f64
        );
    }

    let pre = reports
        .iter()
        .filter(|r| r.matched[L_WORK].worst_ratio() >= PREEXISTING_GATE)
        .count();
    let fixed_by_fit = reports
        .iter()
        .filter(|r| {
            let d2 = r.matched[L_WORK].worst_ratio();
            let d4 = r.matched[L_FULL].worst_ratio();
            d2 >= PREEXISTING_GATE && d4 <= d2 * 0.3
        })
        .count();
    let worsened = reports
        .iter()
        .filter(|r| {
            let d2 = r.matched[L_WORK].worst_ratio();
            let d4 = r.matched[L_FULL].worst_ratio();
            d4 >= PREEXISTING_GATE && d4 >= d2 * 1.3
        })
        .count();
    eprintln!(
        "\n判定：{}/{} 张图在②（施加相机矩阵之后、任何拟合之前）就有 ≥{PREEXISTING_GATE:.2} log2 的\
         色度差（明度配对、有色样本口径）；③④ 把它压到三成以下的有 {fixed_by_fit} 张，\
         反而放大的有 {worsened} 张。",
        pre,
        reports.len()
    );
    if pre == 0 {
        eprintln!("⇒ 色度差不在②之前：它在③④（拟合环节）被引入或放大。");
    } else if worsened == 0 {
        eprintln!(
            "⇒ 出现过的色度差**全部产生在①②（解码 / 相机矩阵色彩转换）**：它在任何拟合之前就已存在，\
             ③④ 只修不造（{}张被压到三成以下，{}张本来就在门槛以下）。",
            fixed_by_fit,
            reports.len() - pre
        );
    } else {
        eprintln!(
            "⇒ 色度差在②之前就已存在（属①②），且至少有一张在③④被进一步放大——两级都有问题。"
        );
    }
}

// ---------------------------------------------------------------------------
// 测试
// ---------------------------------------------------------------------------

#[test]
fn chroma_difference_by_pipeline_stage() {
    let stems: Vec<String> = match std::env::var("NRV_DIAG_STEM") {
        Ok(s) => vec![s],
        Err(_) => vec!["DSC_0001".into(), "DSC_0141".into(), "DSC_8562".into()],
    };

    let mut reports = Vec::new();
    for stem in &stems {
        // analyze 内部对缺席样本返回 None（脚手架已打印提示）；这里不用失败来惩罚缺席。
        if let Some(r) = analyze(stem) {
            print_report(&r);
            reports.push(r);
        }
    }

    if reports.is_empty() {
        eprintln!("三张样本都缺席，分级表未产出（这不是失败）");
        return;
    }
    print_global(&reports);
}
