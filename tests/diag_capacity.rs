//! 诊断 Q：**模型容量扫描（过拟合假设的判定）**。
//!
//! # 要回答的问题
//!
//! 所有「加表达力」的尝试都已证伪（逐通道曲线跨图 −1.2%、输入端增益只占 penalty ~11%、
//! 17³→33³ 更差、覆盖已排除、空间配准已排除）。「加容量更差」这个信号指向**过拟合**：
//! 参数 ≈ 曲线 129×3 + LUT 4913×3 ≈ 15k，而像素在空间上高度相关，有效独立样本远少于 8 万。
//!
//! 本诊断在**留一图**协议下扫两个容量轴，画一条「容量 vs 留出误差」的曲线：
//!
//! - **LUT 分辨率**：33 / 25 / 17 / 13 / 9 / 7 / 5（每轴采样数）
//! - **曲线 bins**：256 / 128 / 64 / 32 / 16
//! - 以及两者的**全交叉**（7×5 = 35 个配置）
//!
//! 每个配置报两个**留出**口径：
//!
//! - **跨图留出**：用其余全部图的**左半**池化拟合 → 在留出图的**右半**评估（生产口径）
//! - **自拟合留出**：用该图**自己**的**左半**拟合 → 在**右半**评估（该容量的上限）
//! - 泛化间隙 = 跨图留出 − 自拟合留出
//!
//! **不报训练误差**（诊断 N 的教训：训练口径的数字会被误读）。
//!
//! # 判定的四个分支（先定后看）
//!
//! | 跨图留出误差的形状 | 判定 |
//! |---|---|
//! | 内部有极小点（U 形） | **过拟合成立**，极小点就是该用的容量 |
//! | 单调下降（容量越大越好） | **过拟合被否** |
//! | 单调上升（容量越小越好） | 过拟合方向，但**最优点在扫描范围之外**（不能叫 U 形） |
//! | 变化小于 [`TREND_EPS`] | **没有观察到明显最优点**（如实报） |
//!
//! # 样本只解码一次
//!
//! 9 张图解码一次，全部 35 个配置复用同一批样本（左半拟合、右半评估），否则跑不完。
//!
//! # 限制
//!
//! 9 张图、两类场景；留一图仍然**共享场景类别**，跨场景这一维未被验证。
//! `NRV_Q_STEMS=...` 可只跑指定图（迭代用）。

mod common;

use common::samples_dir;
use nikonrawview::camera;
use nikonrawview::deltae;
use nikonrawview::fit::{self, Sample};
use nikonrawview::icc;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::mat3::{self, Mat3};

// ---------------------------------------------------------------------------
// 口径与门槛（先定后看）
// ---------------------------------------------------------------------------

const DEFAULT_STEMS: [&str; 9] = [
    "DSC_0001", "DSC_0141", "DSC_0314", "DSC_0507", "DSC_0540", "DSC_0563", "DSC_0567", "DSC_7737",
    "DSC_8562",
];

/// LUT 每轴采样数（容量轴 1）。
const EDGES: [usize; 7] = [33, 25, 17, 13, 9, 7, 5];
/// 曲线箱数（容量轴 2）。
const BINS: [usize; 5] = [256, 128, 64, 32, 16];

const STEP: usize = 17;
const CAM_LO: f32 = 1e-4;
const CAM_HI: f32 = 0.999;
const REF_LO: f32 = 1e-5;
const REF_HI: f32 = 0.999;
const MIN_SAMPLES: usize = 2000;

/// 判定用的可辨趋势阈值（ΔE00）。小于它视为「没有可辨趋势」。
const TREND_EPS: f64 = 0.05;

// ---------------------------------------------------------------------------
// 基础工具
// ---------------------------------------------------------------------------

/// 一个样本：相机空间值、工作空间值（`M0·c`）、参考值。
type Triple = ([f32; 3], [f32; 3], [f32; 3]);

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

fn de_median<F: Fn(&Triple) -> [f32; 3]>(rows: &[Triple], f: F) -> f64 {
    let mut des: Vec<f64> = rows.iter().map(|r| deltae::delta_e_prophoto(f(r), r.2)).collect();
    median_of(&mut des)
}

/// 按生产公式把工作空间值变成 LUT 节点序号（与 `fit::lut_coord` 同一套映射）。
fn lut_index(v: [f32; 3], edge: usize) -> usize {
    let q = |x: f32| ((fit::encode_srgb(x) * (edge - 1) as f32).round() as usize).min(edge - 1);
    (q(v[0]) * edge + q(v[1])) * edge + q(v[2])
}

/// 拟合一个配置：共享色调曲线（`bins` 箱）+ LUT（`edge³`）。
fn build(samples: &[Sample], bins: usize, edge: usize) -> (Vec<f32>, Vec<f32>) {
    let curve = fit::fit_curve(samples, bins);
    let lut = fit::fit_lut(samples, &curve, edge);
    (curve, lut)
}

/// 施加：先曲线，后 LUT（与 `BaseTransform::apply` 同序）。
fn apply(curve: &[f32], lut: &[f32], edge: usize, v: [f32; 3]) -> [f32; 3] {
    fit::apply_lut(lut, edge, fit::forward_curve(curve, v))
}

// ---------------------------------------------------------------------------
// 逐图样本
// ---------------------------------------------------------------------------

struct ImageData {
    stem: String,
    left: Vec<Triple>,
    right: Vec<Triple>,
    /// 相机空间左半（用于估「空间去相关后的独立样本数」）。
    left_sites: usize,
    secs: f64,
}

fn collect(stem: &str, m0: Mat3) -> Option<ImageData> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        common::skip("诊断 Q", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }
    let t0 = std::time::Instant::now();
    let tif_data = std::fs::read(&tif).ok()?;
    let img = fit::read_rgb16(&tif_data).ok()?;
    let plan = icc::plan_for(&tif_data).ok()?;
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let cam = libraw::decode_with_output_for_test(&nef, &opts, OutputColor::Camera).ok()?;
    let mo = camera::read_model(&nef)?;
    let entry = camera::lookup(&mo)?;
    let mg = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = mg.effective(cam.width, cam.height)?;
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
    // 采样步长 STEP 下，左半的空间格点数 —— 作为「空间去相关后的独立样本数」的粗口径
    let left_sites = (cw / 2).div_ceil(STEP) * ch.div_ceil(STEP);
    Some(ImageData { stem: stem.to_string(), left, right, left_sites, secs: t0.elapsed().as_secs_f64() })
}

// ---------------------------------------------------------------------------
// 主测试
// ---------------------------------------------------------------------------

fn env_stems() -> Vec<String> {
    match std::env::var("NRV_Q_STEMS") {
        Ok(v) if !v.trim().is_empty() => v
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect(),
        _ => DEFAULT_STEMS.iter().map(|s| s.to_string()).collect(),
    }
}

/// 一个容量配置的结果。
#[derive(Clone, Copy, Default)]
struct Cell {
    /// 跨图留出：逐图 ΔE00 中位（用于算均值/中位）
    cross: [f64; 32],
    /// 自拟合留出：逐图 ΔE00 中位
    self_fit: [f64; 32],
    n: usize,
    params: usize,
    occupied: usize,
    cal_samples: usize,
}

impl Cell {
    fn cross_mean(&self) -> f64 {
        mean_of(&self.cross[..self.n])
    }
    fn cross_med(&self) -> f64 {
        median_of(&mut self.cross[..self.n].to_vec())
    }
    fn self_mean(&self) -> f64 {
        mean_of(&self.self_fit[..self.n])
    }
}

/// 判定一条「容量递增 vs 误差」曲线的形状。
fn shape_of(errs: &[f64]) -> &'static str {
    let n = errs.len();
    if n < 3 {
        return "样本不足";
    }
    let (mut amin, mut best) = (0usize, f64::INFINITY);
    for (i, e) in errs.iter().enumerate() {
        if *e < best {
            best = *e;
            amin = i;
        }
    }
    let lo = errs[0];
    let hi = errs[n - 1];
    let interior_min = amin > 0 && amin < n - 1 && errs[amin] <= lo - TREND_EPS && errs[amin] <= hi - TREND_EPS;
    if interior_min {
        return "U 形（内部极小）";
    }
    let mut up = true;
    let mut down = true;
    for w in errs.windows(2) {
        if w[1] < w[0] - TREND_EPS {
            up = false;
        }
        if w[1] > w[0] + TREND_EPS {
            down = false;
        }
    }
    match (up, down) {
        (true, false) => "单调上升（容量越小越好）",
        (false, true) => "单调下降（容量越大越好）",
        (true, true) => "无可辨趋势",
        (false, false) => "非单调、但极小点在端点",
    }
}

#[test]
fn capacity_sweep() {
    eprintln!("=== 诊断 Q：模型容量扫描（留出误差 vs 容量）===");
    let Some(cm) = camera::color_matrix("Nikon Z 8") else {
        eprintln!("跳过：机型表里没有 Z8 的色彩矩阵");
        return;
    };
    let m0 = cm.matrix;
    eprintln!("协议：留一图。跨图留出 = 其余图**左半**池化拟合 → 留出图**右半**评估；");
    eprintln!("      自拟合留出 = 该图**左半**拟合 → **右半**评估。**不报训练误差。**");
    eprintln!("判定门槛（先定后看）：趋势小于 {TREND_EPS} ΔE00 视为不可辨；内部极小需比两端都低 {TREND_EPS} 以上。");

    let stems = env_stems();
    eprintln!("\n[0] 图（{} 张）：{}", stems.len(), stems.join(" "));
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
    let n_img = images.len();
    eprintln!(
        "  解码总耗时 {:.0} s · 左半样本 {} · 右半样本 {}（每张中位）",
        images.iter().map(|i| i.secs).sum::<f64>(),
        median_of(&mut images.iter().map(|i| i.left.len() as f64).collect::<Vec<_>>()),
        median_of(&mut images.iter().map(|i| i.right.len() as f64).collect::<Vec<_>>())
    );

    // ---- 样本只解码一次，各容量复用 ----
    // 池化样本：对每个留出图 i，cal_i = 其余图的左半。这里只存每张图自己的 Sample，
    // 每次现拼（避免 35×9 份拷贝）。
    let per_image: Vec<Vec<Sample>> = images
        .iter()
        .map(|im| im.left.iter().map(|r| Sample { ours: r.1, theirs: r.2 }).collect())
        .collect();

    let mut grid: Vec<Vec<Cell>> = vec![vec![Cell::default(); BINS.len()]; EDGES.len()];
    eprintln!("\n[1] 扫描中：{} 个配置 × {} 张留出图 × 2 个口径", EDGES.len() * BINS.len(), n_img);
    for (ei, &edge) in EDGES.iter().enumerate() {
        for (bi, &bins) in BINS.iter().enumerate() {
            let mut cell = Cell { n: 0, ..Default::default() };
            for held in 0..n_img {
                // 跨图：其余图左半池化
                let mut cal: Vec<Sample> = Vec::new();
                for (k, s) in per_image.iter().enumerate() {
                    if k != held {
                        cal.extend_from_slice(s);
                    }
                }
                let (curve, lut) = build(&cal, bins, edge);
                let c = de_median(&images[held].right, |r| apply(&curve, &lut, edge, r.1));
                // 占用节点数（跨图 LUT 在留出图样本上的支撑）
                let n_nodes = edge * edge * edge;
                let mut hit = vec![false; n_nodes];
                for s in &cal {
                    hit[lut_index(fit::forward_curve(&curve, s.ours), edge)] = true;
                }
                let occupied = hit.iter().filter(|v| **v).count();

                // 自拟合：该图左半
                let (curve_s, lut_s) = build(&per_image[held], bins, edge);
                let s_de = de_median(&images[held].right, |r| apply(&curve_s, &lut_s, edge, r.1));

                cell.cross[cell.n] = c;
                cell.self_fit[cell.n] = s_de;
                cell.n += 1;
                cell.occupied = occupied;
                cell.cal_samples = cal.len();
            }
            cell.params = (bins + 1) + 3 * edge * edge * edge;
            grid[ei][bi] = cell;
            eprint!(".");
        }
    }
    eprintln!(" 完成");

    // 粗口径的「空间去相关后独立样本数」：左半格点数（步长 STEP）
    let sites_med = median_of(&mut images.iter().map(|i| i.left_sites as f64).collect::<Vec<_>>());

    // ---- [2] LUT 分辨率单独扫（bins = 128）----
    let bi128 = BINS.iter().position(|b| *b == 128).unwrap();
    eprintln!("\n[2] LUT 分辨率单独扫（曲线固定 {}-箱）", BINS[bi128]);
    eprintln!(
        "  {:>5} {:>12} {:>12} {:>12} {:>12} {:>10} {:>10} {:>12}",
        "LUT", "跨图留出(均值)", "跨图留出(中位)", "自拟合留出(均值)", "泛化间隙", "参数", "占用节点", "样本/参数"
    );
    for (ei, &edge) in EDGES.iter().enumerate() {
        let c = &grid[ei][bi128];
        eprintln!(
            "  {:>4}³ {:>12.3} {:>12.3} {:>12.3} {:>12.3} {:>10} {:>10} {:>12.1}",
            edge,
            c.cross_mean(),
            c.cross_med(),
            c.self_mean(),
            c.cross_mean() - c.self_mean(),
            c.params,
            c.occupied,
            c.cal_samples as f64 / c.params as f64
        );
    }
    let edge_cross: Vec<f64> = EDGES.iter().enumerate().map(|(ei, _)| grid[ei][bi128].cross_mean()).collect();
    eprintln!("  形状判定（容量递增 = LUT 由小到大）：{}", shape_of(&rev(&edge_cross)));

    // ---- [3] 曲线 bins 单独扫（LUT = 17³）----
    let ei17 = EDGES.iter().position(|e| *e == 17).unwrap();
    eprintln!("\n[3] 曲线 bins 单独扫（LUT 固定 {}³）", EDGES[ei17]);
    eprintln!(
        "  {:>5} {:>12} {:>12} {:>12} {:>12} {:>10} {:>10} {:>12}",
        "bins", "跨图留出(均值)", "跨图留出(中位)", "自拟合留出(均值)", "泛化间隙", "参数", "占用节点", "样本/参数"
    );
    for (bi, &bins) in BINS.iter().enumerate() {
        let c = &grid[ei17][bi];
        eprintln!(
            "  {:>5} {:>12.3} {:>12.3} {:>12.3} {:>12.3} {:>10} {:>10} {:>12.1}",
            bins,
            c.cross_mean(),
            c.cross_med(),
            c.self_mean(),
            c.cross_mean() - c.self_mean(),
            c.params,
            c.occupied,
            c.cal_samples as f64 / c.params as f64
        );
    }
    let bins_cross: Vec<f64> = BINS.iter().enumerate().map(|(bi, _)| grid[ei17][bi].cross_mean()).collect();
    eprintln!("  形状判定（容量递增 = bins 由小到大）：{}", shape_of(&rev(&bins_cross)));

    // ---- [4] 全交叉 ----
    eprintln!("\n[4] 全交叉（行 = LUT 每轴采样数，列 = 曲线 bins）");
    eprintln!("  4a) **跨图留出** ΔE00 中位均值");
    eprint!("  {:>5}", "LUT\\bins");
    for b in BINS {
        eprint!(" {:>8}", b);
    }
    eprintln!(" {:>9}", "行最差-最好");
    let mut all_cross: Vec<f64> = Vec::new();
    for (ei, &edge) in EDGES.iter().enumerate() {
        eprint!("  {:>5}", format!("{edge}³"));
        let mut row: Vec<f64> = Vec::new();
        for (bi, _) in BINS.iter().enumerate() {
            let v = grid[ei][bi].cross_mean();
            row.push(v);
            all_cross.push(v);
            eprint!(" {:>8.3}", v);
        }
        let lo = row.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = row.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        eprintln!(" {:>9.3}", hi - lo);
    }
    eprintln!("  4b) **自拟合留出** ΔE00 中位均值（该容量的上限，仍是留出半幅）");
    eprint!("  {:>5}", "LUT\\bins");
    for b in BINS {
        eprint!(" {:>8}", b);
    }
    eprintln!();
    let mut all_self: Vec<f64> = Vec::new();
    for (ei, &edge) in EDGES.iter().enumerate() {
        eprint!("  {:>5}", format!("{edge}³"));
        for (bi, _) in BINS.iter().enumerate() {
            let v = grid[ei][bi].self_mean();
            all_self.push(v);
            eprint!(" {:>8.3}", v);
        }
        eprintln!();
    }
    eprintln!("  4c) 泛化间隙 = 4a − 4b");
    eprint!("  {:>5}", "LUT\\bins");
    for b in BINS {
        eprint!(" {:>8}", b);
    }
    eprintln!();
    for (ei, &edge) in EDGES.iter().enumerate() {
        eprint!("  {:>5}", format!("{edge}³"));
        for (bi, _) in BINS.iter().enumerate() {
            eprint!(" {:>8.3}", grid[ei][bi].cross_mean() - grid[ei][bi].self_mean());
        }
        eprintln!();
    }

    // 全局最优（跨图留出）
    let mut best = (0usize, 0usize, f64::INFINITY);
    for (ei, _) in EDGES.iter().enumerate() {
        for (bi, _) in BINS.iter().enumerate() {
            let v = grid[ei][bi].cross_mean();
            if v < best.2 {
                best = (ei, bi, v);
            }
        }
    }
    let (bei, bbi, bv) = best;
    eprintln!(
        "\n  全表最优：LUT {}³ × 曲线 {}-箱 → 跨图留出 {:.3}（自拟合留出 {:.3}，间隙 {:.3}，参数 {}）",
        EDGES[bei],
        BINS[bbi],
        bv,
        grid[bei][bbi].self_mean(),
        bv - grid[bei][bbi].self_mean(),
        grid[bei][bbi].params
    );
    eprintln!(
        "  全表范围：跨图留出 {:.3} ~ {:.3}（极差 {:.3}）；自拟合留出 {:.3} ~ {:.3}（极差 {:.3}）",
        all_cross.iter().cloned().fold(f64::INFINITY, f64::min),
        all_cross.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        all_cross.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - all_cross.iter().cloned().fold(f64::INFINITY, f64::min),
        all_self.iter().cloned().fold(f64::INFINITY, f64::min),
        all_self.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        all_self.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - all_self.iter().cloned().fold(f64::INFINITY, f64::min)
    );

    eprintln!("  4d) 逐图跨图留出 ΔE00 中位（看逐图波动；三列 = 两端配置与全表最优）");
    let cols = [(0usize, bi128), (EDGES.len() - 1, BINS.len() - 1), (bei, bbi)];
    eprint!("  {:<10}", "图");
    for (ei, bi) in cols {
        eprint!(" {:>14}", format!("{}³×{}箱", EDGES[ei], BINS[bi]));
    }
    eprintln!();
    for (k, im) in images.iter().enumerate() {
        eprint!("  {:<10}", im.stem);
        for (ei, bi) in cols {
            eprint!(" {:>14.3}", grid[ei][bi].cross[k]);
        }
        eprintln!();
    }

    // ---- [5] 参数 / 有效样本 ----
    eprintln!("\n[5] 参数个数 / 有效样本数（口径写明）");
    eprintln!("  · 参数 = (曲线 bins + 1) × 1 条共享曲线 + 3 × LUT 边长³");
    eprintln!("    （曲线是**共享**的 1 条，不是 3 条；强制单调 + 空箱内插会进一步降低其有效自由度）");
    eprintln!("  · 名义样本数 N = 池化后的左半样本数（拟合时实际用到的像素数）");
    eprintln!("  · 占用节点 = 池化样本在 LUT 网格上实际命中的节点数（其余靠邻域填充）");
    eprintln!("  · 空间去相关后独立样本数的**粗口径**：左半格点数 = {}（步长 {}，即把每个采样格视为 1 个独立样本；", sites_med, STEP);
    eprintln!("    真实值介于「格点数」与「像素数」之间——同一格内相邻像素高度相关，所以格点数是**下界口径**）");
    eprintln!(
        "  {:>5} {:>5} {:>10} {:>12} {:>12} {:>14} {:>14}",
        "LUT", "bins", "参数", "名义 N", "N/参数", "占用节点/总数", "N/格点数"
    );
    for (ei, &edge) in EDGES.iter().enumerate() {
        for (bi, &bins) in BINS.iter().enumerate() {
            // 交叉表全打印会重复太多行，这里只打两条轴上的点 + 全局最优
            let interesting = bi == bi128 || ei == ei17 || (ei == bei && bi == bbi);
            if !interesting {
                continue;
            }
            let c = &grid[ei][bi];
            let nodes = edge * edge * edge;
            eprintln!(
                "  {:>4}³ {:>5} {:>10} {:>12} {:>12.1} {:>14} {:>14.1}",
                edge,
                bins,
                c.params,
                c.cal_samples,
                c.cal_samples as f64 / c.params as f64,
                format!("{}/{}", c.occupied, nodes),
                c.cal_samples as f64 / sites_med
            );
        }
    }

    // ---- [6] 判定 ----
    let edge_asc = rev(&edge_cross);
    let bins_asc = rev(&bins_cross);
    let s_edge = shape_of(&edge_asc);
    let s_bins = shape_of(&bins_asc);
    eprintln!("\n[6] 判定（核心产出）");
    eprintln!("  1) LUT 分辨率轴（5³→33³）：跨图留出 {} ⇒ **{s_edge}**", fmt_seq(&edge_asc));
    eprintln!("  2) 曲线 bins 轴（16→256）：跨图留出 {} ⇒ **{s_bins}**", fmt_seq(&bins_asc));
    let spread = all_cross.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - all_cross.iter().cloned().fold(f64::INFINITY, f64::min);
    eprintln!(
        "  3) 全交叉表里跨图留出的极差只有 {spread:.3} ΔE00（最优 {:.3} @ {}³×{}-箱）⇒ {}",
        bv,
        EDGES[bei],
        BINS[bbi],
        if spread < 3.0 * TREND_EPS {
            "**容量在这两个轴上几乎不影响跨图留出误差**"
        } else {
            "容量有可辨影响，按上面的形状判定走"
        }
    );
    let u_like = s_edge.starts_with("U 形") || s_bins.starts_with("U 形");
    let rising = s_edge.starts_with("单调上升") || s_bins.starts_with("单调上升");
    let falling = s_edge.starts_with("单调下降") || s_bins.starts_with("单调下降");
    if u_like {
        eprintln!("  4) ⇒ **过拟合假设成立**（存在内部极小点）；该用的容量见上。");
    } else if falling && !rising {
        eprintln!("  4) ⇒ **过拟合假设被否**（容量越大越好，没有内部极小点）：应回到「空间/色彩区域绑定」另找原因。");
    } else if rising && !falling {
        eprintln!("  4) ⇒ **过拟合方向**（容量越小越好），但**最优点在扫描范围之外**：");
        eprintln!("     扫描下界（{}³ / {}-箱）已经是最优，若继续降容量可能还有收益 —— 这**不是** U 形，不能据此给出最优容量。", EDGES[EDGES.len() - 1], BINS[BINS.len() - 1]);
    } else {
        eprintln!("  4) ⇒ **没有观察到明显最优点**：容量扫描在这两个轴上都不改变结论（极差 {spread:.3}）。");
    }
    eprintln!();
    eprintln!("  限制（必须随结论一起转述）：");
    eprintln!("  · {} 张图、两类场景；留一图仍然**共享场景类别**，跨场景这一维未被验证。", n_img);
    eprintln!("  · 全部是**留出**口径（跨图留出 / 自拟合留出），没有任何训练误差。");
    eprintln!("  · 容量只扫了「曲线 bins × LUT 边长」这两轴；池化方式、损失函数、样本加权都没动。");
    eprintln!("  · 「自拟合留出」是该容量的上限口径，它自己也会过拟合（见诊断 N），所以间隙的解释要小心。");
    eprintln!("  · LUT 空节点靠邻域填充，占用节点数只描述**直接支撑**，不等于有效自由度。");
}

fn rev(v: &[f64]) -> Vec<f64> {
    v.iter().rev().cloned().collect()
}

fn fmt_seq(v: &[f64]) -> String {
    v.iter().map(|x| format!("{x:.3}")).collect::<Vec<_>>().join(" → ")
}

/// 仪器自检：
/// 1. `lut_index` 与 `apply_lut` 的节点映射一致（与诊断 N 同一套防漂移检查）；
/// 2. **仪器能不能看见过拟合**——在合成的「平滑映射 + 噪声」样本上，LUT 分辨率扫描
///    应当出现「太粗欠拟合 / 太细被噪声主导」的形状（33³ 不优于中间容量）。
///    如果这条不成立，那么真实数据上「没看到 U 形」就说明不了任何事。
#[test]
fn instrument_checks() {
    // 1) 节点映射
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

    // 2) 仪器能不能看见过拟合：合成「平滑映射 + 噪声」
    //    样本量固定、噪声固定，只改 LUT 分辨率。
    let mut seed = 12345u64;
    let mut rnd = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) as f32 / (1u64 << 31) as f32) - 0.5
    };
    let truth = |v: [f32; 3]| [v[0].powf(1.4), v[1].powf(1.4), v[2].powf(1.4)];
    let mut samples: Vec<Sample> = Vec::new();
    for _ in 0..6000 {
        let ours = [
            (rnd() + 0.5).clamp(0.02, 0.98),
            (rnd() + 0.5).clamp(0.02, 0.98),
            (rnd() + 0.5).clamp(0.02, 0.98),
        ];
        let t = truth(ours);
        let nz = 0.02;
        samples.push(Sample {
            ours,
            theirs: [t[0] + nz * rnd(), t[1] + nz * rnd(), t[2] + nz * rnd()],
        });
    }
    // 干净的测试网格（16³）
    let mut test: Vec<Sample> = Vec::new();
    for i in 0..16 {
        for j in 0..16 {
            for k in 0..16 {
                let ours = [
                    (i as f32 + 0.5) / 16.0,
                    (j as f32 + 0.5) / 16.0,
                    (k as f32 + 0.5) / 16.0,
                ];
                test.push(Sample { ours, theirs: truth(ours) });
            }
        }
    }
    let mut errs: Vec<(usize, f64)> = Vec::new();
    for edge in [5usize, 9, 13, 17, 25] {
        let (curve, lut) = build(&samples, 32, edge);
        let mut s = 0f64;
        let mut n = 0usize;
        for t in &test {
            let p = apply(&curve, &lut, edge, t.ours);
            for (pk, tk) in p.iter().zip(t.theirs.iter()) {
                s += ((*pk - *tk) as f64).abs();
                n += 1;
            }
        }
        errs.push((edge, s / n as f64));
    }
    let best_edge = errs
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|x| x.0)
        .unwrap_or(0);
    eprintln!(
        "  [自检] 合成「平滑映射 + 噪声」上的 LUT 分辨率扫描（平均绝对误差）：{} ⇒ 最优 {}³",
        errs.iter().map(|(e, v)| format!("{e}³ {v:.4}")).collect::<Vec<_>>().join(" · "),
        best_edge
    );
    eprintln!("  [自检] ⇒ 仪器**能**看见过拟合（若非如此，真实数据上的阴性结论就没有意义）");
    assert!(
        best_edge != 25,
        "在有噪声的合成样本上，LUT 分辨率不应「越细越好」；实得 {:?}——若这条不成立，\
         仪器就看不见过拟合，真实数据上的阴性结论也就没有意义",
        errs
    );
    assert!(
        errs.iter().all(|(_, e)| e.is_finite()),
        "合成检查里出现了非有限误差：{errs:?}"
    );
}
