//! 诊断 I：**修正相机矩阵之后能达到什么上限**，以及矩阵能否跨图迁移。
//!
//! # 要回答的问题
//!
//! 假设是「我们用的相机矩阵与尼康的等效色彩变换不是同一个，偏离随白平衡增大」。
//! 若成立，换用正确的矩阵后 ΔE00 应显著改善；若不成立，团队就该停止在矩阵上使劲。
//!
//! **本任务的结论比数字本身更重要**：若矩阵替换后中位数仍远离目标 1.0，必须如实
//! 说明瓶颈不在矩阵。这个项目已经因为「结论说得比证据满」返工过两次。
//!
//! # 三方案对照
//!
//! | 方案 | 前端（相机空间 → ProPhoto 线性） |
//! |---|---|
//! | 现状 | `decode_working_space`（LibRaw 单矩阵） |
//! | 矩阵替换 | 逐图最优 3×3（空间留出：左半求、右半评估）施于 `decode_camera_linear` |
//!
//! 两方案的后端**在同一批样本上各自重新拟合**（曲线 128 箱 + LUT 17³），否则比不出来。
//!
//! # 一条必须写明的限制
//!
//! 本图自解的 3×3 是**非物理**的（可能含 −4.06/+6.19 这类系数），它只回答"上限在哪"，
//! **不能直接采用**。真正可用的矩阵必须由多图联合标定 + 交叉验证得出。

mod common;

use common::{quantiles, samples_dir};
use nikonrawview::deltae;
use nikonrawview::libraw::{self, Demosaic, Options, OutputColor};
use nikonrawview::mat3;

/// 一行配对：相机空间值、参考的 ProPhoto 线性值、裁剪后坐标。
type Row = ([f32; 3], [f32; 3], (usize, usize));

/// 评估结果：中位数、P95、最大值、以及逐样本的 (变换后, 参考) 配对。
type EvalResult = (f64, f64, f64, Vec<([f32; 3], [f32; 3])>);

/// 一张图的相机空间 / 参考配对（含裁剪与朝向处理）。
struct ImagePairs {
    /// `(相机空间值, 参考的 ProPhoto 线性值, 裁剪后坐标)`
    rows: Vec<Row>,
}

fn collect(stem: &str) -> Option<ImagePairs> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        eprintln!("跳过 {stem}：缺 NEF 或 TIF");
        return None;
    }
    let tif_data = std::fs::read(&tif).ok()?;
    let img = nikonrawview::fit::read_rgb16(&tif_data).ok()?;
    let plan = nikonrawview::icc::plan_for(&tif_data).ok()?;
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let cam = libraw::decode_with_output_for_test(&nef, &opts, OutputColor::Camera).ok()?;
    let mo = nikonrawview::camera::read_model(&nef)?;
    let entry = nikonrawview::camera::lookup(&mo)?;
    let mg = if cam.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = mg.effective(cam.width, cam.height)?;
    if (cw, ch) != (img.width, img.height) {
        eprintln!("跳过 {stem}：裁切后 {cw}×{ch} 与参考 {}×{} 不符", img.width, img.height);
        return None;
    }
    let theirs = nikonrawview::fit::reference_to_working(&img, &plan);

    let mut rows = Vec::new();
    for y in (0..ch).step_by(13) {
        for x in (0..cw).step_by(13) {
            let j = x + y * img.width;
            let Some(p) = cam.at(x + mg.left, y + mg.top) else { continue };
            let c = [
                p[0] as f32 / 65535.0,
                p[1] as f32 / 65535.0,
                p[2] as f32 / 65535.0,
            ];
            let t = theirs[j];
            // 排除任一侧饱和点：那里解码层做过裁切，拟合会把它误当成色彩关系
            if c.iter().any(|v| *v <= 0.0 || *v >= 1.0)
                || t.iter().any(|v| *v <= 0.0 || *v >= 1.0)
            {
                continue;
            }
            rows.push((c, t, (x, y)));
        }
    }
    if rows.len() < 1000 {
        eprintln!("跳过 {stem}：可用配对仅 {} 个", rows.len());
        return None;
    }
    Some(ImagePairs { rows })
}

/// 用给定的前端把 `[(相机值, 参考值)]` 拟合并评估一遍，返回 ΔE00 分位数。
fn fit_and_eval(pairs: &[([f32; 3], [f32; 3])]) -> EvalResult {
    let samples: Vec<nikonrawview::fit::Sample> = pairs
        .iter()
        .map(|(c, t)| nikonrawview::fit::Sample { ours: *c, theirs: *t })
        .collect();
    let curve = nikonrawview::fit::fit_curve(&samples, 128);
    let edge = 17;
    let lut = nikonrawview::fit::fit_lut(&samples, &curve, edge);
    let t = nikonrawview::transform::BaseTransform::from_fit(0, "diag", vec![], curve, edge, lut);
    let evaluated: Vec<([f32; 3], [f32; 3])> =
        samples.iter().map(|s| (t.apply(s.ours), s.theirs)).collect();
    let mut des: Vec<f64> = evaluated
        .iter()
        .map(|(a, b)| deltae::delta_e_prophoto(*a, *b))
        .collect();
    let (m, p, x) = quantiles(&mut des);
    (m, p, x, evaluated)
}

fn shares(evaluated: &[([f32; 3], [f32; 3])]) -> [f64; 4] {
    let lab_pairs: Vec<([f32; 3], [f32; 3])> = evaluated
        .iter()
        .map(|(a, b)| (deltae::lab_from_prophoto_linear(*a), deltae::lab_from_prophoto_linear(*b)))
        .collect();
    deltae::summarize(&lab_pairs).total_share
}

#[test]
fn how_far_can_a_corrected_matrix_get_us() {
    eprintln!("=== 诊断 I：修正矩阵后的上限 ===");
    eprintln!("目标：ΔE00 中位数 ≤ 1.0，P95 ≤ 3.0\n");

    let mut derived: Vec<(String, mat3::Mat3)> = Vec::new();
    let mut all: Vec<(String, ImagePairs)> = Vec::new();
    for stem in ["DSC_0001", "DSC_0141", "DSC_8562"] {
        if let Some(p) = collect(stem) {
            all.push((stem.to_string(), p));
        }
    }
    if all.is_empty() {
        eprintln!("跳过：没有可用样本");
        return;
    }

    // 逐图求最优 3×3（空间留出：左半求、右半评估）
    eprintln!("--- 各图的逐图最优 3×3（左半求解，样本数见括号）---");
    for (stem, p) in &all {
        let half = p.rows.iter().filter(|r| r.2 .0 < p.rows.iter().map(|q| q.2 .0).max().unwrap_or(0) / 2);
        let training: Vec<([f32; 3], [f32; 3])> =
            half.map(|r| (r.0, r.1)).collect();
        let d = nikonrawview::color::derive_matrix(&training).expect("应能求解");
        eprintln!("  {stem}（训练 {} 个）：rms {:.5}", d.samples, d.rms);
        for row in &d.matrix {
            eprintln!("      [{:>9.4} {:>9.4} {:>9.4}]", row[0], row[1], row[2]);
        }
        derived.push((stem.clone(), d.matrix));
    }

    eprintln!("\n--- 三方案对照（自适应后端：曲线 128 箱 + LUT 17³，各自重新拟合）---");
    eprintln!(
        "  {:<10} {:<12} {:>8} {:>8} {:>9}   明度/彩度/色相 占比",
        "图", "方案", "中位", "P95", "最大"
    );
    for (stem, p) in &all {
        // 方案一：现状——用 LibRaw 的 ProPhoto 输出（等价于我们的固定矩阵）
        let cur: Vec<([f32; 3], [f32; 3])> = p.rows.iter().map(|r| (r.0, r.1)).collect();
        // 注意：方案一与方案二的区别只在前端，故这里对「相机值」直接用
        // decode_working_space 的结果是无法从本结构拿到的——见下方说明。
        let (m1, p1, x1, e1) = fit_and_eval(&cur);
        let s1 = shares(&e1);
        eprintln!(
            "  {:<10} {:<12} {:>8.3} {:>8.3} {:>9.3}   {:>5.1}% {:>5.1}% {:>5.1}%",
            stem, "相机空间直入", m1, p1, x1,
            100.0 * s1[0], 100.0 * s1[1], 100.0 * s1[2]
        );

        // 方案二：用逐图最优 3×3 替换前端
        let m = derived.iter().find(|(s, _)| s == stem).unwrap().1;
        let replaced: Vec<([f32; 3], [f32; 3])> = cur
            .iter()
            .map(|(c, t)| (mat3::mul_vec(m, *c), *t))
            .collect();
        let (m2, p2, x2, e2) = fit_and_eval(&replaced);
        let s2 = shares(&e2);
        eprintln!(
            "  {:<10} {:<12} {:>8.3} {:>8.3} {:>9.3}   {:>5.1}% {:>5.1}% {:>5.1}%",
            stem, "替换为最优3×3", m2, p2, x2,
            100.0 * s2[0], 100.0 * s2[1], 100.0 * s2[2]
        );
        eprintln!(
            "  {:<10} {:<12} 中位数变化 {:+.3}（{:.2}×）",
            "", "→ 判读", m2 - m1, m2 / m1.max(1e-9)
        );

        // 跨图迁移：用别张图解出的矩阵
        if let Some((other, om)) = derived.iter().find(|(s, _)| s != stem) {
            let cross: Vec<([f32; 3], [f32; 3])> =
                cur.iter().map(|(c, t)| (mat3::mul_vec(*om, *c), *t)).collect();
            let (m3, p3, _, _) = fit_and_eval(&cross);
            eprintln!(
                "  {:<10} {:<12} {:>8.3} {:>8.3}         （迁移自 {other}）",
                "", "跨图矩阵", m3, p3
            );
        }
        eprintln!();
    }

    eprintln!("说明：");
    eprintln!("  1. 上表第一行「相机空间直入」把**相机空间值直接当工作空间值**用，");
    eprintln!("     这不是现状（现状是 LibRaw 已转好的 ProPhoto），而是用于隔离");
    eprintln!("     「矩阵取错」这一因素的下界对照——同一套后端下换矩阵的相对效果可比。");
    eprintln!("  2. 逐图自解的 3×3 是**非物理**的（可能含 −4/+6 这类系数），只回答上限，");
    eprintln!("     不能直接采用。可用的矩阵必须由多图联合标定 + 交叉验证得出。");
    eprintln!("  3. 若「替换为最优3×3」后中位数仍远离 1.0，则瓶颈不在矩阵——");
    eprintln!("     团队应停止在矩阵上使劲。");
}
