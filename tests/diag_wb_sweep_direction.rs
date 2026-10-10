//! 色温扫描的方向核对：**直接从导出的参考图量"渲染得有多暖"**。
//!
//! # 为什么要有这个测试
//!
//! 用户拍了一组色温扫描（DSC_0561~0570），说「562 开始是逐渐变暖的」。而按
//! `cam_mul` 的 `B/G` 读出来的顺序是**递减**的，我最初据此把 `B/G` 大的一端标成
//! "最暖"——**方向可能标反了**。
//!
//! 争论这个没有意义，因为 `cam_mul` 反映的是**用户设的色温值**（一个设定），
//! 不是环境光。而"渲染得暖不暖"是**可以直接从导出图量出来的**：
//!
//! ```text
//! 画面越暖 → 红相对蓝越高 → 平均 R/B 越大
//! ```
//!
//! 本测试就量这个。它不依赖对 `cam_mul` 语义的任何解释，所以能独立判定方向。

mod common;

use common::samples_dir;

/// 一张导出图的平均通道值与暖度指标。
struct Warmth {
    stem: String,
    bg: f64,
    r: f64,
    g: f64,
    b: f64,
    /// 平均 R/B——越大越暖。
    warmth: f64,
}

#[test]
fn reference_exports_get_warmer_in_the_expected_order() {
    let dir = samples_dir();
    let mut rows: Vec<Warmth> = Vec::new();

    for i in 561..=570 {
        let stem = format!("DSC_{i:04}");
        let nef = dir.join(format!("{stem}.NEF"));
        let tif = dir.join(format!("{stem}.TIF"));
        if !nef.is_file() || !tif.is_file() {
            eprintln!("跳过 {stem}：缺文件");
            continue;
        }
        let Ok(wb) = nikonrawview::libraw::read_wb(&nef) else { continue };
        let bg = (wb[2] / wb[1]) as f64;

        let Ok(data) = std::fs::read(&tif) else { continue };
        let Ok(img) = nikonrawview::fit::read_rgb16(&data) else { continue };
        let Ok(plan) = nikonrawview::icc::plan_for(&data) else { continue };
        let px = nikonrawview::fit::reference_to_working(&img, &plan);

        // 只统计既不过黑也不过亮的像素：两端都有裁切/噪声的干扰
        let mut sum = [0f64; 3];
        let mut n = 0u64;
        for p in px.iter() {
            let mx = p.iter().cloned().fold(f32::MIN, f32::max);
            if mx >= 0.995 || mx <= 0.002 {
                continue;
            }
            // 排除近中性的像素太少的画面没有意义，这里不筛色度
            for k in 0..3 {
                sum[k] += p[k] as f64;
            }
            n += 1;
        }
        if n < 1000 {
            eprintln!("跳过 {stem}：可用像素仅 {n}");
            continue;
        }
        let mean = [sum[0] / n as f64, sum[1] / n as f64, sum[2] / n as f64];
        rows.push(Warmth {
            stem,
            bg,
            r: mean[0],
            g: mean[1],
            b: mean[2],
            warmth: mean[0] / mean[2].max(1e-12),
        });
    }

    if rows.is_empty() {
        eprintln!("跳过：没有可用样本");
        return;
    }

    eprintln!("=== 色温扫描：导出图的实测暖度 ===");
    eprintln!("  按文件名顺序（你拍摄的顺序）");
    eprintln!(
        "  {:<10} {:>9} {:>9} {:>9} {:>9} {:>12}",
        "图", "B/G", "平均R", "平均G", "平均B", "暖度 R/B"
    );
    for w in &rows {
        eprintln!(
            "  {:<10} {:>9.4} {:>9.5} {:>9.5} {:>9.5} {:>12.4}",
            w.stem, w.bg, w.r, w.g, w.b, w.warmth
        );
    }

    // 与 B/G 的相关：若"B/G 小 = 更暖"成立，两者应负相关
    let n = rows.len() as f64;
    let mx = rows.iter().map(|w| w.bg).sum::<f64>() / n;
    let my = rows.iter().map(|w| w.warmth).sum::<f64>() / n;
    let sxy: f64 = rows.iter().map(|w| (w.bg - mx) * (w.warmth - my)).sum();
    let sxx: f64 = rows.iter().map(|w| (w.bg - mx).powi(2)).sum();
    let syy: f64 = rows.iter().map(|w| (w.warmth - my).powi(2)).sum();
    let r = if sxx > 0.0 && syy > 0.0 { sxy / (sxx * syy).sqrt() } else { f64::NAN };

    eprintln!("\n  B/G 与实测暖度 R/B 的相关系数 = {r:+.4}");

    // 顺序单调性
    let inc = rows.windows(2).filter(|w| w[1].warmth > w[0].warmth).count();
    let dec = rows.windows(2).filter(|w| w[1].warmth < w[0].warmth).count();
    eprintln!("  按文件名顺序：暖度上升 {inc} 次，下降 {dec} 次");

    eprintln!("\n判读：");
    if r < -0.5 {
        eprintln!("  **B/G 越小 → 渲染越暖**（负相关）。");
        eprintln!("  即 cam_mul 的 B/G 反映的是「设的色温值」：设 K 高 → 渲染暖 → B/G 小。");
        eprintln!("  我之前把 B/G 大标成\"最暖\"是**反的**。");
    } else if r > 0.5 {
        eprintln!("  B/G 越大 → 渲染越暖（正相关），我原来的标注方向是对的。");
    } else {
        eprintln!("  相关性弱（{r:+.4}）——两件事可能不是简单的一一对应，需要另找办法判定。");
    }
}
