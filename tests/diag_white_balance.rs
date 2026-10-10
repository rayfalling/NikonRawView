//! 白平衡一致性诊断：参考导出的白平衡是否与我们的解码一致？
//!
//! # 为什么这条最该先查
//!
//! 诊断 D 留下一个异常：`DSC_0141` 的原始解码 `B/G` 比参考低 **30%**，而另两张只有
//! −7% / +6%，且该偏差与亮度近乎无关。诊断 F 进一步指出这个亏损**不随彩度变化**、
//! 是**逐图性质**。
//!
//! 一个白平衡差异**恰好就是逐图、与亮度无关、且表现为通道增益**——三条特征全中。
//! 而白平衡在 NX Studio 里是可以被编辑的：边车若改了白平衡，导出就会用改后的值，
//! 而我们的解码用的是 **NEF 里记录的相机原始白平衡**。
//!
//! 若属实，`DSC_0141` 就不是一个有效的参考样本——它测的是"两份不同的白平衡之差"，
//! 而不是"同一个白平衡下两种渲染之差"。**那会把整个跨图标定带偏。**
//!
//! 本测试把三张图的 `cam_mul` 与由参考反推的白平衡摆在一起看。

mod common;

use common::{fitted_fixture_for, samples_dir};

/// 由一对 (我们, 参考) 反推每个通道需要的增益——即"参考相对我们"的白平衡比。
///
/// 取高亮度、低彩度的像素：那里三个通道都远离噪声底，且色调曲线接近线性，
/// 增益估计最干净。
fn implied_gain(samples: &[common::DiagSample]) -> Option<[f64; 3]> {
    let mut sum = [0f64; 3];
    let mut n = 0u64;
    for s in samples {
        let lum = s.lum_ours();
        // 取中间调（曲线在两端非线性最强）
        if !(0.05..0.50).contains(&lum) {
            continue;
        }
        for (k, slot) in sum.iter_mut().enumerate() {
            let o = s.ours_linear[k] as f64;
            let t = s.theirs_linear[k] as f64;
            if o > 1e-4 && t > 1e-4 {
                *slot += t / o;
                if k == 0 {
                    n += 1;
                }
            }
        }
    }
    if n < 100 {
        return None;
    }
    let g = [sum[0] / n as f64, sum[1] / n as f64, sum[2] / n as f64];
    // 归一化到绿通道
    Some([g[0] / g[1], 1.0, g[2] / g[1]])
}

#[test]
fn white_balance_is_consistent_between_decode_and_export() {
    let dir = samples_dir();
    let mut rows: Vec<(String, [f32; 3], Option<[f64; 3]>)> = Vec::new();

    for stem in ["DSC_0001", "DSC_0141", "DSC_8562", "DSC_0002"] {
        let nef = dir.join(format!("{stem}.NEF"));
        if !nef.is_file() {
            continue;
        }
        let wb = nikonrawview::libraw::read_wb(&nef).ok();
        let Some(fx) = fitted_fixture_for(stem, false) else { continue };
        let gain = implied_gain(&fx.eval);
        rows.push((stem.to_string(), wb.unwrap_or([f32::NAN; 3]), gain));
    }

    if rows.is_empty() {
        eprintln!("跳过：没有可用样本");
        return;
    }

    eprintln!("=== 白平衡一致性 ===");
    eprintln!("  相机记录的 cam_mul（已归一化到 G=1）与「参考/我们」反推的增益");
    eprintln!();
    eprintln!(
        "  {:<10} {:>10} {:>10} {:>10}   {:>10} {:>10} {:>10}   判定",
        "图", "R", "G", "B", "反推R/G", "反推G", "反推B/G"
    );
    for (stem, wb, gain) in &rows {
        let b = wb[2] / wb[1];
        let rr = wb[0] / wb[1];
        let (gs, verdict) = match gain {
            None => ("—".to_string(), "样本不足".to_string()),
            Some(g) => {
                // 反推增益应接近 1（同一白平衡）。偏离 10% 以上值得警惕。
                let dev = (g[0] - 1.0).abs().max((g[2] - 1.0).abs());
                let v = if dev < 0.05 {
                    "一致".to_string()
                } else if dev < 0.15 {
                    format!("轻度偏离 {:.1}%", dev * 100.0)
                } else {
                    format!("**显著偏离 {:.1}%**", dev * 100.0)
                };
                (format!("{:>10.4} {:>10.4} {:>10.4}", g[0], g[1], g[2]), v)
            }
        };
        eprintln!("  {:<10} {:>10.4} {:>10.4} {:>10.4}   {}   {}",
            stem, rr, 1.0, b, gs, verdict);
    }
    eprintln!("\n判读：");
    eprintln!("  反推增益是「参考 ÷ 我们」的逐通道比值。若三分量都在 1.0 附近，");
    eprintln!("  说明两边用的是同一套白平衡；若某张图显著偏离，则该图的参考导出");
    eprintln!("  与我们的解码白平衡不一致——**它不是有效参考样本**。");
}
