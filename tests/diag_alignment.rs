//! 空间配准检查：**我们与参考的像素配对，最优偏移是不是 `(0,0)`？**
//!
//! # 为什么先做这个
//!
//! 诊断 P（stages）挖到一条比"跨图不迁移"更强的现象：
//!
//! ```text
//! 完整曲线神谕（左半拟合、右半评估）：
//!   左半（训练侧）均值改善  +0.982
//!   右半（留出侧）均值改善  −0.210      ← 同一张图内部都不迁移
//! ```
//!
//! 它的结论是"需要的修正与『画面哪里是什么颜色』绑定"。**但还有一条更朴素、更该先排除的
//! 解释：空间配准本身错了。**
//!
//! 我们的管线是**纯逐像素**的，位置配对靠裁剪边距（Z8 是 `12/8/12/8`）。
//! **若边距差一个像素，那么每一个"配对"其实都是两个相邻但不同的像素。** 后果与上面
//! 观察到的现象完全吻合：
//!
//! - 拟合会去拟合**边缘结构**而不是颜色关系 → 逐图特有、无法迁移
//! - 左右半的**边缘分布不同** → 同一张图内也不迁移
//! - 跨图更不用说 → 之前的 −1.2%
//!
//! 而它也很容易解释另一个旧观察：LUT 网格 **95.52% 为空**（若配对错了，样本云会散开）。
//!
//! # 判据
//!
//! 在 `dx, dy ∈ [-3, +3]` 上扫，找让**残差最小**的偏移：
//!
//! - 最优是 `(0,0)` → 配准正确，排除这条线，回到诊断 P 的结论
//! - 最优不是 `(0,0)` → **配准 bug**，且能直接读出差了几个像素

mod common;

use common::samples_dir;
use nikonrawview::libraw::{self, Demosaic, Options};

/// 在给定偏移下，两张图重叠区域的平均绝对差（0..1 域，三通道）。
fn residual_at(ours: &[[f32; 3]], w: usize, h: usize, theirs: &[[f32; 3]], dx: i32, dy: i32) -> Option<f64> {
    let mut sum = 0f64;
    let mut n = 0u64;
    // 逐 4 像素抽样，够快且不失代表性
    let y0 = dy.max(0) as usize;
    let y1 = (h as i32 + dy.min(0)) as usize;
    let x0 = dx.max(0) as usize;
    let x1 = (w as i32 + dx.min(0)) as usize;
    if y0 >= y1 || x0 >= x1 {
        return None;
    }
    for y in (y0..y1).step_by(4) {
        let ty = (y as i32 - dy) as usize;
        for x in (x0..x1).step_by(4) {
            let tx = (x as i32 - dx) as usize;
            let a = ours[y * w + x];
            let Some(b) = theirs.get(ty * w + tx) else { continue };
            for k in 0..3 {
                sum += (a[k] as f64 - b[k] as f64).abs();
            }
            n += 3;
        }
    }
    (n > 1000).then(|| sum / n as f64)
}

#[test]
fn pixel_pairing_is_aligned() {
    let dir = samples_dir();
    let stems = ["DSC_0001", "DSC_0141", "DSC_0567", "DSC_8562"];
    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };

    eprintln!("=== 空间配准：最优偏移是不是 (0,0)？ ===");
    eprintln!("（若最优不是 (0,0)，则裁剪边距差了几个像素——那会让每一个『配对』都是两个不同的像素）\n");

    let mut any = false;
    for stem in stems {
        let nef = dir.join(format!("{stem}.NEF"));
        let tif = dir.join(format!("{stem}.TIF"));
        if !nef.is_file() || !tif.is_file() {
            eprintln!("  跳过 {stem}：缺文件");
            continue;
        }
        let Ok(dec) = libraw::decode_working_space(&nef, &opts) else {
            eprintln!("  跳过 {stem}：解码失败");
            continue;
        };
        let Ok(tif_data) = std::fs::read(&tif) else { continue };
        let Ok(img) = nikonrawview::fit::read_rgb16(&tif_data) else { continue };
        let Ok(plan) = nikonrawview::icc::plan_for(&tif_data) else { continue };
        let theirs = nikonrawview::fit::reference_to_working(&img, &plan);

        let Some(model) = nikonrawview::camera::read_model(&nef) else { continue };
        let Some(entry) = nikonrawview::camera::lookup(&model) else { continue };
        let mg = if dec.rotated { entry.margins.rotated() } else { entry.margins };
        let Some((cw, ch)) = mg.effective(dec.width, dec.height) else { continue };
        if (cw, ch) != (img.width, img.height) {
            eprintln!("  跳过 {stem}：裁切后 {cw}×{ch} 与参考 {}×{} 不符", img.width, img.height);
            continue;
        }
        // 我们的裁切结果
        let ours: Vec<[f32; 3]> = (0..ch)
            .flat_map(|y| (0..cw).map(move |x| (x, y)))
            .map(|(x, y)| {
                let p = dec.at(x + mg.left, y + mg.top).unwrap_or([0, 0, 0]);
                [
                    p[0] as f32 / 65535.0,
                    p[1] as f32 / 65535.0,
                    p[2] as f32 / 65535.0,
                ]
            })
            .collect();

        any = true;
        let mut best = (0i32, 0i32, f64::MAX);
        let mut table = Vec::new();
        for dy in -3i32..=3 {
            for dx in -3i32..=3 {
                if let Some(r) = residual_at(&ours, cw, ch, &theirs, dx, dy) {
                    table.push((dx, dy, r));
                    if r < best.2 {
                        best = (dx, dy, r);
                    }
                }
            }
        }
        let zero = table.iter().find(|(x, y, _)| *x == 0 && *y == 0).map(|t| t.2).unwrap_or(f64::NAN);
        eprintln!(
            "  {stem:<11} 最优偏移 ({:>2},{:>2}) 残差 {:.6}    (0,0) 残差 {:.6}   改善 {:.2}%",
            best.0, best.1, best.2, zero,
            100.0 * (zero - best.2) / zero.max(1e-12)
        );
        // 打印最优那一行/列，看形状是否明显
        if best.0 != 0 || best.1 != 0 {
            eprintln!("       注：最优非 (0,0)，请检查裁剪边距");
        }
    }

    if !any {
        eprintln!("跳过：没有可用样本");
        return;
    }

    eprintln!("\n判读：");
    eprintln!("  若各图最优偏移均为 (0,0) → **配准正确**，诊断 P 的结论成立（修正与空间/色彩区域绑定，");
    eprintln!("  需要另找原因）。");
    eprintln!("  若有图最优非 (0,0) → **配准 bug**，且偏差的像素数就直接读出来了；");
    eprintln!("  那此前所有基于逐像素配对的拟合结论都要在修正配准后重测。");
}
