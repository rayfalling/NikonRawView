//! 跨图泛化诊断：**「尼康的渲染」是不是一个良定义的目标？**
//!
//! # 为什么这条诊断与另外三条不同
//!
//! 另外三条（残差随亮度的分布、ΔE00 分量拆解、暗部信噪比）都在问"**我们的拟合哪里没做好**"。
//! 这一条问的是更底层的问题：**拿 A 图拟合出来的变换，套到 B 图上还成立吗？**
//!
//! 若成立，说明基准渲染确实是机型级的固定变换，残差只是拟合精度问题，可以靠改进
//! 拟合方法来解决。
//!
//! 若不成立——即同机型、同基准、同设置的不同照片之间就无法用一个变换覆盖——那么
//! **继续打磨拟合方法是白费力气**，因为要拟合的目标本身随图变化。那时该做的是
//! 找出变化的来源（是相机的场景自适应？是 NX Studio 的逐图处理？还是我们的解码
//! 在不同画面下与尼康的偏差不同？）。
//!
//! # 判定口径
//!
//! 与同图留出评估的 ΔE00 直接对比：若跨图的 ΔE00 与同图**同量级**，则变换是有效的
//! 机型级变换；若显著更差，则目标随图变化。

mod common;

use common::{fitted_fixture_for, quantiles};

#[test]
fn does_the_transform_generalise_across_images() {
    // 用 DSC_0001 拟合
    let Some(base) = fitted_fixture_for("DSC_0001", true) else {
        eprintln!("跳过：缺少 DSC_0001 样本");
        return;
    };

    let mut same_img: Vec<f64> = base.eval.iter().map(|s| s.delta_e()).collect();
    let (m0, p0, x0) = quantiles(&mut same_img);
    eprintln!("=== 参照：同图留出评估（DSC_0001 拟合 → DSC_0001 评估）===");
    eprintln!("  拟合样本 {}，评估样本 {}", base.fit_count, base.eval.len());
    eprintln!("  ΔE00 中位数 {m0:.3}  P95 {p0:.3}  最大 {x0:.3}");

    // 拿这份变换去套别的照片
    eprintln!("\n=== 跨图：DSC_0001 拟合 → 其他照片评估 ===");
    let mut any = false;
    for stem in ["DSC_0002", "DSC_0003", "DSC_0004"] {
        let Some(other) = fitted_fixture_for(stem, false) else {
            eprintln!("  {stem}：样本缺席，跳过");
            continue;
        };
        any = true;
        let mut des: Vec<f64> = other
            .eval
            .iter()
            .map(|s| {
                nikonrawview::deltae::delta_e_prophoto(base.transform.apply(s.ours_linear), s.theirs_linear)
            })
            .collect();
        let (m, p, x) = quantiles(&mut des);
        eprintln!("  {stem}（{} 个评估样本）", other.eval.len());
        eprintln!("    ΔE00 中位数 {m:.3}  P95 {p:.3}  最大 {x:.3}");
        eprintln!("    相对同图：中位数 ×{:.2}，P95 ×{:.2}", m / m0.max(1e-9), p / p0.max(1e-9));
    }

    if !any {
        eprintln!("跳过跨图部分：simple/ 中只有 DSC_0001 一对样本");
        return;
    }

    eprintln!("\n判读：跨图的中位数若与同图同量级（比如相差不到 1.5 倍），");
    eprintln!("      说明基准渲染是机型级固定变换，残差属拟合精度问题；");
    eprintln!("      若显著更差，则目标本身随图变化，打磨拟合方法无用。");
}
