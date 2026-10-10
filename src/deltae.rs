//! 色差评估：CIEDE2000。
//!
//! # 为什么不能看线性域的绝对值
//!
//! 拟合残差原先用「ProPhoto 线性的绝对差」衡量（实测 0.008522）。那个数字**不能当
//! 感知差异用**：线性域的 0.008 在暗部与在亮部对应的视觉差异相差很多，而人眼对暗部
//! 远比对亮部敏感。要判断"够不够好"，必须换算到感知均匀的空间。
//!
//! # 为什么是 ΔE00
//!
//! CIEDE2000 是目前 CIE 推荐的色差公式，它修正了 CIELAB 在蓝色区、低彩度区与明度
//! 方向的已知偏差——那三处恰好都是基准渲染关心的区域。**JND（刚可察觉差异）约 2.3**，
//! 本项目据此把目标定为中位数 ≤ 1.0、P95 ≤ 3.0。
//!
//! # 实现来源
//!
//! 公式按 CIE 技术报告 224:2017 实现。正确性用 Sharma 等人的公开测试向量校验——
//! 那组数据专门用来暴露 ΔE00 实现里最容易写错的几处（色相角环绕、`T` 项的度数、
//! `RT` 的符号）。

use crate::mat3;

/// D50 白点（ICC PCS 的标准光源）。
pub const WHITE_D50: [f32; 3] = [0.964_22, 1.0, 0.825_21];

/// ProPhoto 线性 → CIELAB（D50）。
pub fn lab_from_prophoto_linear(rgb: [f32; 3]) -> [f32; 3] {
    let xyz = mat3::mul_vec(prophoto_to_xyz(), rgb);
    xyz_to_lab(xyz, WHITE_D50)
}

/// ProPhoto 线性 → XYZ（D50）。矩阵按需取逆，避免手抄一份转置写错。
fn prophoto_to_xyz() -> mat3::Mat3 {
    mat3::inverse(mat3::XYZ_TO_PROPHOTO).unwrap_or(mat3::IDENTITY)
}

/// XYZ → CIELAB。
pub fn xyz_to_lab(xyz: [f32; 3], white: [f32; 3]) -> [f32; 3] {
    const EPS: f32 = 216.0 / 24389.0;
    const KAPPA: f32 = 24389.0 / 27.0;
    let f = |t: f32| {
        if t > EPS {
            t.cbrt()
        } else {
            (KAPPA * t + 16.0) / 116.0
        }
    };
    let fx = f(xyz[0] / white[0]);
    let fy = f(xyz[1] / white[1]);
    let fz = f(xyz[2] / white[2]);
    [116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

/// CIEDE2000 色差（kL = kC = kH = 1）。
pub fn ciede2000(lab1: [f32; 3], lab2: [f32; 3]) -> f64 {
    let (l1, a1, b1) = (lab1[0] as f64, lab1[1] as f64, lab1[2] as f64);
    let (l2, a2, b2) = (lab2[0] as f64, lab2[1] as f64, lab2[2] as f64);

    let c1 = (a1 * a1 + b1 * b1).sqrt();
    let c2 = (a2 * a2 + b2 * b2).sqrt();
    let c_bar = (c1 + c2) / 2.0;
    let c7 = c_bar.powi(7);
    let g = 0.5 * (1.0 - (c7 / (c7 + 25f64.powi(7))).sqrt());

    let a1p = (1.0 + g) * a1;
    let a2p = (1.0 + g) * a2;
    let c1p = (a1p * a1p + b1 * b1).sqrt();
    let c2p = (a2p * a2p + b2 * b2).sqrt();

    // 色相角（度），两分量皆零时按 0 处理
    let hp = |b: f64, ap: f64| -> f64 {
        if b.abs() < 1e-12 && ap.abs() < 1e-12 {
            0.0
        } else {
            let h = b.atan2(ap).to_degrees();
            if h < 0.0 {
                h + 360.0
            } else {
                h
            }
        }
    };
    let h1p = hp(b1, a1p);
    let h2p = hp(b2, a2p);

    let dlp = l2 - l1;
    let dcp = c2p - c1p;

    let dhp = if c1p * c2p == 0.0 {
        0.0
    } else {
        let d = h2p - h1p;
        if d.abs() <= 180.0 {
            d
        } else if d > 180.0 {
            d - 360.0
        } else {
            d + 360.0
        }
    };
    let dbig_hp = 2.0 * (c1p * c2p).sqrt() * (dhp.to_radians() / 2.0).sin();

    let l_bar_p = (l1 + l2) / 2.0;
    let c_bar_p = (c1p + c2p) / 2.0;

    let h_bar_p = if c1p * c2p == 0.0 {
        h1p + h2p
    } else {
        let d = (h1p - h2p).abs();
        if d <= 180.0 {
            (h1p + h2p) / 2.0
        } else if h1p + h2p < 360.0 {
            (h1p + h2p + 360.0) / 2.0
        } else {
            (h1p + h2p - 360.0) / 2.0
        }
    };

    let t = 1.0 - 0.17 * (h_bar_p - 30.0).to_radians().cos()
        + 0.24 * (2.0 * h_bar_p).to_radians().cos()
        + 0.32 * (3.0 * h_bar_p + 6.0).to_radians().cos()
        - 0.20 * (4.0 * h_bar_p - 63.0).to_radians().cos();

    let d_theta = 30.0 * (-(((h_bar_p - 275.0) / 25.0).powi(2))).exp();
    let rc = 2.0 * ((c_bar_p.powi(7)) / (c_bar_p.powi(7) + 25f64.powi(7))).sqrt();

    let sl = 1.0 + (0.015 * (l_bar_p - 50.0).powi(2)) / (20.0 + (l_bar_p - 50.0).powi(2)).sqrt();
    let sc = 1.0 + 0.045 * c_bar_p;
    let sh = 1.0 + 0.015 * c_bar_p * t;
    let rt = -(2.0 * d_theta.to_radians()).sin() * rc;

    let x = dlp / sl;
    let y = dcp / sc;
    let z = dbig_hp / sh;

    (x * x + y * y + z * z + rt * y * z).max(0.0).sqrt()
}

/// 对 ProPhoto 线性的一对颜色求 ΔE00。
pub fn delta_e_prophoto(a: [f32; 3], b: [f32; 3]) -> f64 {
    ciede2000(lab_from_prophoto_linear(a), lab_from_prophoto_linear(b))
}

/// 一份质量报告。
#[derive(Debug, Clone)]
pub struct Report {
    /// 参与统计的样本数。
    pub count: usize,
    /// ΔE00 中位数。
    pub median: f64,
    /// ΔE00 的 95 百分位。
    pub p95: f64,
    /// ΔE00 最大值。
    pub max: f64,
    /// 误差最大的若干样本：`(样本序号, ΔE00, 我们的线性值, 参考的线性值)`。
    pub worst: Vec<(usize, f64, [f32; 3], [f32; 3])>,
    /// 目标：中位数 ≤ 1.0 且 P95 ≤ 3.0。
    pub median_target: f64,
    pub p95_target: f64,
}

impl Report {
    /// 是否达标。
    pub fn meets_target(&self) -> bool {
        self.median <= self.median_target && self.p95 <= self.p95_target
    }

    /// 未达标时说明差在哪。
    pub fn verdict(&self) -> String {
        if self.meets_target() {
            format!(
                "达标：中位数 {:.3} ≤ {:.1}，P95 {:.3} ≤ {:.1}（JND 约 2.3）",
                self.median, self.median_target, self.p95, self.p95_target
            )
        } else {
            let mut v = String::from("**未达标**：");
            if self.median > self.median_target {
                v.push_str(&format!(
                    "中位数 {:.3} 超出目标 {:.1}；",
                    self.median, self.median_target
                ));
            }
            if self.p95 > self.p95_target {
                v.push_str(&format!("P95 {:.3} 超出目标 {:.1}；", self.p95, self.p95_target));
            }
            v
        }
    }
}

/// 由一对对「我们算出的值, 参考值」生成报告。
pub fn report(pairs: &[([f32; 3], [f32; 3])], worst_n: usize) -> Report {
    let mut all: Vec<(usize, f64, [f32; 3], [f32; 3])> = pairs
        .iter()
        .enumerate()
        .map(|(i, (a, b))| (i, delta_e_prophoto(*a, *b), *a, *b))
        .collect();

    let mut des: Vec<f64> = all.iter().map(|x| x.1).collect();
    des.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // 百分位取法：nearest-rank——升序排列后取第 `ceil(q·n)` 个（1 基），即下标
    // `ceil(q·n) - 1`。
    //
    // 曾经用 `round(q·(n-1))`：n=100、q=0.5 时得 round(49.5)=50，落在**第 51 个**
    // 元素上，越过中点一个。对双峰分布（如一半样本完全无差异、一半有差异）会把中位数
    // 直接报到有差异那一侧，把结论从"达标"翻成"不达标"。
    let pick = |q: f64| -> f64 {
        if des.is_empty() {
            return f64::NAN;
        }
        let rank = (q * des.len() as f64).ceil().max(1.0) as usize;
        des[rank.saturating_sub(1).min(des.len() - 1)]
    };

    all.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let worst = all.into_iter().take(worst_n).collect();

    Report {
        count: des.len(),
        median: pick(0.5),
        p95: pick(0.95),
        max: des.last().copied().unwrap_or(f64::NAN),
        worst,
        median_target: 1.0,
        p95_target: 3.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sharma 等人的公开测试向量。
    ///
    /// 这组数据专门用来暴露 ΔE00 实现里最容易写错的几处：色相角环绕、`T` 项的度数、
    /// `RT` 的符号、以及 `C1'·C2' = 0` 的退化分支。只用"两个相同颜色差为 0"这种
    /// 自明用例是测不出来的。
    #[test]
    fn matches_published_test_vectors() {
        let cases: [([f32; 3], [f32; 3], f64); 16] = [
            ([50.0, 2.6772, -79.7751], [50.0, 0.0, -82.7485], 2.0425),
            ([50.0, 3.1571, -77.2803], [50.0, 0.0, -82.7485], 2.8615),
            ([50.0, 2.8361, -74.0200], [50.0, 0.0, -82.7485], 3.4412),
            ([50.0, -1.3802, -84.2814], [50.0, 0.0, -82.7485], 1.0000),
            ([50.0, -1.1848, -84.8006], [50.0, 0.0, -82.7485], 1.0000),
            ([50.0, -0.9009, -85.5211], [50.0, 0.0, -82.7485], 1.0000),
            ([50.0, 0.0, 0.0], [50.0, -1.0, 2.0], 2.3669),
            ([50.0, -1.0, 2.0], [50.0, 0.0, 0.0], 2.3669),
            ([50.0, 2.49, -0.001], [50.0, -2.49, 0.0009], 7.1792),
            ([50.0, 2.49, -0.001], [50.0, -2.49, 0.001], 7.1792),
            ([50.0, 2.49, -0.001], [50.0, -2.49, 0.0011], 7.2195),
            ([50.0, -0.001, 2.49], [50.0, 0.0009, -2.49], 4.8045),
            ([50.0, 2.5, 0.0], [50.0, 0.0, -2.5], 4.3065),
            ([50.0, 2.5, 0.0], [73.0, 25.0, -18.0], 27.1492),
            ([50.0, 2.5, 0.0], [61.0, -5.0, 29.0], 22.8977),
            ([50.0, 2.5, 0.0], [56.0, -27.0, -3.0], 31.9030),
        ];
        let mut worst = 0f64;
        for (a, b, want) in cases {
            let got = ciede2000(a, b);
            worst = worst.max((got - want).abs());
            assert!(
                (got - want).abs() < 1e-4,
                "ΔE00({a:?}, {b:?}) = {got:.4}，期望 {want:.4}"
            );
        }
        eprintln!("16 个公开测试向量全部通过，最大偏差 {worst:.2e}");
    }

    #[test]
    fn identical_colors_have_zero_difference() {
        for lab in [[50.0, 0.0, 0.0], [10.0, 30.0, -40.0], [90.0, -20.0, 60.0]] {
            assert!(ciede2000(lab, lab) < 1e-9);
        }
    }

    #[test]
    fn symmetry_holds() {
        let a = [50.0f32, 2.5, 0.0];
        let b = [73.0f32, 25.0, -18.0];
        assert!((ciede2000(a, b) - ciede2000(b, a)).abs() < 1e-9);
    }

    #[test]
    fn lab_of_d50_white_is_100_0_0() {
        let lab = xyz_to_lab(WHITE_D50, WHITE_D50);
        eprintln!("D50 白点的 Lab = {lab:?}");
        assert!((lab[0] - 100.0).abs() < 1e-3, "L 应为 100，得到 {}", lab[0]);
        assert!(lab[1].abs() < 1e-3 && lab[2].abs() < 1e-3, "a/b 应为 0，得到 {lab:?}");
    }

    #[test]
    fn prophoto_white_maps_to_l_100() {
        // ProPhoto 线性 (1,1,1) 应是无彩度的白
        let lab = lab_from_prophoto_linear([1.0, 1.0, 1.0]);
        eprintln!("ProPhoto 线性 (1,1,1) 的 Lab = {lab:?}");
        assert!((lab[0] - 100.0).abs() < 0.5, "L 应接近 100，得到 {}", lab[0]);
        assert!(lab[1].abs() < 1.0 && lab[2].abs() < 1.0, "应接近无彩度：{lab:?}");
    }

    #[test]
    fn report_percentiles_and_verdict() {
        // 造一组已知误差：一半 0，一半 10
        let mut pairs = Vec::new();
        for i in 0..100 {
            let v = [0.2f32, 0.3, 0.4];
            let other = if i < 50 { v } else { [0.5f32, 0.6, 0.7] };
            pairs.push((v, other));
        }
        let r = report(&pairs, 5);
        eprintln!("样本 ΔE00 抽样：{:?}", (0..100).step_by(10).map(|i| {
            let (a, b) = pairs[i];
            (i, delta_e_prophoto(a, b))
        }).collect::<Vec<_>>());
        assert_eq!(r.count, 100);
        assert!(r.median < 1e-6, "一半样本无差异，中位数应为 0：{}", r.median);
        assert!(r.p95 > 1.0, "P95 应落在有差异的那一半：{}", r.p95);
        assert_eq!(r.worst.len(), 5);
        eprintln!("{}", r.verdict());
        assert!(!r.meets_target(), "这组数据不该达标");
        assert!(r.verdict().contains("未达标"), "未达标必须明说：{}", r.verdict());
    }

    #[test]
    fn perfect_fit_meets_target() {
        let pairs: Vec<([f32; 3], [f32; 3])> =
            (0..50).map(|i| { let v = [i as f32 / 50.0; 3]; (v, v) }).collect();
        let r = report(&pairs, 3);
        assert!(r.meets_target());
        assert!(r.verdict().contains("达标"));
        eprintln!("{}", r.verdict());
    }
}
