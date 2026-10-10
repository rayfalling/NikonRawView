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
//!
//! # 为什么还要把 ΔE00 拆开
//!
//! 报告只给出「中位数多少、达没达标」，却指不出下一步该改什么：误差几乎全在明度上，
//! 说明色调曲线的定义域还不够细；几乎全在彩度与色相上，说明该动三维 LUT 或色相混合器。
//! 因此 [`decompose`] 把 ΔE00 按来源拆成明度、彩度、色相与旋转交叉项四项，并给出各自
//! 占总平方的比例，[`summarize`] 再把一组样本汇总成「合计占比 / 中位数占比 / 主因票数」
//! 三个口径，供诊断报告直接引用。

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

/// CIEDE2000 的中间量：三项归一化差值，加一个旋转交叉项系数。
///
/// # 为什么抽成一个内部函数
///
/// [`ciede2000`] 与 [`decompose`] 必须共用同一份中间量。拆解的全部意义是解释总差值，若
/// 两条路径各算一遍，任何一处公式漂移都会让「各项之和」与总数悄悄对不上——而只看总数，
/// 这种不一致是发现不了的。
struct Terms {
    /// 明度项 `ΔL'/S_L`（带符号）。
    dl: f64,
    /// 彩度项 `ΔC'/S_C`（带符号）。
    dc: f64,
    /// 色相项 `ΔH'/S_H`（带符号）。
    dh: f64,
    /// 旋转交叉项系数 `R_T`。它本身不是差值，只以 `R_T·(ΔC'/S_C)·(ΔH'/S_H)` 参与总平方和。
    rt: f64,
}

impl Terms {
    fn new(lab1: [f32; 3], lab2: [f32; 3]) -> Self {
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

        Self {
            dl: dlp / sl,
            dc: dcp / sc,
            dh: dbig_hp / sh,
            rt,
        }
    }

    /// `ΔE00 = √(ΔL'² + ΔC'² + ΔH'² + R_T·ΔC'·ΔH')`。
    ///
    /// 括号内理论上非负，但 `R_T` 的幅度最大可到 2，极端输入下交叉项能把总和拉到零以下；
    /// 开方前取 `max(0)`，与公式的通行实现一致。
    fn delta_e(&self) -> f64 {
        let sum = self.dl * self.dl
            + self.dc * self.dc
            + self.dh * self.dh
            + self.rt * self.dc * self.dh;
        sum.max(0.0).sqrt()
    }
}

/// CIEDE2000 色差（kL = kC = kH = 1）。
pub fn ciede2000(lab1: [f32; 3], lab2: [f32; 3]) -> f64 {
    Terms::new(lab1, lab2).delta_e()
}

/// ΔE00 的拆解：这个总差值由哪几部分构成。
///
/// ```text
/// ΔE00² = (ΔL'/S_L)² + (ΔC'/S_C)² + (ΔH'/S_H)² + R_T·(ΔC'/S_C)·(ΔH'/S_H)
/// ```
///
/// 最后一项是**旋转交叉项**，只在蓝区（色相角约 275°）附近且彩度足够时明显，可正可负。
/// ΔE00 因此不是一个欧氏距离：少了它，三项平方和并不等于总平方和，各项占比也加不到 1。
///
/// # 为什么需要它
///
/// 只知道「中位数 2.1、未达标」决定不了下一步改什么。把总数按来源拆开，报告才能回答
/// 「残差是明度还是色度造成的」——那两种结论对应的修法完全不同。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Decomposition {
    /// 明度项 `ΔL'/S_L`。带符号：正表示 `lab2` 更亮。
    pub dl: f64,
    /// 彩度项 `ΔC'/S_C`。带符号：正表示 `lab2` 更饱和。
    pub dc: f64,
    /// 色相项 `ΔH'/S_H`。带符号：正负取决于色相角沿短弧往哪边绕。
    pub dh: f64,
    /// 旋转交叉项系数 `R_T`。**不是差值**，实际贡献见 [`Decomposition::rotation_contribution`]。
    pub rt: f64,
    /// 总色差 ΔE00。
    pub total: f64,
}

impl Decomposition {
    /// 旋转交叉项的实际贡献 `R_T·(ΔC'/S_C)·(ΔH'/S_H)`，可为负。
    ///
    /// 与 `rt` 分开是因为 `R_T` 只是系数：要在平方和里相加的是这一项，而它只有在彩度与
    /// 色相**同时**有差时才非零。
    pub fn rotation_contribution(&self) -> f64 {
        self.rt * self.dc * self.dh
    }

    /// 明度项平方占总平方的比例。
    pub fn lightness_share(&self) -> f64 {
        self.share(self.dl * self.dl)
    }

    /// 彩度项平方占总平方的比例。
    pub fn chroma_share(&self) -> f64 {
        self.share(self.dc * self.dc)
    }

    /// 色相项平方占总平方的比例。
    pub fn hue_share(&self) -> f64 {
        self.share(self.dh * self.dh)
    }

    /// 旋转交叉项占总平方的比例，可为负——负值意味着它**抵消**了一部分色差。
    pub fn rotation_share(&self) -> f64 {
        self.share(self.rotation_contribution())
    }

    /// 色度侧合计占比 = 彩度 + 色相 + 交叉项 = `1 − 明度占比`。
    ///
    /// 整个 ΔE00 只在明度与色度两侧之间二分，交叉项属于色度侧：只有色度存在时它才非零。
    pub fn chromatic_share(&self) -> f64 {
        1.0 - self.lightness_share()
    }

    /// 四项占比，顺序 `[明度, 彩度, 色相, 旋转]`，合计为 1。
    pub fn shares(&self) -> [f64; 4] {
        [
            self.lightness_share(),
            self.chroma_share(),
            self.hue_share(),
            self.rotation_share(),
        ]
    }

    fn share(&self, part: f64) -> f64 {
        let total_sq = self.total * self.total;
        // 两个颜色完全相同时没有可归因的部分：返回 0 而不是 0/0 得到的 NaN。一个 NaN 混进
        // 成千上万个样本里会把中位数整个污染掉，而且不会报错。
        if total_sq <= 0.0 {
            0.0
        } else {
            part / total_sq
        }
    }
}

/// 把 ΔE00 拆成明度 / 彩度 / 色相 / 旋转四项。
pub fn decompose(lab1: [f32; 3], lab2: [f32; 3]) -> Decomposition {
    let t = Terms::new(lab1, lab2);
    Decomposition {
        dl: t.dl,
        dc: t.dc,
        dh: t.dh,
        rt: t.rt,
        total: t.delta_e(),
    }
}

/// 对 ProPhoto 线性的一对颜色求 ΔE00。
pub fn delta_e_prophoto(a: [f32; 3], b: [f32; 3]) -> f64 {
    ciede2000(lab_from_prophoto_linear(a), lab_from_prophoto_linear(b))
}

/// 对 ProPhoto 线性的一对颜色求 ΔE00 并拆解。
pub fn delta_e_prophoto_decomposed(a: [f32; 3], b: [f32; 3]) -> Decomposition {
    decompose(lab_from_prophoto_linear(a), lab_from_prophoto_linear(b))
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

/// 升序排序；无法比较的值（NaN）按相等处理。
///
/// 与 `report` 里对 ΔE00 排序的做法一致：`partial_cmp` 返回 `None` 时退化成保持原序，
/// 而不是 panic。统计量可能因此偏一点，但不会因为一个异常样本让整份报告崩掉。
fn sorted_ascending(values: impl Iterator<Item = f64>) -> Vec<f64> {
    let mut v: Vec<f64> = values.collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v
}

/// 百分位：nearest-rank——升序排列后取第 `ceil(q·n)` 个（1 基），即下标 `ceil(q·n) - 1`。
///
/// 曾经用 `round(q·(n-1))`：n=100、q=0.5 时得 round(49.5)=50，落在**第 51 个**
/// 元素上，越过中点一个。对双峰分布（如一半样本完全无差异、一半有差异）会把中位数
/// 直接报到有差异那一侧，把结论从"达标"翻成"不达标"。
fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let rank = (q * sorted.len() as f64).ceil().max(1.0) as usize;
    sorted[rank.saturating_sub(1).min(sorted.len() - 1)]
}

/// 由一对对「我们算出的值, 参考值」生成报告。
pub fn report(pairs: &[([f32; 3], [f32; 3])], worst_n: usize) -> Report {
    let mut all: Vec<(usize, f64, [f32; 3], [f32; 3])> = pairs
        .iter()
        .enumerate()
        .map(|(i, (a, b))| (i, delta_e_prophoto(*a, *b), *a, *b))
        .collect();

    let des = sorted_ascending(all.iter().map(|x| x.1));

    all.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let worst = all.into_iter().take(worst_n).collect();

    Report {
        count: des.len(),
        median: percentile(&des, 0.5),
        p95: percentile(&des, 0.95),
        max: des.last().copied().unwrap_or(f64::NAN),
        worst,
        median_target: 1.0,
        p95_target: 3.0,
    }
}

/// 拆解的四项，顺序与 [`Decomposition::shares`]、[`DecompositionSummary`] 各数组一致。
pub const COMPONENTS: [&str; 4] = ["明度", "彩度", "色相", "旋转交叉项"];

/// 一组配对的拆解汇总——诊断报告直接用它。
///
/// # 三个口径各回答什么，必须一起看
///
/// 1. [`Self::total_share`]（合计占比）：全部样本的平方误差里每一项占了多少。这是唯一
///    可加、合计恰为 1 的口径，回答「误差量主要由谁承担」。
/// 2. 四个「占比中位数」字段与 [`Self::median_delta_e`]：典型样本长什么样。
/// 3. [`Self::primary_fraction`]：每个成分在多大比例的样本里是最大的一项。
///
/// # 为什么中位数占比不能单独当作结论
///
/// 占比是**逐样本**算的：每个样本自己的四项加起来才是 1，而各样本的主因并不相同。误差
/// 在样本间"轮流"由不同成分主导时，四项的中位数会同时偏小——实测留出集上四项中位数之和
/// 只有 0.657，而合计占比之和恒为 1；两者差多少，正说明主因在样本间换了多少。所以判断
/// 「残差是明度还是色度造成的」要用合计占比或主因票数，中位数只用来看典型样本。
#[derive(Debug, Clone)]
pub struct DecompositionSummary {
    /// 参与统计的样本数。
    pub count: usize,
    /// ΔE00 的中位数，应与 [`report`] 给出的中位数一致。
    pub median_delta_e: f64,
    /// 全部样本的 `ΔE00²` 之和。为 0 表示这组样本毫无差异，各项占比无从谈起。
    pub total_squared_error: f64,
    /// 明度项占总平方比例的**中位数**。
    pub lightness_share: f64,
    /// 彩度项占总平方比例的中位数。
    pub chroma_share: f64,
    /// 色相项占总平方比例的中位数。
    pub hue_share: f64,
    /// 旋转交叉项占总平方比例的中位数，可为负。
    pub rotation_share: f64,
    /// 全部样本平方误差的**合计占比**，顺序 `[明度, 彩度, 色相, 旋转]`，合计恰为 1。
    ///
    /// 按每个样本的 `ΔE00²` 加权，所以误差大的样本说话更响——这正是「误差量由谁承担」
    /// 该有的口径，与中位数口径回答的不是同一个问题。
    pub total_share: [f64; 4],
    /// 各成分成为单样本主因（四项占比最大）的样本比例，顺序同上，合计为 1。
    ///
    /// 旋转项只有为正时才可能胜出：它一旦为负就是在**抵消**色差，不该算作主因。
    pub primary_fraction: [f64; 4],
}

impl DecompositionSummary {
    /// 四项占比的中位数，顺序 `[明度, 彩度, 色相, 旋转]`。
    pub fn median_shares(&self) -> [f64; 4] {
        [
            self.lightness_share,
            self.chroma_share,
            self.hue_share,
            self.rotation_share,
        ]
    }

    /// 色度侧合计的**中位数**占比（含交叉项），由 `1 − 明度占比` 得到。
    pub fn chromatic_share(&self) -> f64 {
        1.0 - self.lightness_share
    }

    /// 色度侧合计的**合计**占比（含交叉项），由 `1 − 明度合计占比` 得到。
    pub fn total_chromatic_share(&self) -> f64 {
        1.0 - self.total_share[0]
    }

    /// 误差量更大的那一侧：`"明度"` 或 `"色度"`；无样本或毫无差异时是 `"未知"`。
    ///
    /// 判断用合计占比而不是中位数占比：后者会因主因在样本间轮换而同时偏小。
    pub fn dominant(&self) -> &'static str {
        if self.count == 0 || self.total_squared_error <= 0.0 {
            return "未知";
        }
        if self.total_share[0] >= self.total_chromatic_share() {
            "明度"
        } else {
            "色度"
        }
    }

    /// 一句话结论，可直接写进诊断报告。
    pub fn verdict(&self) -> String {
        if self.count == 0 {
            return String::from("样本为空，无法判断残差来源");
        }
        if self.total_squared_error <= 0.0 {
            return format!("{} 个样本全部无差异，没有残差可拆", self.count);
        }
        format!(
            "残差主因：{}。全部误差量的构成：{}；占比中位数（典型样本，四项不可加）：{}；\
             单样本主因分布：{}；ΔE00 中位数 {:.3}",
            self.dominant(),
            shares_text(&self.total_share),
            shares_text(&self.median_shares()),
            shares_text(&self.primary_fraction),
            self.median_delta_e,
        )
    }
}

/// 把四项占比排成一行。旋转项带符号显示——负值意味着它在抵消色差。
fn shares_text(shares: &[f64; 4]) -> String {
    COMPONENTS
        .iter()
        .zip(shares)
        .enumerate()
        .map(|(i, (name, v))| {
            if i == 3 {
                format!("{name} {:+.1}%", v * 100.0)
            } else {
                format!("{name} {:.1}%", v * 100.0)
            }
        })
        .collect::<Vec<_>>()
        .join("、")
}

/// 对一组 Lab 配对汇总拆解的各项占比（中位数 / 合计 / 主因票数三个口径）。
pub fn summarize(pairs: &[([f32; 3], [f32; 3])]) -> DecompositionSummary {
    summarize_decompositions(pairs.iter().map(|(a, b)| decompose(*a, *b)))
}

/// 对一组 ProPhoto 线性配对汇总拆解的各项占比（口径同 [`summarize`]）。
pub fn summarize_prophoto(pairs: &[([f32; 3], [f32; 3])]) -> DecompositionSummary {
    summarize_decompositions(pairs.iter().map(|(a, b)| delta_e_prophoto_decomposed(*a, *b)))
}

/// 逐样本拆解，再按三个口径分别汇总。
fn summarize_decompositions(items: impl Iterator<Item = Decomposition>) -> DecompositionSummary {
    let all: Vec<Decomposition> = items.collect();
    let median =
        |f: fn(&Decomposition) -> f64| percentile(&sorted_ascending(all.iter().map(f)), 0.5);

    // 合计占比：按每个样本的 total² 加权，因此各分量可加、合计恰为 1
    let total_squared_error: f64 = all.iter().map(|d| d.total * d.total).sum();
    let mut total_share = [0.0f64; 4];
    let mut primary_fraction = [0.0f64; 4];
    for d in &all {
        let shares = d.shares();
        for (k, part) in shares.iter().enumerate() {
            total_share[k] += part * d.total * d.total;
        }
        let mut best = 0usize;
        for k in 1..shares.len() {
            if shares[k] > shares[best] {
                best = k;
            }
        }
        primary_fraction[best] += 1.0;
    }
    // 全无差异时保持全 0：除零会得到 NaN，而 NaN 会一路污染到报告里
    if total_squared_error > 0.0 {
        for v in &mut total_share {
            *v /= total_squared_error;
        }
    }
    if !all.is_empty() {
        for v in &mut primary_fraction {
            *v /= all.len() as f64;
        }
    }

    DecompositionSummary {
        count: all.len(),
        median_delta_e: median(|d| d.total),
        total_squared_error,
        lightness_share: median(Decomposition::lightness_share),
        chroma_share: median(Decomposition::chroma_share),
        hue_share: median(Decomposition::hue_share),
        rotation_share: median(Decomposition::rotation_share),
        total_share,
        primary_fraction,
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
    ///
    /// 单独成一个函数是为了让总量校验与拆解校验用**同一份**数据：抄成两份，
    /// 迟早有一份会被改歪。
    fn sharma_vectors() -> [([f32; 3], [f32; 3], f64); 16] {
        [
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
        ]
    }

    #[test]
    fn matches_published_test_vectors() {
        let cases = sharma_vectors();
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

    #[test]
    fn decomposition_total_matches_ciede2000_on_published_vectors() {
        for (a, b, want) in sharma_vectors() {
            let d = decompose(a, b);
            assert!(
                (d.total - ciede2000(a, b)).abs() < 1e-9,
                "拆解的总值必须与 ciede2000 逐位一致：{a:?} vs {b:?}"
            );
            assert!(
                (d.total - want).abs() < 1e-4,
                "ΔE00({a:?}, {b:?}) = {:.4}，期望 {want:.4}",
                d.total
            );
            // 四项占比必须凑成 1：少了旋转交叉项就会差一点，而那小一截正是最容易漏掉的
            let sum: f64 = d.shares().iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "占比合计应为 1，得到 {sum}");
        }
    }

    #[test]
    fn decomposition_fields_rebuild_the_total() {
        // 公开字段若与总值的算法脱节（例如把 R_T 的贡献漏在结构体外），这里就崩
        for (a, b, _) in sharma_vectors() {
            let d = decompose(a, b);
            let rebuilt = (d.dl * d.dl
                + d.dc * d.dc
                + d.dh * d.dh
                + d.rotation_contribution())
            .max(0.0)
            .sqrt();
            assert!(
                (rebuilt - d.total).abs() < 1e-12,
                "{a:?} vs {b:?}：由字段重建的 {rebuilt} 与原值 {} 不符",
                d.total
            );
        }
    }

    #[test]
    fn achromatic_pair_has_no_chroma_or_hue_term() {
        // a = b = 0：C' 恒为 0，彩度与色相两项都该是零，色差全部由明度解释
        let d = decompose([50.0, 0.0, 0.0], [70.0, 0.0, 0.0]);
        assert!(d.dc.abs() < 1e-15, "无彩度对不应有彩度项：{}", d.dc);
        assert!(d.dh.abs() < 1e-15, "无彩度对不应有色相项：{}", d.dh);
        assert!(d.rotation_contribution().abs() < 1e-15, "交叉项应为零");
        assert!(d.dl.abs() > 1.0, "明度项应显著非零：{}", d.dl);
        assert!((d.lightness_share() - 1.0).abs() < 1e-12);
        assert!(d.chromatic_share().abs() < 1e-12);
        assert!((d.total - d.dl.abs()).abs() < 1e-12);
    }

    #[test]
    fn pure_lightness_difference_is_all_lightness() {
        // 同一 a、b 只改 L：C' 与色相角都不变，色度侧必须整块为零
        let d = decompose([40.0, 20.0, -30.0], [70.0, 20.0, -30.0]);
        assert!(d.dc.abs() < 1e-15, "纯明度差不应有彩度项：{}", d.dc);
        assert!(d.dh.abs() < 1e-15, "纯明度差不应有色相项：{}", d.dh);
        assert!(d.dl.abs() > 1.0, "明度项应显著非零：{}", d.dl);
        assert!((d.lightness_share() - 1.0).abs() < 1e-12);
        assert!(d.chromatic_share().abs() < 1e-12);
        assert!((d.total - d.dl.abs()).abs() < 1e-12);
    }

    #[test]
    fn pure_chroma_difference_is_all_chroma() {
        // 同 L、同色相角（都在 a 轴上）只改彩度：明度侧必须为零
        let d = decompose([50.0, 10.0, 0.0], [50.0, 30.0, 0.0]);
        assert!(d.dl.abs() < 1e-15, "纯彩度差不应有明度项：{}", d.dl);
        assert!(d.dh.abs() < 1e-15, "同色相角不应有色相项：{}", d.dh);
        assert!(d.lightness_share().abs() < 1e-12);
        assert!((d.chroma_share() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rotation_term_is_reported_in_the_blue_region() {
        // 蓝区（色相角接近 275°）且彩度足够时 R_T 才明显。没有这一条，旋转占比即使恒为
        // 零也不会被发现——而那等于拆解里少了一块。
        let d = decompose([50.0, 0.0, -30.0], [55.0, 10.0, -25.0]);
        eprintln!("蓝区样本的拆解 = {d:?}，占比 = {:?}", d.shares());
        assert!(d.rt.abs() > 0.5, "R_T 在蓝区应明显非零：{}", d.rt);
        assert!(
            d.rotation_share().abs() > 0.01,
            "旋转项应占到可观比例：{}",
            d.rotation_share()
        );
        let sum: f64 = d.shares().iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn prophoto_decomposition_matches_prophoto_delta_e() {
        let (a, b) = ([0.2f32, 0.3, 0.4], [0.35f32, 0.25, 0.5]);
        let d = delta_e_prophoto_decomposed(a, b);
        assert!((d.total - delta_e_prophoto(a, b)).abs() < 1e-9);

        // 完全相同的颜色：占比应是 0 而不是 0/0 得到的 NaN
        let same = delta_e_prophoto_decomposed(a, a);
        assert!(same.total < 1e-9);
        assert_eq!(same.shares(), [0.0; 4], "无差异时各项占比都应为 0");
    }

    #[test]
    fn summarize_separates_lightness_from_chroma() {
        // 明度主导：只有 L 在动
        let lightness: Vec<([f32; 3], [f32; 3])> = (0..20)
            .map(|i| {
                let l = 20.0 + i as f32 * 3.0;
                ([l, 12.0, -25.0], [l + 8.0, 12.0, -25.0])
            })
            .collect();
        let s = summarize(&lightness);
        eprintln!("{}", s.verdict());
        assert_eq!(s.count, 20);
        assert!((s.lightness_share - 1.0).abs() < 1e-12, "{}", s.lightness_share);
        assert_eq!(s.dominant(), "明度");
        // 汇总的中位数与 report 的中位数必须是同一个数——但要注意两者的输入空间不同：
        // `summarize` 收的是 Lab，`report` 收的是 ProPhoto 线性，所以要比就都用 ProPhoto。
        let p = summarize_prophoto(&lightness);
        let r = report(&lightness, 3);
        assert_eq!(p.count, r.count);
        assert!(
            (p.median_delta_e - r.median).abs() < 1e-12,
            "{} vs {}",
            p.median_delta_e,
            r.median
        );
        eprintln!("同一组 ProPhoto 配对：{}", p.verdict());

        // 色度主导：只有彩度在动
        let chroma: Vec<([f32; 3], [f32; 3])> = (0..20)
            .map(|i| {
                let c = 5.0 + i as f32 * 2.0;
                ([50.0, c, 0.0], [50.0, c + 10.0, 0.0])
            })
            .collect();
        let s = summarize(&chroma);
        eprintln!("{}", s.verdict());
        assert!(s.lightness_share.abs() < 1e-12, "{}", s.lightness_share);
        assert_eq!(s.dominant(), "色度");
        assert!((s.chroma_share - 1.0).abs() < 1e-12);
        assert!(s.verdict().contains("色度"));
    }

    #[test]
    fn summarize_handles_empty_input() {
        let s = summarize(&[]);
        assert_eq!(s.count, 0);
        assert_eq!(s.dominant(), "未知");
        assert!(s.verdict().contains("样本为空"), "{}", s.verdict());
    }

    #[test]
    fn summary_without_any_difference_reports_no_residual() {
        let pairs: Vec<([f32; 3], [f32; 3])> =
            (0..10).map(|i| { let v = [i as f32 / 20.0, 0.3, 0.4]; (v, v) }).collect();
        let s = summarize_prophoto(&pairs);
        assert_eq!(s.count, 10);
        assert!(s.total_squared_error <= 0.0);
        assert_eq!(s.total_share, [0.0; 4], "毫无差异时占比应全为 0 而不是 NaN");
        assert_eq!(s.dominant(), "未知");
        assert!(s.verdict().contains("全部无差异"), "{}", s.verdict());
    }

    #[test]
    fn summary_distinguishes_aggregate_from_median() {
        // 一半样本只有明度差，一半只有彩度差。主因在样本间"轮流"出现时，四项的中位数会
        // 同时趋近 0，而合计占比仍各占一半——这正是文档里警告的那个陷阱，必须被钉住。
        let mut pairs = Vec::new();
        for i in 0..20 {
            if i % 2 == 0 {
                pairs.push(([40.0f32, 10.0, -20.0], [60.0f32, 10.0, -20.0]));
            } else {
                pairs.push(([50.0f32, 10.0, 0.0], [50.0f32, 40.0, 0.0]));
            }
        }
        let s = summarize(&pairs);
        eprintln!("{}", s.verdict());

        let total_sum: f64 = s.total_share.iter().sum();
        assert!((total_sum - 1.0).abs() < 1e-12, "合计占比必须可加：{total_sum}");
        assert!(
            (s.total_share[0] + s.total_chromatic_share() - 1.0).abs() < 1e-12,
            "明度 + 色度两侧必须恰好构成全部"
        );
        assert!(s.total_share[0] > 0.0 && s.total_share[1] > 0.0, "两侧都该有份额");

        let primary_sum: f64 = s.primary_fraction.iter().sum();
        assert!((primary_sum - 1.0).abs() < 1e-12);
        assert!((s.primary_fraction[0] - 0.5).abs() < 1e-12, "{:?}", s.primary_fraction);
        assert!((s.primary_fraction[1] - 0.5).abs() < 1e-12, "{:?}", s.primary_fraction);

        let median_sum: f64 = s.median_shares().iter().sum();
        assert!(
            median_sum < 0.5,
            "中位数占比不可加，主因轮换时应远小于 1：{median_sum}"
        );
    }
}
