//! 色彩管线：相机线性 RGB → ProPhoto 线性工作空间。
//!
//! # 矩阵从哪来
//!
//! 相机 → 输出空间的 3×3 矩阵由解码层给出（`libraw_get_rgb_cam`，输出空间设为
//! ProPhoto）。**这个矩阵已经包含了按拍摄白平衡做的光源适配**——实测同一台相机
//! 不同白平衡的两张照片得到的矩阵不同，因此本模块不再自己做光源插值。
//!
//! 相机色彩矩阵**并不在 NEF 里**（实测 DNG 标签为空），它来自解码器内置的机型
//! 标定。所以这里 MUST NOT 用公开的通用矩阵替代——通用矩阵缺少机型标定，会引入
//! 难以定位的偏色。
//!
//! # 为什么由本模块施加矩阵，而不是让解码层直接输出 ProPhoto
//!
//! 解码层在转换时会就地裁到 0..65535。高饱和颜色与负值分量会在那里被截断，
//! 之后再想处理就晚了。本模块把矩阵施加在 `f32` 上，**中途不裁切、不截断负值**。

use crate::error::{Error, Result};
use crate::libraw::Decoded;
use crate::mat3::{self, Mat3};

/// 一次色彩变换的方案。
#[derive(Debug, Clone, PartialEq)]
pub struct ColorPlan {
    matrix: Mat3,
    /// 把解码层的 16 位整数换算到 0..1 的比例。
    scale: f32,
    /// 解码层报告的白电平（原始域），仅用于报告与自检。
    color_maximum: i32,
}

/// 解码层输出的位深满量程。
///
/// 我们请求 `output_bps = 16`，解码层据此把**相机白电平映射到 65535**——
/// 实测印证：输出最大值确实是 65535，而 `libraw_get_color_maximum` 返回的
/// `15375` 是**原始域**的白电平，不是输出域的刻度。
///
/// 这一点曾经搞错过：原先按 `1 / color_maximum` 归一化，结果白电平被映射到
/// 4.26 而不是 1.0。是任务 3.3 的交叉验证把它抓出来的。
const OUTPUT_FULL_SCALE: f32 = 65535.0;

impl ColorPlan {
    /// 由解码结果构造。
    pub fn from_decoded(d: &Decoded) -> Result<Self> {
        if d.color_maximum <= 0 {
            return Err(Error::InvalidInput(format!(
                "色彩上限非正（{}），无法确认解码层的电平刻度",
                d.color_maximum
            )));
        }
        if !d.rgb_cam.iter().flatten().all(|v| v.is_finite()) {
            return Err(Error::InvalidInput("相机矩阵含非有限分量".into()));
        }
        if mat3::inverse(d.rgb_cam).is_none() {
            return Err(Error::InvalidInput("相机矩阵不可逆，说明取到的不是有效标定".into()));
        }
        Ok(Self {
            matrix: d.rgb_cam,
            scale: 1.0 / OUTPUT_FULL_SCALE,
            color_maximum: d.color_maximum,
        })
    }

    /// 由**推导**出的矩阵构造（见 [`derive_matrix`]）。
    ///
    /// 与 [`Self::from_decoded`] 的区别：那条走解码层的相机矩阵访问接口，而实测该
    /// 接口返回的矩阵与输出色彩空间无关、并非实际参与转换的那个；这条用的是自己
    /// 测出来的矩阵。
    pub fn from_derived(matrix: Mat3, full_scale: f32) -> Self {
        Self { matrix, scale: 1.0 / full_scale, color_maximum: full_scale as i32 }
    }

    /// 矩阵本身，供结果记录与手工复核。
    pub fn matrix(&self) -> Mat3 {
        self.matrix
    }

    pub fn scale(&self) -> f32 {
        self.scale
    }

    /// 解码层报告的原始域白电平。
    pub fn color_maximum(&self) -> i32 {
        self.color_maximum
    }

    /// 把 0..1 归一化值换算回 16 位输出刻度。
    pub fn to_u16(&self, v: f32) -> u16 {
        (v * OUTPUT_FULL_SCALE).clamp(0.0, OUTPUT_FULL_SCALE) as u16
    }

    /// 完整 9 分量的可读表示。
    pub fn describe(&self) -> String {
        let m = self.matrix;
        format!(
            "[{:>9.5} {:>9.5} {:>9.5}] [{:>9.5} {:>9.5} {:>9.5}] [{:>9.5} {:>9.5} {:>9.5}]",
            m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]
        )
    }

    /// 把一个 16 位相机像素换算到 0..1 并施加矩阵。
    ///
    /// **不裁切、不截断**：结果可以为负、可以大于 1，交由最终编码阶段处理。
    pub fn apply_u16(&self, px: [u16; 3]) -> [f32; 3] {
        let s = self.scale;
        self.apply([px[0] as f32 * s, px[1] as f32 * s, px[2] as f32 * s])
    }

    /// 对 0..1 的相机线性值施加矩阵。
    pub fn apply(&self, cam: [f32; 3]) -> [f32; 3] {
        mat3::mul_vec(self.matrix, cam)
    }

    /// 对整幅图像施加。输出为行主序的三通道 `f32`。
    pub fn apply_all(&self, pixels: &[u16]) -> Vec<[f32; 3]> {
        debug_assert_eq!(pixels.len() % 3, 0);
        pixels.chunks_exact(3).map(|c| self.apply_u16([c[0], c[1], c[2]])).collect()
    }
}

// ---------------------------------------------------------------------------
// 相机矩阵的推导
// ---------------------------------------------------------------------------

/// 最小二乘推导出的相机矩阵及其拟合质量。
#[derive(Debug, Clone, PartialEq)]
pub struct DerivedMatrix {
    /// 相机 → 工作空间的 3×3 矩阵。
    pub matrix: Mat3,
    /// 参与拟合的样本数。
    pub samples: usize,
    /// 残差的均方根（单位与样本同为 0..1）。
    pub rms: f64,
    /// 最大绝对残差。
    pub max_abs: f64,
}

/// 由「相机空间像素 ↔ 工作空间像素」配对最小二乘求出相机矩阵。
///
/// # 为什么这样求是可行的
///
/// 实测相机矩阵**按机型固定、不随拍摄白平衡变化**（`DSC_4569` 与 `DSC_5157` 的
/// 白平衡系数差异和 0.6152，而矩阵差异和 0.000000）。因此每个机型只需测一次，
/// 之后所有同机型照片共用——这也正是"自己推导"能成立的前提。
///
/// # 为什么不能直接问解码层要
///
/// 解码层确实有一个「相机 → 输出空间」的访问接口，但它返回的矩阵与输出色彩空间
/// **无关**（Camera / sRGB / ProPhoto 三种输出下逐位相同），而像素却随输出色彩空间
/// 变化——说明那不是实际参与转换的矩阵。与其复刻其内部数据，不如**测**出来。
///
/// # 求解方式
///
/// 对每个输出通道解一个三元最小二乘：`min Σ (p_k − m_k·c)²`，正规方程
/// `(Σ c cᵀ) m_k = Σ c p_k`。三个通道共用同一个 3×3 正规矩阵，只需求一次逆。
///
/// 样本应**只保留双方都未饱和的点**——饱和处解码层做过裁切，把它们算进拟合会把
/// 裁切行为误当成矩阵的一部分。
pub fn derive_matrix(pairs: &[([f32; 3], [f32; 3])]) -> Option<DerivedMatrix> {
    if pairs.len() < 3 {
        return None;
    }

    // 正规矩阵 S = Σ c cᵀ（对称），以及右端项 B[k] = Σ c · p_k
    let mut s = [[0f64; 3]; 3];
    let mut b = [[0f64; 3]; 3];
    for (c, p) in pairs {
        let c = [c[0] as f64, c[1] as f64, c[2] as f64];
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += c[i] * c[j];
            }
            for k in 0..3 {
                b[k][i] += c[i] * p[k] as f64;
            }
        }
    }

    // 用 f64 版的高斯-约当求逆；这里不复用 f32 的 mat3::inverse 以免精度损失
    let inv = invert3_f64(s)?;
    let mut m = [[0f32; 3]; 3];
    for k in 0..3 {
        for i in 0..3 {
            let v: f64 = (0..3).map(|j| inv[i][j] * b[k][j]).sum();
            m[k][i] = v as f32;
        }
    }

    // 残差
    let mut sum_sq = 0f64;
    let mut max_abs = 0f64;
    let mut n = 0usize;
    for (c, p) in pairs {
        let pred = mat3::mul_vec(m, *c);
        for k in 0..3 {
            let e = (pred[k] - p[k]).abs() as f64;
            sum_sq += e * e;
            max_abs = max_abs.max(e);
            n += 1;
        }
    }
    let rms = (sum_sq / n.max(1) as f64).sqrt();

    Some(DerivedMatrix { matrix: m, samples: pairs.len(), rms, max_abs })
}

fn invert3_f64(m: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-18 {
        return None;
    }
    let d = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d,
        ],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libraw::{Decoded, Demosaic};

    #[test]
    fn derive_recovers_a_known_matrix() {
        // 造一组由已知矩阵生成的配对，验证能解回来
        let truth: Mat3 = [[1.3, -0.2, -0.1], [-0.15, 1.4, -0.25], [-0.02, -0.3, 1.32]];
        let mut pairs = Vec::new();
        for i in 0..40 {
            for j in 0..40 {
                let c = [i as f32 / 40.0, j as f32 / 40.0, ((i + j) % 40) as f32 / 40.0];
                let p = mat3::mul_vec(truth, c);
                pairs.push((c, p));
            }
        }
        let d = derive_matrix(&pairs).expect("应能求解");
        assert_eq!(d.samples, 1600);
        assert!(d.rms < 1e-6, "rms={}", d.rms);
        for (i, row) in d.matrix.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!(
                    (v - truth[i][j]).abs() < 1e-4,
                    "[{i}][{j}] {v} vs {}",
                    truth[i][j]
                );
            }
        }
    }

    #[test]
    fn derive_needs_enough_samples() {
        assert!(derive_matrix(&[]).is_none());
        assert!(derive_matrix(&[([0.1, 0.2, 0.3], [0.1, 0.2, 0.3])]).is_none());
    }

    #[test]
    fn derive_rejects_degenerate_samples() {
        // 所有样本共线 → 正规矩阵奇异
        let pairs: Vec<_> = (0..50)
            .map(|i| {
                let c = [i as f32 / 50.0, 0.0, 0.0];
                (c, c)
            })
            .collect();
        assert!(derive_matrix(&pairs).is_none(), "退化样本应返回 None 而非胡乱给一个矩阵");
    }

    fn fake(matrix: Mat3, color_maximum: i32) -> Decoded {
        Decoded {
            width: 1,
            height: 1,
            pixels: vec![0; 3],
            wb: [1.0, 1.0, 1.0],
            rgb_cam: matrix,
            color_maximum,
            demosaic: Demosaic::Dht,
            output: crate::libraw::OutputColor::Camera,
            raw_width: 1,
            raw_height: 1,
        }
    }

    #[test]
    fn identity_matrix_preserves_values() {
        let d = fake(mat3::IDENTITY, 15375);
        let p = ColorPlan::from_decoded(&d).unwrap();
        // 满量程 65535 应映射到 1.0——白电平在输出域就是满量程
        let got = p.apply_u16([65535, 32768, 16384]);
        assert!((got[0] - 1.0).abs() < 1e-6, "{got:?}");
        assert!((got[1] - 0.5).abs() < 2e-4, "{got:?}");
        assert!((got[2] - 0.25).abs() < 2e-4, "{got:?}");
    }

    #[test]
    fn full_scale_is_output_range_not_color_maximum() {
        // 回归：曾按 1/color_maximum 归一化，导致白电平被映射到 65535/15375 ≈ 4.26
        let d = fake(mat3::IDENTITY, 15375);
        let p = ColorPlan::from_decoded(&d).unwrap();
        assert!((p.scale() - 1.0 / 65535.0).abs() < 1e-12, "scale={}", p.scale());
        assert_eq!(p.color_maximum(), 15375, "原始域白电平应被保留供报告");
        assert_eq!(p.to_u16(1.0), 65535);
        assert_eq!(p.to_u16(0.5), 32767);
        assert_eq!(p.to_u16(2.0), 65535, "超量程应夹到满量程");
        assert_eq!(p.to_u16(-1.0), 0);
    }

    #[test]
    fn negatives_are_not_clipped() {
        // 相机矩阵常含负系数——高饱和色会产生负分量，必须原样保留
        let m: Mat3 = [[1.0, -0.8, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let d = fake(m, 1000);
        let p = ColorPlan::from_decoded(&d).unwrap();
        let got = p.apply_u16([0, 65535, 0]);
        assert!(got[0] < 0.0, "应产生负分量而不是被截断为 0，得到 {got:?}");
        assert!((got[0] + 0.8).abs() < 1e-6);
    }

    #[test]
    fn values_above_one_are_not_clipped() {
        let m: Mat3 = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let d = fake(m, 1000);
        let p = ColorPlan::from_decoded(&d).unwrap();
        let got = p.apply_u16([65535, 0, 0]);
        assert!((got[0] - 3.0).abs() < 1e-5, "超色域值不应被裁到 1.0，得到 {got:?}");
    }

    #[test]
    fn zero_color_maximum_is_rejected() {
        let d = fake(mat3::IDENTITY, 0);
        let e = ColorPlan::from_decoded(&d).unwrap_err();
        assert!(matches!(e, Error::InvalidInput(_)), "得到 {e:?}");
    }

    #[test]
    fn singular_matrix_is_rejected() {
        let m: Mat3 = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 1.0, 1.0]];
        let d = fake(m, 1000);
        let e = ColorPlan::from_decoded(&d).unwrap_err();
        assert!(matches!(e, Error::InvalidInput(_)), "得到 {e:?}");
    }

    #[test]
    fn nan_matrix_is_rejected() {
        let m: Mat3 = [[f32::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let d = fake(m, 1000);
        assert!(ColorPlan::from_decoded(&d).is_err());
    }

    #[test]
    fn describe_lists_all_nine_components() {
        let d = fake(mat3::IDENTITY, 1000);
        let p = ColorPlan::from_decoded(&d).unwrap();
        let s = p.describe();
        assert_eq!(s.matches('[').count(), 3, "{s}");
        assert_eq!(s.matches(']').count(), 3, "{s}");
    }

    #[test]
    fn apply_all_matches_per_pixel() {
        let m: Mat3 = [[1.1, -0.1, 0.0], [0.0, 1.0, 0.05], [-0.02, 0.0, 0.9]];
        let d = fake(m, 4095);
        let p = ColorPlan::from_decoded(&d).unwrap();
        let px: Vec<u16> = vec![0, 100, 4095, 2048, 512, 77, 4095, 4095, 4095];
        let all = p.apply_all(&px);
        assert_eq!(all.len(), 3);
        for (i, c) in px.chunks_exact(3).enumerate() {
            assert_eq!(all[i], p.apply_u16([c[0], c[1], c[2]]));
        }
    }
}
