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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libraw::{Decoded, Demosaic};

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
