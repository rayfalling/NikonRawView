//! 3×3 矩阵与色彩变换的公共数学。
//!
//! 本模块只做矩阵运算，不含色彩策略——策略在 [`crate::color`]。

/// 行主序 3×3 矩阵。
pub type Mat3 = [[f32; 3]; 3];

pub const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

pub fn mul(a: Mat3, b: Mat3) -> Mat3 {
    let mut o = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    o
}

pub fn mul_vec(m: Mat3, v: [f32; 3]) -> [f32; 3] {
    let mut o = [0.0f32; 3];
    for i in 0..3 {
        o[i] = (0..3).map(|k| m[i][k] * v[k]).sum();
    }
    o
}

/// 行列式。
pub fn det(m: Mat3) -> f32 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// 逆矩阵。奇异时返回 `None`——调用方必须处理，不得静默返回单位阵。
pub fn inverse(m: Mat3) -> Option<Mat3> {
    let d = det(m);
    if d.abs() < 1e-12 {
        return None;
    }
    let inv = 1.0 / d;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv,
        ],
    ])
}

/// 由长度为 9 的行主序切片构造矩阵。
pub fn from_slice(s: &[f32]) -> Option<Mat3> {
    if s.len() < 9 {
        return None;
    }
    Some([[s[0], s[1], s[2]], [s[3], s[4], s[5]], [s[6], s[7], s[8]]])
}

/// 逐元素线性插值：`t = 0` 取 `a`，`t = 1` 取 `b`。
pub fn lerp(a: Mat3, b: Mat3, t: f32) -> Mat3 {
    let mut o = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = a[i][j] * (1.0 - t) + b[i][j] * t;
        }
    }
    o
}

/// Bradford 色适应矩阵（XYZ 空间）。
pub const BRADFORD: Mat3 = [
    [0.8951, 0.2664, -0.1614],
    [-0.7502, 1.7135, 0.0367],
    [0.0389, -0.0685, 1.0296],
];

/// Bradford 的逆矩阵（直接给出，避免每次求逆带来的误差累积）。
pub const BRADFORD_INV: Mat3 = [
    [0.986_992_9, -0.147_054_3, 0.159_962_7],
    [0.432_305_3, 0.518_360_3, 0.049_291_2],
    [-0.008_528_7, 0.040_042_8, 0.968_486_7],
];

/// XYZ(D50) → ProPhoto RGB 线性。
///
/// ProPhoto（ROMM）的原色与 D50 白点，取自规范。
pub const XYZ_TO_PROPHOTO: Mat3 = [
    [1.345_943_3, -0.255_607_5, -0.051_111_8],
    [-0.544_598_9, 1.508_167_3, 0.020_535_1],
    [0.0, 0.0, 1.211_812_8],
];

/// 计算把 `src_white` 适应到 `dst_white` 的 Bradford 色适应矩阵。
///
/// 两者均为 XYZ 三刺激值；内部按 Bradford 锥体响应空间对角缩放。
pub fn bradford_adapt(src_white: [f32; 3], dst_white: [f32; 3]) -> Mat3 {
    let s = mul_vec(BRADFORD, src_white);
    let d = mul_vec(BRADFORD, dst_white);
    let scale: Mat3 = [
        [d[0] / s[0], 0.0, 0.0],
        [0.0, d[1] / s[1], 0.0],
        [0.0, 0.0, d[2] / s[2]],
    ];
    mul(BRADFORD_INV, mul(scale, BRADFORD))
}

/// D50 白点的 XYZ 三刺激值（ProPhoto 的参考白）。
pub const WHITE_D50: [f32; 3] = [0.964_295_8, 1.0, 0.825_104_6];

/// 归一化到 Y = 1。
pub fn normalize_y(v: [f32; 3]) -> [f32; 3] {
    if v[1].abs() < 1e-12 {
        return v;
    }
    [v[0] / v[1], 1.0, v[2] / v[1]]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn inverse_roundtrips() {
        let m: Mat3 = [[0.5, 0.2, 0.1], [0.1, 0.9, 0.3], [0.2, 0.1, 0.7]];
        let i = inverse(m).unwrap();
        let p = mul(m, i);
        for r in 0..3 {
            for c in 0..3 {
                let want = if r == c { 1.0 } else { 0.0 };
                assert!(approx(p[r][c], want, 1e-4), "{p:?}");
            }
        }
    }

    #[test]
    fn singular_inverse_is_none() {
        let m: Mat3 = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 1.0, 1.0]];
        assert!(inverse(m).is_none(), "奇异矩阵必须返回 None 而非静默给单位阵");
    }

    #[test]
    fn bradford_identity_on_same_white() {
        let m = bradford_adapt(WHITE_D50, WHITE_D50);
        for r in 0..3 {
            for c in 0..3 {
                let want = if r == c { 1.0 } else { 0.0 };
                assert!(approx(m[r][c], want, 1e-5), "{m:?}");
            }
        }
    }

    #[test]
    fn bradford_maps_src_white_to_dst_white() {
        let d65 = [0.950_455_9, 1.0, 1.088_754_2];
        let m = bradford_adapt(d65, WHITE_D50);
        let got = mul_vec(m, d65);
        let want = normalize_y(WHITE_D50);
        let got = normalize_y(got);
        for i in 0..3 {
            assert!(approx(got[i], want[i], 1e-3), "{got:?} vs {want:?}");
        }
    }

    #[test]
    fn prophoto_matrix_maps_d50_to_neutral() {
        let rgb = mul_vec(XYZ_TO_PROPHOTO, WHITE_D50);
        assert!(approx(rgb[0] / rgb[1], 1.0, 1e-3), "{rgb:?}");
        assert!(approx(rgb[2] / rgb[1], 1.0, 1e-3), "{rgb:?}");
    }

    #[test]
    fn lerp_endpoints() {
        let a = IDENTITY;
        let b = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        assert_eq!(lerp(a, b, 0.0), a);
        assert_eq!(lerp(a, b, 1.0), b);
        assert_eq!(lerp(a, b, 0.5)[0][0], 1.5);
    }
}
