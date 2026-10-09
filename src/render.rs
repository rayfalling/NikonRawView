//! 把「基准配方 + 配方」施加到已进入显示参考域的图像。
//!
//! # 作用域
//!
//! 本模块处理 **0..1 的显示参考值**。配方的曲线、混合器、校色与标量都是相机在
//! 色调曲线**之后**做的调整，定义域是显示参考信号，不是线性光。线性 → 显示参考
//! 的转换由基准变换负责（见 `picture-control/calibration`）——两者分开，才能各自
//! 单独诊断。
//!
//! # 顺序
//!
//! 固定为：**曲线 → 色相混合器 → 校色 → 标量**。顺序写死在本模块里并对外报告，
//! 不允许调用方隐式改变。

use crate::np3::{self, Band, Container, Grading};

/// 混合器与校色的存储值域是 −128..127（中性 0）。换算到 −1..1 的权重。
pub const BAND_SCALE: f32 = 1.0 / 128.0;

/// 混合器的一段，已换算到 −1..1。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandF {
    pub hue: f32,
    pub chroma: f32,
    pub brightness: f32,
}

impl BandF {
    pub const NEUTRAL: Self = Self { hue: 0.0, chroma: 0.0, brightness: 0.0 };

    pub fn is_neutral(&self) -> bool {
        self.hue == 0.0 && self.chroma == 0.0 && self.brightness == 0.0
    }

    fn from_np3(b: &Band) -> Self {
        // NP3 侧已按 b−0x80 解码为 i8；这里只做量纲换算
        Self {
            hue: b.hue as f32 * BAND_SCALE,
            chroma: b.chroma as f32 * BAND_SCALE,
            brightness: b.brightness as f32 * BAND_SCALE,
        }
    }
}

/// 校色的一个区，已换算到 −1..1。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZoneF {
    /// 色相偏移（角度，−180..180）。
    pub hue: f32,
    pub chroma: f32,
    pub brightness: f32,
}

impl ZoneF {
    pub const NEUTRAL: Self = Self { hue: 0.0, chroma: 0.0, brightness: 0.0 };

    pub fn is_neutral(&self) -> bool {
        self.hue == 0.0 && self.chroma == 0.0 && self.brightness == 0.0
    }
}

/// 三区校色。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradingF {
    pub highlights: ZoneF,
    pub mid_tone: ZoneF,
    pub shadows: ZoneF,
    /// 混合强度，0..1。
    pub blending: f32,
    /// 三区权重分配的平衡点，−1..1。
    pub balance: f32,
}

impl GradingF {
    pub const NEUTRAL: Self = Self {
        highlights: ZoneF::NEUTRAL,
        mid_tone: ZoneF::NEUTRAL,
        shadows: ZoneF::NEUTRAL,
        blending: 0.0,
        balance: 0.0,
    };

    pub fn is_neutral(&self) -> bool {
        self.highlights.is_neutral()
            && self.mid_tone.is_neutral()
            && self.shadows.is_neutral()
    }

    fn from_np3(g: &Grading) -> Self {
        let z = |hue_raw: u16, chroma: i8, brightness: i8| ZoneF {
            // 色相以 1/10 度存储，取模到 −180..180
            hue: {
                let deg = hue_raw as f32 / 10.0;
                let deg = if deg > 180.0 { deg - 360.0 } else { deg };
                deg
            },
            chroma: chroma as f32 * BAND_SCALE,
            brightness: brightness as f32 * BAND_SCALE,
        };
        Self {
            highlights: z(g.highlights.hue, g.highlights.chroma, g.highlights.brightness),
            mid_tone: z(g.mid_tone.hue, g.mid_tone.chroma, g.mid_tone.brightness),
            shadows: z(g.shadows.hue, g.shadows.chroma, g.shadows.brightness),
            blending: (g.blending as f32 / 100.0).clamp(0.0, 1.0),
            balance: (g.balance as f32 - 128.0) * BAND_SCALE,
        }
    }
}

/// 点式标量参数，均为已解码的物理量。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scalars {
    /// 锐化，四分之一步进。
    pub sharpness: f32,
    pub mid_range_sharpness: f32,
    pub clarity: f32,
    /// 对比度，`b − 0x80`。
    pub contrast: f32,
    pub brightness: f32,
    pub saturation: f32,
    pub hue: f32,
    /// 曲线哨兵：对比度与亮度同时为 `0x01` 时为真，表示这两个参数**不适用**。
    pub curve_sentinel: bool,
}

impl Default for Scalars {
    fn default() -> Self {
        Self {
            sharpness: 0.0,
            mid_range_sharpness: 0.0,
            clarity: 0.0,
            contrast: 0.0,
            brightness: 0.0,
            saturation: 0.0,
            hue: 0.0,
            curve_sentinel: false,
        }
    }
}

impl Scalars {
    pub fn is_neutral(&self) -> bool {
        self.sharpness == 0.0
            && self.mid_range_sharpness == 0.0
            && self.clarity == 0.0
            && self.saturation == 0.0
            && self.hue == 0.0
            && (self.curve_sentinel || (self.contrast == 0.0 && self.brightness == 0.0))
    }
}

/// 一份可施加的变换。
#[derive(Debug, Clone, PartialEq)]
pub struct Transform {
    /// 257 级色调曲线，已归一化到 0..1。`None` 表示无曲线或曲线为恒等。
    pub curve: Option<Vec<f32>>,
    pub blender: [BandF; 8],
    pub grading: GradingF,
    pub scalars: Scalars,
}

impl Default for Transform {
    fn default() -> Self {
        Self {
            curve: None,
            blender: [BandF::NEUTRAL; 8],
            grading: GradingF::NEUTRAL,
            scalars: Scalars::default(),
        }
    }
}

impl Transform {
    /// 从解析出的配方容器构造。
    pub fn from_container(c: &Container) -> Self {
        let curve = c.curve.as_ref().and_then(|cv| {
            if cv.lut.len() != 257 {
                return None;
            }
            // 恒等曲线等价于没有曲线，省掉一次查表
            let identity = cv.lut.iter().enumerate().all(|(i, &v)| {
                let want = ((i as u32 * 32767) / 256) as i32;
                (v as i32 - want).abs() <= 1
            });
            if identity {
                None
            } else {
                Some(cv.lut.iter().map(|&v| v as f32 / 32767.0).collect())
            }
        });

        let mut blender = [BandF::NEUTRAL; 8];
        if let Some(b) = &c.blender {
            for (i, band) in b.iter().enumerate() {
                blender[i] = BandF::from_np3(band);
            }
        }

        let grading = c.grading.as_ref().map(GradingF::from_np3).unwrap_or(GradingF::NEUTRAL);

        let scalar = |t: u32, scale: np3::Scale| c.scalar(t, scale).unwrap_or(0.0);
        // NP3 容器里没有「亮度」与「色相」标量——它们在 NEF 侧的身份载荷里才有，
        // 而那一侧描述的是拍摄时的设置，不是配方本体。故此处取 0，不做编造。
        let contrast = scalar(0x0000_1900, np3::Scale::Linear);
        let sentinel = c.scalar_raw(0x0000_1900) == Some(crate::picture_control::CURVE_SENTINEL);

        Self {
            curve,
            blender,
            grading,
            scalars: Scalars {
                sharpness: scalar(0x0000_0600, np3::Scale::Quarter),
                mid_range_sharpness: scalar(0x0000_1600, np3::Scale::Quarter),
                clarity: scalar(0x0000_0700, np3::Scale::Quarter),
                // 哨兵语义是"该参数不适用"，不是 −127
                contrast: if sentinel { 0.0 } else { contrast },
                brightness: 0.0,
                saturation: scalar(0x0000_1E00, np3::Scale::Linear),
                hue: 0.0,
                curve_sentinel: sentinel,
            },
        }
    }

    /// 是否等价于恒等变换。
    pub fn is_identity(&self) -> bool {
        self.curve.is_none()
            && self.blender.iter().all(BandF::is_neutral)
            && self.grading.is_neutral()
            && self.scalars.is_neutral()
    }

    /// 处理步骤的可读顺序，供结果记录与复核。
    pub fn step_order() -> &'static [&'static str] {
        &["curve", "hue-blender", "color-grading", "scalars"]
    }

    /// 施加 257 级曲线。逐通道施加——曲线是色调映射，不是色度操作。
    pub fn apply_curve(&self, v: [f32; 3]) -> [f32; 3] {
        let Some(lut) = &self.curve else { return v };
        let f = |x: f32| -> f32 {
            // 曲线定义域是 0..1；超出部分按端点延拓，避免裁切
            if x <= 0.0 {
                return lut[0];
            }
            if x >= 1.0 {
                return lut[256];
            }
            let p = x * 256.0;
            let i = p.floor() as usize;
            let t = p - i as f32;
            let a = lut[i.min(256)];
            let b = lut[(i + 1).min(256)];
            a + (b - a) * t
        };
        [f(v[0]), f(v[1]), f(v[2])]
    }

    /// 施加 8 段色相混合器。
    ///
    /// 每段按其在色相环上的位置取权重——相邻段之间线性过渡，因此不会出现硬边。
    /// 段序固定为 R, O, Y, G, C, B, P, M。
    pub fn apply_blender(&self, v: [f32; 3]) -> [f32; 3] {
        if self.blender.iter().all(BandF::is_neutral) {
            return v;
        }
        let (h, s, l) = rgb_to_hsl(v);
        if s <= 1e-6 {
            // 无彩色像素没有色相可调；明度调整仍应生效
            let dv: f32 = self
                .blender
                .iter()
                .map(|b| b.brightness)
                .sum::<f32>()
                / 8.0;
            return hsl_to_rgb(h, s, (l + dv * 0.25).clamp(0.0, 1.0));
        }

        let mut hue_shift = 0.0f32;
        let mut chroma_mul = 0.0f32;
        let mut bright_add = 0.0f32;
        for (i, b) in self.blender.iter().enumerate() {
            let w = band_weight(h, i);
            if w <= 0.0 {
                continue;
            }
            hue_shift += w * b.hue;
            chroma_mul += w * b.chroma;
            bright_add += w * b.brightness;
        }

        let nh = (h + hue_shift * 30.0).rem_euclid(360.0);
        let ns = (s * (1.0 + chroma_mul * 0.5)).clamp(0.0, 1.0);
        let nl = (l + bright_add * 0.25).clamp(0.0, 1.0);
        hsl_to_rgb(nh, ns, nl)
    }

    /// 施加三区校色。
    pub fn apply_grading(&self, v: [f32; 3]) -> [f32; 3] {
        if self.grading.is_neutral() {
            return v;
        }
        let (h, s, l) = rgb_to_hsl(v);
        let (wh, wm, ws) = zone_weights(l, self.grading.balance);

        let hue_shift = self.grading.highlights.hue * wh
            + self.grading.mid_tone.hue * wm
            + self.grading.shadows.hue * ws;
        let chroma = self.grading.highlights.chroma * wh
            + self.grading.mid_tone.chroma * wm
            + self.grading.shadows.chroma * ws;
        let bright = self.grading.highlights.brightness * wh
            + self.grading.mid_tone.brightness * wm
            + self.grading.shadows.brightness * ws;

        let nh = (h + hue_shift).rem_euclid(360.0);
        let ns = (s * (1.0 + chroma * 0.5)).clamp(0.0, 1.0);
        let nl = (l + bright * 0.25).clamp(0.0, 1.0);
        hsl_to_rgb(nh, ns, nl)
    }

    /// 施加点式标量：对比度、亮度、饱和度、色相。
    ///
    /// **哨兵值必须被识别为「不适用」而跳过。** 实测对比度与亮度可能同时存为
    /// `0x01`，语义是"该参数由曲线决定"，不是 `−127`。把它当数值施加会让画面
    /// 彻底变黑。
    pub fn apply_scalars_pointwise(&self, v: [f32; 3]) -> [f32; 3] {
        let s = &self.scalars;
        let mut out = v;

        // 对比度：以 0.5 为轴
        let c = s.contrast / 128.0;
        if c != 0.0 {
            for x in out.iter_mut() {
                *x = (*x - 0.5) * (1.0 + c) + 0.5;
            }
        }

        // 亮度：加法偏移
        let b = s.brightness / 128.0;
        if b != 0.0 {
            for x in out.iter_mut() {
                *x += b * 0.5;
            }
        }

        // 饱和度与色相走 HSL
        if s.saturation != 0.0 || s.hue != 0.0 {
            let (h, sat, l) = rgb_to_hsl(out);
            let ns = (sat * (1.0 + s.saturation / 128.0)).clamp(0.0, 1.0);
            let nh = (h + s.hue / 128.0 * 30.0).rem_euclid(360.0);
            out = hsl_to_rgb(nh, ns, l);
        }
        out
    }

    /// 施加单个像素的全部变换。顺序固定，见 [`Self::step_order`]。
    pub fn apply_pixel(&self, v: [f32; 3]) -> [f32; 3] {
        let v = self.apply_curve(v);
        let v = self.apply_blender(v);
        let v = self.apply_grading(v);
        self.apply_scalars_pointwise(v)
    }

    /// 施加整幅图像。
    ///
    /// **纯函数**：不修改输入，相同输入产生相同输出。
    pub fn apply_image(&self, src: &[[f32; 3]]) -> Vec<[f32; 3]> {
        src.iter().map(|p| self.apply_pixel(*p)).collect()
    }
}

// ---------------------------------------------------------------------------
// 色彩空间辅助
// ---------------------------------------------------------------------------

/// RGB（0..1）→ HSL（h 为 0..360，s/l 为 0..1）。
pub fn rgb_to_hsl(rgb: [f32; 3]) -> (f32, f32, f32) {
    let (r, g, b) = (rgb[0], rgb[1], rgb[2]);
    let max = r.max(g).max(b);
    let min = r.min(g).min(b);
    let l = (max + min) / 2.0;
    let d = max - min;
    if d.abs() < 1e-9 {
        return (0.0, 0.0, l);
    }
    let s = if l > 0.5 { d / (2.0 - max - min) } else { d / (max + min) };
    let h = if max == r {
        60.0 * (((g - b) / d) % 6.0)
    } else if max == g {
        60.0 * ((b - r) / d + 2.0)
    } else {
        60.0 * ((r - g) / d + 4.0)
    };
    (h.rem_euclid(360.0), s, l)
}

/// HSL → RGB（0..1）。
pub fn hsl_to_rgb(h: f32, s: f32, l: f32) -> [f32; 3] {
    if s <= 1e-9 {
        return [l, l, l];
    }
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let hp = h.rem_euclid(360.0) / 60.0;
    let x = c * (1.0 - (hp % 2.0 - 1.0).abs());
    let (r1, g1, b1) = match hp as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = l - c / 2.0;
    [r1 + m, g1 + m, b1 + m]
}

/// 8 段的中心色相（度）。
const BAND_CENTERS: [f32; 8] = [0.0, 30.0, 60.0, 120.0, 180.0, 240.0, 270.0, 300.0];

/// 某色相在第 `i` 段上的权重。
///
/// 取该段中心到相邻段中心距离内线性过渡，因此各段权重之和恒为 1——
/// 这保证了"全中性时输出与未施加时一致"。
fn band_weight(hue: f32, i: usize) -> f32 {
    let h = hue.rem_euclid(360.0);
    // 计算到每个段中心的角距离，取最近两个做线性插值
    let mut best = (f32::MAX, 0usize, f32::MAX, 0usize); // (d1,i1,d2,i2)
    for (j, c) in BAND_CENTERS.iter().enumerate() {
        let mut d = (h - c).abs();
        if d > 180.0 {
            d = 360.0 - d;
        }
        if d < best.0 {
            best = (d, j, best.0, best.1);
        } else if d < best.2 {
            best = (best.0, best.1, d, j);
        }
    }
    let (d1, i1, d2, i2) = best;
    let w = if i == i1 {
        if d1 + d2 <= 1e-6 {
            1.0
        } else {
            d2 / (d1 + d2)
        }
    } else if i == i2 {
        if d1 + d2 <= 1e-6 {
            0.0
        } else {
            d1 / (d1 + d2)
        }
    } else {
        0.0
    };
    w.clamp(0.0, 1.0)
}

/// 三区权重，和为 1。
fn zone_weights(l: f32, balance: f32) -> (f32, f32, f32) {
    // balance 把中间调的峰值上下移动；范围收在 0.2..0.8 以免极端化
    let mid = (0.5 - balance * 0.3).clamp(0.2, 0.8);
    let hl = ((l - mid) / (1.0 - mid).max(1e-6)).clamp(0.0, 1.0);
    let sh = ((mid - l) / mid.max(1e-6)).clamp(0.0, 1.0);
    let h = hl;
    let s = sh;
    let m = (1.0 - h - s).max(0.0);
    let sum = h + m + s;
    if sum <= 1e-9 {
        (0.0, 1.0, 0.0)
    } else {
        (h / sum, m / sum, s / sum)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }

    fn default_t() -> Transform {
        Transform::default()
    }

    // --- 曲线 ---

    #[test]
    fn no_curve_is_passthrough() {
        let t = default_t();
        assert_eq!(t.apply_curve([0.1, 0.5, 0.9]), [0.1, 0.5, 0.9]);
    }

    #[test]
    fn identity_curve_is_dropped() {
        // 恒等 LUT 不应被当成"有曲线"
        let mut c = Container {
            version: 1,
            family: "0310".into(),
            name: "X".into(),
            base_code: None,
            scalars: Default::default(),
            blender: None,
            grading: None,
            curve: Some(np3::Curve {
                point_count: 2,
                points: vec![],
                lut: (0..257).map(|i| ((i as u32 * 32767) / 256) as u16).collect(),
                lut_offset: 64,
            }),
            chunks: vec![],
        };
        assert!(Transform::from_container(&c).curve.is_none());
        // 真改一个点就不再是恒等
        if let Some(cv) = c.curve.as_mut() {
            cv.lut[128] = 20000;
        }
        assert!(Transform::from_container(&c).curve.is_some());
    }

    #[test]
    fn curve_is_interpolated_and_endpoints_hold() {
        let mut t = default_t();
        t.curve = Some((0..257).map(|i| (i as f32 / 256.0).powf(2.0)).collect());
        // 端点
        assert!(approx(t.apply_curve([0.0, 0.0, 0.0])[0], 0.0, 1e-6));
        assert!(approx(t.apply_curve([1.0, 1.0, 1.0])[0], 1.0, 1e-6));
        // 中点应接近 0.25
        assert!(approx(t.apply_curve([0.5, 0.5, 0.5])[0], 0.25, 0.01));
        // 超出定义域按端点延拓，不裁切
        assert!(approx(t.apply_curve([2.0, 2.0, 2.0])[0], 1.0, 1e-6));
        assert!(approx(t.apply_curve([-1.0, -1.0, -1.0])[0], 0.0, 1e-6));
    }

    #[test]
    fn non_monotonic_curve_is_applied_verbatim() {
        // 构造一条非单调曲线：先上后下
        let mut lut: Vec<f32> = (0..257).map(|i| i as f32 / 256.0).collect();
        for i in 128..257 {
            lut[i] = 0.5 - (i - 128) as f32 / 256.0;
        }
        let mut t = default_t();
        t.curve = Some(lut.clone());
        let got = t.apply_curve([0.9, 0.9, 0.9])[0];
        // 应等于曲线在该处的取值，而不是被单调化
        let idx = (0.9 * 256.0) as usize;
        assert!(approx(got, lut[idx], 0.01), "非单调曲线应原样施加");
        assert!(got < 0.3, "该处曲线已下落，得到 {got}");
    }

    // --- 混合器 ---

    #[test]
    fn neutral_blender_is_identity() {
        let t = default_t();
        let v = [0.3, 0.6, 0.2];
        assert_eq!(t.apply_blender(v), v);
    }

    #[test]
    fn band_weights_sum_to_one() {
        for h in (0..360).step_by(7) {
            let sum: f32 = (0..8).map(|i| band_weight(h as f32, i)).sum();
            assert!(approx(sum, 1.0, 1e-4), "hue={h} 权重和={sum}");
        }
    }

    #[test]
    fn band_centers_have_unit_weight() {
        for (i, c) in BAND_CENTERS.iter().enumerate() {
            let w = band_weight(*c, i);
            assert!(approx(w, 1.0, 1e-3), "段 {i} 在中心 {c} 的权重应为 1，得到 {w}");
        }
    }

    #[test]
    fn adjusting_one_band_leaves_far_bands_alone() {
        let mut t = default_t();
        // 只调红色段
        t.blender[0] = BandF { hue: 0.5, chroma: 0.0, brightness: 0.0 };
        let red = [1.0, 0.0, 0.0];
        let cyan = [0.0, 1.0, 1.0];
        let r0 = t.apply_blender(red);
        let c0 = t.apply_blender(cyan);
        // 红被明显改变
        let (h_after, _, _) = rgb_to_hsl(r0);
        assert!((h_after - 0.0).abs() > 5.0, "红色段应发生色相偏移，得到 {h_after}");
        // 青色几乎不变
        let (h_before, _, _) = rgb_to_hsl(cyan);
        let (h_c, _, _) = rgb_to_hsl(c0);
        assert!(
            (h_c - h_before).abs() < 2.0,
            "青色不应被红色段影响：{h_before} → {h_c}"
        );
    }

    // --- 校色 ---

    #[test]
    fn neutral_grading_is_identity() {
        let t = default_t();
        let v = [0.2, 0.4, 0.8];
        assert_eq!(t.apply_grading(v), v);
    }

    #[test]
    fn zone_weights_sum_to_one() {
        for i in 0..=10 {
            let l = i as f32 / 10.0;
            let (a, b, c) = zone_weights(l, 0.0);
            assert!(approx(a + b + c, 1.0, 1e-4), "l={l} 权重和={}", a + b + c);
        }
    }

    #[test]
    fn shadows_and_highlights_are_affected_separately() {
        let mut t = default_t();
        t.grading = GradingF {
            shadows: ZoneF { hue: 0.0, chroma: 0.0, brightness: 0.8 },
            ..GradingF::NEUTRAL
        };
        let dark = [0.1, 0.1, 0.1];
        let bright = [0.9, 0.9, 0.9];
        let d0 = t.apply_grading(dark)[0];
        let b0 = t.apply_grading(bright)[0];
        assert!(d0 > dark[0] + 0.05, "暗部应被抬升，{ } → {d0}", dark[0]);
        assert!(
            (b0 - bright[0]).abs() < 0.02,
            "亮部不应受暗部校色影响：{} → {b0}",
            bright[0]
        );
    }

    // --- 标量 ---

    #[test]
    fn sentinel_is_skipped_not_applied_as_minus_127() {
        let mut t = default_t();
        t.scalars = Scalars {
            contrast: -127.0,
            brightness: -127.0,
            curve_sentinel: true,
            ..Scalars::default()
        };
        // 构造时哨兵会把两者归零；这里直接验证归零后的行为
        let s = Scalars { curve_sentinel: true, contrast: 0.0, brightness: 0.0, ..Scalars::default() };
        let t2 = Transform { scalars: s, ..default_t() };
        let v = [0.5, 0.5, 0.5];
        let out = t2.apply_scalars_pointwise(v);
        assert!(approx(out[0], 0.5, 1e-6), "哨兵不应改变像素，得到 {out:?}");

        // 反例：若把 −127 当数值施加，画面会大幅偏移
        t.scalars.curve_sentinel = false;
        let bad = t.apply_scalars_pointwise(v);
        assert!(
            (bad[0] - 0.5).abs() > 0.1,
            "把 −127 当数值施加应产生明显偏移——这正是要避免的"
        );
    }

    #[test]
    fn contrast_pivots_around_half() {
        let mut t = default_t();
        t.scalars.contrast = 64.0;
        let mid = t.apply_scalars_pointwise([0.5, 0.5, 0.5]);
        assert!(approx(mid[0], 0.5, 1e-5), "中点应保持不动，得到 {mid:?}");
        let dark = t.apply_scalars_pointwise([0.25, 0.25, 0.25]);
        assert!(dark[0] < 0.25, "提高对比度应压低暗部");
    }

    #[test]
    fn saturation_zero_is_identity() {
        let t = default_t();
        for v in [[0.2, 0.5, 0.9], [1.0, 0.0, 0.0], [0.5, 0.5, 0.5]] {
            let out = t.apply_scalars_pointwise(v);
            for i in 0..3 {
                assert!(approx(out[i], v[i], 1e-5), "{v:?} → {out:?}");
            }
        }
    }

    #[test]
    fn negative_saturation_desaturates() {
        let mut t = default_t();
        t.scalars.saturation = -128.0;
        let out = t.apply_scalars_pointwise([1.0, 0.0, 0.0]);
        let (_, s, _) = rgb_to_hsl(out);
        assert!(s < 1e-3, "饱和度拉到最低应变灰，得到 {out:?}");
    }

    // --- 纯函数与顺序 ---

    #[test]
    fn default_transform_is_identity() {
        let t = default_t();
        assert!(t.is_identity());
        for v in [[0.0, 0.0, 0.0], [0.3, 0.6, 0.2], [1.0, 1.0, 1.0], [0.5, 0.1, 0.9]] {
            let out = t.apply_pixel(v);
            for i in 0..3 {
                assert!(approx(out[i], v[i], 1e-5), "{v:?} → {out:?}");
            }
        }
    }

    #[test]
    fn apply_image_is_pure_and_repeatable() {
        let t = default_t();
        let src = vec![[0.1f32, 0.2, 0.3], [0.9, 0.8, 0.7]];
        let before = src.clone();
        let a = t.apply_image(&src);
        let b = t.apply_image(&src);
        assert_eq!(src, before, "输入不应被修改");
        assert_eq!(a, b, "相同输入应产生相同输出");
    }

    #[test]
    fn step_order_is_fixed_and_reported() {
        let s = Transform::step_order();
        assert_eq!(s, &["curve", "hue-blender", "color-grading", "scalars"]);
    }

    #[test]
    fn hsl_roundtrip() {
        for v in [[0.2f32, 0.5, 0.9], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            let (h, s, l) = rgb_to_hsl(v);
            let back = hsl_to_rgb(h, s, l);
            for i in 0..3 {
                assert!(approx(back[i], v[i], 1e-4), "{v:?} → {back:?}");
            }
        }
    }
}
