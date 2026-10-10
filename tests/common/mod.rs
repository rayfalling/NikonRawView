//! 测试夹具：可提交的**期望值**，以及合成文件的构造器。
//!
//! 原始素材（NEF/JPEG 样张）不入库，仓库内只保留期望值。需要真实样本的测试
//! 在样本缺席时打印提示并跳过，而不是失败——这样克隆仓库后 `cargo test` 仍然干净。

#![allow(dead_code)]

use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// 期望值（fixture）
// ---------------------------------------------------------------------------

/// `simple/DSC_4143.NEF` 的实测期望值。
pub struct Expected {
    pub file: &'static str,
    pub name: &'static str,
    pub base: &'static str,
    pub payload_len: usize,
    pub curve_sentinel: bool,
    pub sharpness: f32,
    pub clarity: f32,
    pub saturation: i16,
    pub hue: i16,
}

pub const DSC4143: Expected = Expected {
    file: "DSC_4143.NEF",
    name: "LINKS-Nature",
    base: "NEUTRAL",
    payload_len: 108,
    curve_sentinel: true,
    sharpness: 0.5,
    clarity: 0.5,
    saturation: 5,
    hue: -8,
};

/// 本机样本目录（默认 `simple/`，可用 `NIKONRAWVIEW_SAMPLES` 覆盖）。
pub fn samples_dir() -> PathBuf {
    std::env::var_os("NIKONRAWVIEW_SAMPLES")
        .map(PathBuf::from)
        .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("simple"))
}

/// 取一个样本文件；不存在时返回 `None`（调用方据此跳过）。
pub fn sample(name: &str) -> Option<Vec<u8>> {
    std::fs::read(samples_dir().join(name)).ok()
}

/// 样本缺席时的统一提示。
pub fn skip(test: &str, name: &str) {
    eprintln!(
        "跳过 {test}：未找到本机样本 {}（把样张放到 {} 或设置 NIKONRAWVIEW_SAMPLES）",
        name,
        samples_dir().display()
    );
}

// ---------------------------------------------------------------------------
// 诊断脚手架
// ---------------------------------------------------------------------------

/// 一个诊断样本。
#[derive(Debug, Clone, Copy)]
pub struct DiagSample {
    /// 我们的**原始**线性工作空间值（未经基准变换）。
    pub ours_linear: [f32; 3],
    /// 施加基准变换之后的值。
    pub ours_transformed: [f32; 3],
    /// 参考导出的线性工作空间值。
    pub theirs_linear: [f32; 3],
    /// 像素在图像中的位置（裁剪后坐标），便于回看具体区域。
    pub xy: (usize, usize),
}

/// 拟合产物 + 评估样本。
pub struct FittedFixture {
    pub transform: nikonrawview::transform::BaseTransform,
    pub eval: Vec<DiagSample>,
    pub fit_count: usize,
    /// 参考图与解码图的尺寸，便于诊断报告写清上下文。
    pub size: (usize, usize),
}

impl DiagSample {
    /// 亮度（Rec.709 权重，作用于我们的原始线性值）。
    pub fn lum_ours(&self) -> f32 {
        0.2126 * self.ours_linear[0] + 0.7152 * self.ours_linear[1] + 0.0722 * self.ours_linear[2]
    }
    /// 参考的亮度。
    pub fn lum_theirs(&self) -> f32 {
        0.2126 * self.theirs_linear[0] + 0.7152 * self.theirs_linear[1] + 0.0722 * self.theirs_linear[2]
    }
    /// 该样本的 ΔE00（变换后 vs 参考）。
    pub fn delta_e(&self) -> f64 {
        nikonrawview::deltae::delta_e_prophoto(self.ours_transformed, self.theirs_linear)
    }
}

/// 构建 `DSC_0001` 的拟合产物与留出评估集。
///
/// **拟合与评估用不同的像素子集**——即便正式划分（5.6）尚未落地，也不该拿"背下来的
/// 答案"来评估自己。
///
/// 三条诊断共用这一份，避免各自重跑一遍解码（一张 45 MP 的 NEF 解码约 20 秒，
/// 参考 TIF 有 260 MB）。
pub fn fitted_neutral_fixture() -> Option<FittedFixture> {
    fitted_fixture_for("DSC_0001", true)
}

/// 构建 `stem` 的拟合产物与留出评估集。
///
/// `do_fit` 为 false 时只收集样本、不拟合——用于「拿 A 图拟合出的变换去套 B 图」的
/// 跨图验证。返回值里的 `transform` 此时是恒等，调用方应自行提供变换。
pub fn fitted_fixture_for(stem: &str, do_fit: bool) -> Option<FittedFixture> {
    let dir = samples_dir();
    let nef = dir.join(format!("{stem}.NEF"));
    let tif = dir.join(format!("{stem}.TIF"));
    if !nef.is_file() || !tif.is_file() {
        skip("诊断脚手架", &format!("{stem}.NEF / {stem}.TIF"));
        return None;
    }

    let tif_data = std::fs::read(&tif).ok()?;
    let img = nikonrawview::fit::read_rgb16(&tif_data).ok()?;
    let plan = nikonrawview::icc::plan_for(&tif_data).ok()?;
    let opts = nikonrawview::libraw::Options {
        demosaic: nikonrawview::libraw::Demosaic::Dht,
        user_mul: None,
    };
    let dec = nikonrawview::libraw::decode_working_space(&nef, &opts).ok()?;
    let model = nikonrawview::camera::read_model(&nef)?;
    let entry = nikonrawview::camera::lookup(&model)?;
    let mg = if dec.rotated { entry.margins.rotated() } else { entry.margins };
    let (cw, ch) = mg.effective(dec.width, dec.height)?;
    let theirs = nikonrawview::fit::reference_to_working(&img, &plan);

    let mut fit_samples = Vec::new();
    let mut heldout = Vec::new();
    let mut k = 0usize;
    for y in (0..ch).step_by(13) {
        for x in (0..cw).step_by(13) {
            let j = x + y * img.width;
            let Some(p) = dec.at(x + mg.left, y + mg.top) else { continue };
            let ours = [
                p[0] as f32 / 65535.0,
                p[1] as f32 / 65535.0,
                p[2] as f32 / 65535.0,
            ];
            k += 1;
            if k % 2 == 0 {
                fit_samples.push(nikonrawview::fit::Sample { ours, theirs: theirs[j] });
            } else {
                heldout.push((ours, theirs[j], (x, y)));
            }
        }
    }

    let curve = nikonrawview::fit::fit_curve(&fit_samples, 128);
    let edge = 17;
    let lut = nikonrawview::fit::fit_lut(&fit_samples, &curve, edge);
    let transform = if do_fit {
        nikonrawview::transform::BaseTransform::from_fit(
            0x0000,
            "NEUTRAL",
            vec![format!("simple/{stem}.NEF + simple/{stem}.TIF")],
            curve,
            edge,
            lut,
        )
    } else {
        // 只收样本，不拟合——跨图验证时调用方会提供在别处拟合出的变换
        nikonrawview::transform::BaseTransform::new(0x0000, "未拟合")
    };

    let eval = heldout
        .into_iter()
        .map(|(o, t, xy)| DiagSample {
            ours_linear: o,
            ours_transformed: transform.apply(o),
            theirs_linear: t,
            xy,
        })
        .collect();

    Some(FittedFixture {
        transform,
        eval,
        fit_count: fit_samples.len(),
        size: (cw, ch),
    })
}

/// 计算一组 ΔE00 的中位数 / P95 / 最大值。
pub fn quantiles(values: &mut [f64]) -> (f64, f64, f64) {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if values.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let at = |q: f64| {
        let rank = (q * values.len() as f64).ceil().max(1.0) as usize;
        values[rank.saturating_sub(1).min(values.len() - 1)]
    };
    (at(0.5), at(0.95), *values.last().unwrap())
}

// ---------------------------------------------------------------------------
// 合成文件构造器
// ---------------------------------------------------------------------------

const TIFF_BASE_IFD0: usize = 8;

/// 构造一个自洽的 TIFF（NEF 形态）：IFD0 → ExifIFD → MakerNote → tag 0x0023。
///
/// 所有偏移都相对 TIFF 头，因此把同一份字节放进 JPEG 的 APP1 里时，
/// **只有加上基址才能解析成功**——这正是要回归验证的性质。
pub fn build_tiff(pc_payload: &[u8], marker_bytes: &[u8; 8]) -> Vec<u8> {
    let ifd0_off = TIFF_BASE_IFD0;
    let exif_off = ifd0_off + 2 + 12 + 4;
    let mn_off = exif_off + 2 + 12 + 4;
    // MakerNote: "Nikon\0" + 2 字节版本 + 2 字节填充，其后内嵌 TIFF 头
    let mn_tiff = mn_off + 10;
    let mn_ifd = mn_tiff + 8;
    let marker_off = mn_ifd + 2 + 12 + 4;
    let pc_off = marker_off + marker_bytes.len();

    let mut d = vec![0u8; pc_off + pc_payload.len()];
    d[0..2].copy_from_slice(b"II");
    d[2..4].copy_from_slice(&42u16.to_le_bytes());
    d[4..8].copy_from_slice(&(ifd0_off as u32).to_le_bytes());

    // IFD0：一个条目，指向 ExifIFD
    d[ifd0_off..ifd0_off + 2].copy_from_slice(&1u16.to_le_bytes());
    put_entry(&mut d, ifd0_off + 2, 0x8769, 4, 1, exif_off as u32);

    // ExifIFD：一个条目，指向 MakerNote
    d[exif_off..exif_off + 2].copy_from_slice(&1u16.to_le_bytes());
    put_entry(&mut d, exif_off + 2, 0x927C, 7, marker_bytes.len() as u32, mn_off as u32);

    // MakerNote 头 + 内嵌 TIFF 头
    d[mn_off..mn_off + 6].copy_from_slice(b"Nikon\x00");
    d[mn_off + 6] = 0x02;
    d[mn_off + 7] = 0x11;
    d[mn_tiff..mn_tiff + 2].copy_from_slice(b"II");
    d[mn_tiff + 2..mn_tiff + 4].copy_from_slice(&42u16.to_le_bytes());
    d[mn_tiff + 4..mn_tiff + 8].copy_from_slice(&8u32.to_le_bytes());

    // MakerNote IFD：一个条目 tag 0x0023
    d[mn_ifd..mn_ifd + 2].copy_from_slice(&1u16.to_le_bytes());
    put_entry(&mut d, mn_ifd + 2, 0x0023, 7, pc_payload.len() as u32,
              (pc_off - mn_tiff) as u32);

    // 载荷之前放上那段**可变的**前导字节，以验证实现不依赖它
    d[marker_off..marker_off + marker_bytes.len()].copy_from_slice(marker_bytes);
    d[pc_off..].copy_from_slice(pc_payload);
    d
}

fn put_entry(d: &mut [u8], at: usize, tag: u16, typ: u16, count: u32, value: u32) {
    d[at..at + 2].copy_from_slice(&tag.to_le_bytes());
    d[at + 2..at + 4].copy_from_slice(&typ.to_le_bytes());
    d[at + 4..at + 8].copy_from_slice(&count.to_le_bytes());
    d[at + 8..at + 12].copy_from_slice(&value.to_le_bytes());
    // 条目之后一个 IFD 的偏移
    let next = at + 12;
    if next + 4 <= d.len() {
        d[next..next + 4].copy_from_slice(&0u32.to_le_bytes());
    }
}

/// 把一份 TIFF 包进 JPEG 的 Exif APP1 段。
pub fn wrap_jpeg(tiff: &[u8]) -> Vec<u8> {
    let mut app1 = Vec::new();
    app1.extend_from_slice(b"Exif\x00\x00");
    app1.extend_from_slice(tiff);

    let mut j = Vec::new();
    j.extend_from_slice(&[0xFF, 0xD8]); // SOI
    j.extend_from_slice(&[0xFF, 0xE1]); // APP1
    j.extend_from_slice(&((app1.len() + 2) as u16).to_be_bytes());
    j.extend_from_slice(&app1);
    j.extend_from_slice(&[0xFF, 0xD9]); // EOI
    j
}

/// 构造一份符合 PictureControl3 布局的合成载荷。
pub fn synth_pc(name: &str, base: &str, contrast: u8, brightness: u8) -> Vec<u8> {
    let mut b = vec![0u8; 108];
    b[..8].copy_from_slice(b"03100310");
    b[8..8 + name.len()].copy_from_slice(name.as_bytes());
    b[28..28 + base.len()].copy_from_slice(base.as_bytes());
    b[54] = 1;
    b[55] = 0xFF;
    b[57] = 0x82;
    b[59] = 0x80;
    b[61] = 0x82;
    b[63] = contrast;
    b[65] = brightness;
    b[67] = 0x85;
    b[69] = 0x78;
    b[71] = 0xFF;
    b[72] = 0xFF;
    b[73] = 0xFF;
    b
}

/// 那段实测会变化的前导字节（`01 00 ?? 00 02 00 00 00`）的两种取值。
pub const MARKER_V1: [u8; 8] = [0x01, 0x00, 0x01, 0x00, 0x02, 0x00, 0x00, 0x00];
pub const MARKER_V3: [u8; 8] = [0x01, 0x00, 0x03, 0x00, 0x02, 0x00, 0x00, 0x00];
