//! 解码验证：对着真实的 NEF 检查尺寸、白平衡、相机矩阵与线性性。
//!
//! 样张不入库（见 `.gitignore`），缺席时打印提示并跳过。

use std::path::PathBuf;

use nikonrawview::libraw::{self, Demosaic, Options};
use nikonrawview::camera;

fn samples_dir() -> PathBuf {
    std::env::var_os("NIKONRAWVIEW_SAMPLES")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("simple"))
}

#[test]
fn links_and_reports_version() {
    eprintln!("LibRaw {}", libraw::version_string());
    assert!(libraw::version_number() > 0);
}

#[test]
fn decodes_z8_nef_to_linear_camera_rgb() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }

    let t = std::time::Instant::now();
    let d = libraw::decode_camera_linear(&path, &Options::default())
        .expect("应能解码 DSC_4143.NEF");
    let ms = t.elapsed().as_millis();

    eprintln!("解码耗时 {ms} ms");
    eprintln!(
        "  输出 {}×{}   传感器全幅 {}×{}   已裁到有效区={}",
        d.width, d.height, d.raw_width, d.raw_height, d.is_cropped()
    );
    eprintln!("  去马赛克={}", d.demosaic.name());
    eprintln!("  白平衡 R/G/B = {:?}", d.wb);
    eprintln!("  色彩上限 = {}", d.color_maximum);
    eprintln!("  相机 → ProPhoto 矩阵：");
    for row in &d.rgb_cam {
        eprintln!("    [{:>9.5} {:>9.5} {:>9.5}]", row[0], row[1], row[2]);
    }

    // 尺寸：Z8 的传感器全幅为 8280×5520，有效像素区为 8256×5504。
    // LibRaw 默认输出全幅（含光学黑/掩蔽边框），是否已裁到有效区由 is_cropped() 报告。
    assert!(
        d.width == 8280 || d.width == 8256,
        "宽度应为 8280（全幅）或 8256（有效区），得到 {}",
        d.width
    );
    assert!(
        d.height == 5520 || d.height == 5504,
        "高度应为 5520（全幅）或 5504（有效区），得到 {}",
        d.height
    );
    assert_eq!(d.pixels.len(), d.width * d.height * 3);

    // 白平衡：与早前 LibRaw 独立测得的一致（归一化到 G=1）
    let g = d.wb[1];
    assert!(g > 0.0, "绿通道系数应为正");
    let r = d.wb[0] / g;
    let b = d.wb[2] / g;
    eprintln!("  归一化后 R/G/B = {r:.4} / 1.0 / {b:.4}");
    assert!((r - 1.8125).abs() < 0.05, "R 系数应约 1.8125，得到 {r}");
    assert!((b - 1.5703).abs() < 0.05, "B 系数应约 1.5703，得到 {b}");

    // 相机矩阵：不应是单位阵——否则说明取到的不是相机标定
    let mut is_identity = true;
    for (i, row) in d.rgb_cam.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let want = if i == j { 1.0 } else { 0.0 };
            if (v - want).abs() > 1e-3 {
                is_identity = false;
            }
        }
    }
    assert!(!is_identity, "rgb_cam 不应是单位阵");
    assert!(
        d.rgb_cam.iter().flatten().all(|v| v.is_finite()),
        "矩阵分量应全为有限值"
    );
}

/// 线性性：输出不应含 gamma。
///
/// 判据是**相邻档位的比值**——线性数据的分布会被场景决定，但若含 sRGB gamma，
/// 暗部会被显著抬升，中位值相对最大值的比例会明显偏高。
#[test]
fn output_is_linear_not_gamma_encoded() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let d = libraw::decode_camera_linear(&path, &Options::default()).unwrap();

    // 取绿通道（CFA 中数量最多，噪声最小）
    let mut g: Vec<u16> = d.pixels.iter().skip(1).step_by(3).copied().collect();
    g.sort_unstable();
    let max = *g.last().unwrap() as f64;
    assert!(max > 0.0, "不应全黑");
    let median = g[g.len() / 2] as f64;

    let ratio = median / max;
    eprintln!("绿通道 中位/最大 = {ratio:.4}（线性应明显低于 gamma 编码）");

    // 若被 sRGB gamma 编码，中位/最大 会显著偏高。这里给一个宽松但有判别力的上界：
    // 线性下该比值通常在 0.05–0.35；gamma 编码会把中暗部抬到 0.55 以上。
    assert!(
        ratio < 0.5,
        "中位/最大 = {ratio:.4} 过高，输出可能含 gamma 曲线"
    );
}

/// 不同去马赛克算法应产生不同结果，且算法被如实报告。
#[test]
fn demosaic_choice_is_reported_and_effective() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let lin = libraw::decode_camera_linear(&path, &Options { demosaic: Demosaic::Linear, user_mul: None })
        .expect("linear 去马赛克应可用");
    assert_eq!(lin.demosaic, Demosaic::Linear);
    assert_eq!(lin.demosaic.name(), "linear");

    let dht = libraw::decode_camera_linear(&path, &Options::default()).unwrap();
    assert_eq!(dht.demosaic, Demosaic::Dht);

    assert_eq!(lin.pixels.len(), dht.pixels.len());
    assert_ne!(
        lin.pixels, dht.pixels,
        "两种去马赛克算法不应产出逐位相同的结果"
    );
}

/// 可复现：相同参数两次解码逐位相同。
#[test]
fn decoding_is_reproducible() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let o = Options { demosaic: Demosaic::Ahd, user_mul: None };
    let a = libraw::decode_camera_linear(&path, &o).unwrap();
    let b = libraw::decode_camera_linear(&path, &o).unwrap();
    assert_eq!(a.pixels, b.pixels, "两次解码应逐位相同");
    assert_eq!(a.rgb_cam, b.rgb_cam);
    assert_eq!(a.wb, b.wb);
}

/// 裁切：机型从文件读出，边距取自机型表，开关两种状态都要验证。
#[test]
fn crop_toggle_yields_expected_sizes() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }

    let model = camera::read_model(&path).expect("应能从 IFD0 读出机型");
    eprintln!("从 IFD0 读出的机型：{model:?}");
    assert_eq!(camera::lookup(&model).map(|m| m.display), Some("Nikon Z 8"));

    let d = libraw::decode_camera_linear(&path, &Options::default()).unwrap();
    assert_eq!((d.width, d.height), (8280, 5520), "解码层给出传感器全幅");

    // 开关关闭：保留全幅
    let off = camera::decide(&model, d.width, d.height, false);
    assert_eq!(off, camera::CropDecision::Disabled);
    let (px, w, h) = camera::apply(&d.pixels, d.width, d.height, &off);
    assert_eq!((w, h), (8280, 5520));
    assert_eq!(px.len(), d.pixels.len());
    eprintln!("  开关关闭 → {w}×{h}  {}", off.describe());

    // 开关开启：裁到有效像素区
    let on = camera::decide(&model, d.width, d.height, true);
    assert!(matches!(on, camera::CropDecision::Cropped { .. }), "得到 {on:?}");
    let (px, w, h) = camera::apply(&d.pixels, d.width, d.height, &on);
    assert_eq!((w, h), (8256, 5504), "裁切后应为 Z8 的有效像素区");
    assert_eq!(px.len(), w * h * 3);
    eprintln!("  开关开启 → {w}×{h}  {}", on.describe());

    // 裁切后的首像素应等于原图 (12, 8) 处
    let src = d.at(12, 8).unwrap();
    assert_eq!([px[0], px[1], px[2]], src, "裁切起点应为左上角边距处");
}

/// 机型不在边距表中时不猜边距。
#[test]
fn unknown_model_does_not_guess_margins() {
    let d = camera::decide("Canon EOS R5", 8280, 5520, true);
    assert!(matches!(d, camera::CropDecision::UnknownModel { .. }), "得到 {d:?}");
    assert!(d.describe().contains("机型不在边距表中"));
    let px = vec![1u16; 4 * 4 * 3];
    let (out, w, h) = camera::apply(&px, 4, 4, &d);
    assert_eq!((w, h), (4, 4), "未知机型应原样返回而非裁切");
    assert_eq!(out.len(), px.len());
}
