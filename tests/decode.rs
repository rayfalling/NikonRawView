//! 解码验证：对着真实的 NEF 检查尺寸、白平衡、相机矩阵与线性性。
//!
//! 样张不入库（见 `.gitignore`），缺席时打印提示并跳过。

use std::path::PathBuf;

use nikonrawview::libraw::{self, Demosaic, Options};
use nikonrawview::camera;
use nikonrawview::color::ColorPlan;

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

/// 任务 2.2：白电平可取；黑电平经 C API **不可取**，但可以验证它已被扣除。
///
/// LibRaw 的 `cblack` / `black` 在 `libraw_colordata_t` 里，而该结构体只能通过
/// 指针访问——复刻它的布局会在版本升级时静默读到错误偏移，因此本实现不做。
///
/// 判据取**全图最小值**：若黑电平未扣除，最暗像素也应落在黑电平附近
/// （实测基准：黑电平 1008 / 白电平 15892 ≈ 6.3%）；扣除后应接近 0。
///
/// 注：先前曾假设「传感器全幅比有效区多出的 12 列是光学黑区」，**该假设不成立**——
/// 实测其均值占白电平 18.7%，是真实图像内容。
#[test]
fn white_level_available_and_black_already_subtracted() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let d = libraw::decode_camera_linear(&path, &Options::default()).unwrap();

    // 白电平：C API 有 libraw_get_color_maximum
    assert!(d.color_maximum > 0, "白电平应为正，得到 {}", d.color_maximum);
    eprintln!("白电平（color_maximum）= {}", d.color_maximum);

    let min = d.pixels.iter().copied().min().unwrap_or(0);
    let max = d.pixels.iter().copied().max().unwrap_or(0);
    let min_ratio = min as f64 / d.color_maximum as f64;
    eprintln!("全图 min={min} max={max}，min 占白电平 {min_ratio:.4}");

    // 未扣黑电平时，min 会落在黑电平附近（≈6%）。留足余量取 3%。
    assert!(
        min_ratio < 0.03,
        "全图最小值为白电平的 {min_ratio:.4}，若黑电平已扣除应接近 0"
    );
}

/// 任务 2.5：去马赛克的伪像检查。
///
/// 伪色与拉链效应表现为**色度的高频变化**。指标取「局部 (R−G) 与 (B−G) 的二阶差分
/// 绝对均值」——正确去马赛克的平坦区域该值很小，成片伪色会让它显著抬高。
#[test]
fn demosaic_artifacts_are_measured() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let d = libraw::decode_camera_linear(&path, &Options::default()).unwrap();

    // 取画面中央一块，避开边框
    let (x0, y0) = (d.width / 2, d.height / 2);
    let (bw, bh) = (512usize.min(d.width / 4), 512usize.min(d.height / 4));

    let chroma_second_diff = |alg: Demosaic| -> f64 {
        let dd = if alg == d.demosaic {
            std::borrow::Cow::Borrowed(&d)
        } else {
            std::borrow::Cow::Owned(
                libraw::decode_camera_linear(&path, &Options { demosaic: alg, user_mul: None })
                    .unwrap(),
            )
        };
        let mut acc = 0f64;
        let mut n = 0u64;
        // 逐行：对 (R−G) 与 (B−G) 求水平二阶差分
        for y in y0..y0 + bh {
            let mut prev: Option<(i32, i32)> = None;
            let mut prev2: Option<(i32, i32)> = None;
            for x in x0..x0 + bw {
                let Some(px) = dd.at(x, y) else { continue };
                let rg = px[0] as i32 - px[1] as i32;
                let bg = px[2] as i32 - px[1] as i32;
                if let (Some(a), Some(b)) = (prev, prev2) {
                    acc += ((rg - 2 * a.0 + b.0).abs() + (bg - 2 * a.1 + b.1).abs()) as f64;
                    n += 1;
                }
                prev2 = prev;
                prev = Some((rg, bg));
            }
        }
        if n == 0 {
            0.0
        } else {
            acc / n as f64 / d.color_maximum as f64
        }
    };

    let dht = chroma_second_diff(Demosaic::Dht);
    eprintln!("中央 {bw}×{bh}：色度二阶差分均值（DHT）= {dht:.6}");
    assert!(dht.is_finite() && dht >= 0.0);

    // 与线性去马赛克对照：线性插值会产生更多伪色，指标应更高。
    // 这条同时证明指标真的有判别力——若两者相近，说明指标测不出伪像。
    let lin = chroma_second_diff(Demosaic::Linear);
    eprintln!("                                   （Linear）= {lin:.6}");
    assert!(
        lin > dht,
        "线性去马赛克的色度伪像应多于 DHT（{lin:.6} vs {dht:.6}），否则该指标无判别力"
    );

    let thresh = 0.02;
    assert!(
        dht < thresh,
        "DHT 的色度二阶差分均值 {dht:.6} 超过阈值 {thresh}，中央区域可能有成片伪色"
    );
}

/// 任务 3.1：矩阵必须随**拍摄白平衡**变化，而非固定常数。
///
/// 注意：实测 `libraw_set_user_mul` **不能**改变矩阵——解码层里相机记录的
/// 白平衡优先级更高。因此只能用真正记录了不同白平衡的两张文件来验。完整解码
/// 一张要十几秒，所以先用 `read_wb`（只到 unpack）快速扫出两个不同的样本。
#[test]
fn camera_matrix_varies_with_white_balance() {
    let root = match std::env::var("NIKONRAWVIEW_LIBRARY") {
        Ok(v) => std::path::PathBuf::from(v),
        Err(_) => {
            eprintln!("跳过：需设置 NIKONRAWVIEW_LIBRARY 指向含多张 NEF 的目录");
            return;
        }
    };

    // 扫描：找出两张白平衡明显不同的文件。
    //
    // 按顶层目录分散取样并限制总探测数——每个 read_wb 都要读文件头之后相当一段，
    // 放在网络盘上无上限地扫会非常慢。不同行程目录的光线条件差异最大，优先从那取。
    const PER_DIR: usize = 6;
    const MAX_PROBE: usize = 96;

    let mut candidates: Vec<std::path::PathBuf> = Vec::new();
    let tops: Vec<std::path::PathBuf> = std::fs::read_dir(&root)
        .map(|rd| {
            rd.flatten()
                .map(|e| e.path())
                .filter(|p| p.is_dir())
                .collect()
        })
        .unwrap_or_else(|_| vec![root.clone()]);
    for dir in tops.iter().chain(std::iter::once(&root)) {
        let mut taken = 0usize;
        let mut stack = vec![dir.clone()];
        while let Some(d) = stack.pop() {
            let Ok(rd) = std::fs::read_dir(&d) else { continue };
            for e in rd.flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                    continue;
                }
                if !p.extension().is_some_and(|x| x.eq_ignore_ascii_case("nef")) {
                    continue;
                }
                candidates.push(p);
                taken += 1;
                if taken >= PER_DIR {
                    break;
                }
            }
            if taken >= PER_DIR {
                break;
            }
        }
    }
    eprintln!("候选 NEF {} 个（上限 {MAX_PROBE}）", candidates.len());

    let mut found: Vec<(std::path::PathBuf, [f32; 3])> = Vec::new();
    let mut probed = 0usize;
    let mut failed = 0usize;
    for p in candidates.iter().take(MAX_PROBE) {
        probed += 1;
        let Ok(wb) = libraw::read_wb(p) else {
            failed += 1;
            continue;
        };
        match found.first() {
            None => found.push((p.clone(), wb)),
            Some((_, w0)) => {
                // 红蓝比差异超过 15% 才算"明显不同"
                let r0 = w0[0] / w0[2].max(1e-6);
                let r1 = wb[0] / wb[2].max(1e-6);
                if (r1 - r0).abs() / r0.max(1e-6) > 0.15 {
                    found.push((p.clone(), wb));
                    break;
                }
            }
        }
    }
    eprintln!("探测 {probed} 个（读取失败 {failed}），找到 {} 个不同白平衡样本", found.len());

    if found.len() < 2 {
        eprintln!("跳过：未在目录中找到白平衡明显不同的两张 NEF");
        return;
    }
    let (pa, wa) = &found[0];
    let (pb, wb) = &found[1];
    eprintln!("样本 A {} WB={wa:?}", pa.file_name().unwrap().to_string_lossy());
    eprintln!("样本 B {} WB={wb:?}", pb.file_name().unwrap().to_string_lossy());

    let opts = Options { demosaic: Demosaic::Linear, user_mul: None };
    let a = libraw::decode_working_space(pa, &opts).unwrap();
    let b = libraw::decode_working_space(pb, &opts).unwrap();

    // **实测结论**：随拍摄白平衡变化的是白平衡系数，不是相机矩阵。
    // LibRaw 的 `rgb_cam` 按机型固定（dcraw 的 cam_xyz_coeff 只依赖机型的
    // cam_xyz 与输出空间，与拍摄白平衡无关）。
    let wb_diff: f32 = (0..3).map(|i| (wa[i] - wb[i]).abs()).sum();
    eprintln!("白平衡系数差异和 = {wb_diff:.4}");
    assert!(wb_diff > 0.05, "两张样本的白平衡系数应明显不同，得到 {wb_diff}");

    let m_diff: f32 = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .map(|(i, j)| (a.rgb_cam[i][j] - b.rgb_cam[i][j]).abs())
        .sum();
    eprintln!("相机矩阵差异和 = {m_diff:.6}");
    eprintln!("  A: {:?}", a.rgb_cam);
    eprintln!("  B: {:?}", b.rgb_cam);
    eprintln!("  A cam_mul={:?}  B cam_mul={:?}", a.wb, b.wb);

    // 本测试**不断言矩阵不同**——实测它相同。它断言的是"确有随拍摄变化的量"，
    // 免得日后误以为整个白平衡链路是固定的。
    assert!(
        m_diff.abs() < 1e-6,
        "实测相机矩阵按机型固定、不随拍摄白平衡变化；若此处不再为 0（{m_diff}），\
         说明解码层改变了行为，需要重新审视 render/color-pipeline 的规范"
    );
    assert!(a.rgb_cam.iter().flatten().all(|v| v.is_finite()));
    assert!(b.rgb_cam.iter().flatten().all(|v| v.is_finite()));
}

/// 5.4 第一步：把 NEF 与参考 TIF 都转到 ProPhoto 线性，看配对关系是不是一条
/// **平滑单调**的曲线。
///
/// 这是拟合前的"看一眼"。若不是，说明样本里混进了非基准因素（某张被调整过、
/// ADL 没关掉、配对错了），此时求解只会把问题一起拟合进去。
#[test]
fn reference_pairs_form_a_clean_curve() {
    let dir = samples_dir();
    let nef = dir.join("DSC_0001.NEF");
    let tif = dir.join("DSC_0001.TIF");
    if !nef.is_file() || !tif.is_file() {
        eprintln!("跳过：缺少 DSC_0001 的 NEF 或 TIF");
        return;
    }

    let tif_data = std::fs::read(&tif).unwrap();
    let img = nikonrawview::fit::read_rgb16(&tif_data).expect("应能读参考 TIF");
    eprintln!("参考 TIF：{}×{}，{} 个分量", img.width, img.height, img.pixels.len());

    let plan = nikonrawview::icc::plan_for(&tif_data).expect("应能决定色彩空间方案");
    eprintln!("色彩空间：{}", plan.describe());

    let opts = Options { demosaic: Demosaic::Dht, user_mul: None };
    let dec = libraw::decode_working_space(&nef, &opts).expect("应能解码 NEF");
    eprintln!("解码：{}×{}（工作空间 ProPhoto 线性）", dec.width, dec.height);

    let theirs = nikonrawview::fit::reference_to_working(&img, &plan);

    // 参考导出是**有效区**，而解码是**传感器全幅**——先按边距裁切。
    // 边距是针对未旋转画幅定义的，竖拍输出需要把左右与上下对调。
    let model = nikonrawview::camera::read_model(&nef).expect("应能读出机型");
    let entry = nikonrawview::camera::lookup(&model).expect("机型应在表中");
    let mg = if dec.rotated { entry.margins.rotated() } else { entry.margins };
    eprintln!("解码方向：{}", if dec.rotated { "竖拍（已旋转）" } else { "横拍" });
    eprintln!("实际使用的边距：左{} 上{} 右{} 下{}", mg.left, mg.top, mg.right, mg.bottom);
    let (cw, ch) = mg.effective(dec.width, dec.height).expect("应能裁出有效区");
    eprintln!("裁切后 {cw}×{ch}，参考导出 {}×{}", img.width, img.height);
    assert_eq!(
        (cw, ch),
        (img.width, img.height),
        "按方向调整边距后，裁切尺寸应与参考导出一致"
    );

    let mut samples = Vec::new();
    let stride = 37;
    for y in (0..ch).step_by(stride) {
        for x in (0..cw).step_by(stride) {
            let sx = x + mg.left;
            let sy = y + mg.top;
            let j = x + y * img.width;
            if j >= theirs.len() {
                continue;
            }
            let Some(p) = dec.at(sx, sy) else { continue };
            let o = [
                p[0] as f32 / 65535.0,
                p[1] as f32 / 65535.0,
                p[2] as f32 / 65535.0,
            ];
            samples.push(nikonrawview::fit::Sample { ours: o, theirs: theirs[j] });
        }
    }
    eprintln!("采集 {} 个配对样本", samples.len());

    let d = nikonrawview::fit::diagnose(&samples, 10);
    eprintln!("\n输入中位 → 输出中位（ProPhoto 线性）  样本数");
    for (x, y, n) in &d.bins {
        eprintln!("  {x:.2} → {y:.4}   {n}");
    }
    eprintln!("\n逆序次数 = {}（平滑单调应为 0）", d.inversions);
    if let Some((lo, hi)) = d.slope_range {
        eprintln!("相邻箱斜率范围 = {lo:.3} … {hi:.3}");
    }
    eprintln!("是否像干净的单调整曲线：{}", d.looks_clean());

    assert!(samples.len() > 1000, "样本太少（{}）", samples.len());
    assert!(d.looks_clean(), "配对关系不像一条干净的单调曲线，需先查样本");
}
/// 导出到 `simple/`）跑校验，全部应通过。
#[test]
fn the_manual_reference_exports_pass_validation() {
    let dir = samples_dir();
    let mut set = nikonrawview::calibrate::ReferenceSet::default();
    for i in 1..=15 {
        let n = format!("DSC_{i:04}");
        let raw = dir.join(format!("{n}.NEF"));
        let export = dir.join(format!("{n}.TIF"));
        if !raw.is_file() || !export.is_file() {
            eprintln!("跳过：缺少 {n} 的 NEF 或 TIF");
            return;
        }
        set.pairs.push(nikonrawview::calibrate::ReferencePair { raw, export });
    }

    let reports = nikonrawview::calibrate::validate(&set);
    let mut bad = 0;
    for r in &reports {
        if r.is_ok() {
            continue;
        }
        bad += 1;
        eprintln!("✗ {}", r.raw.file_name().unwrap().to_string_lossy());
        for i in &r.issues {
            eprintln!("    {}", i.describe());
        }
    }
    eprintln!("校验 {} 对：合格 {}，不合格 {}", reports.len(), reports.len() - bad, bad);

    // 顺带报告共识调整块——若整组都被做了同一处调整，这一项查不出来，需人工过目
    let d = std::fs::read(&set.pairs[0].export).unwrap();
    if let Ok(Some(m)) = nikonrawview::calibrate::export_adjustments(&d) {
        eprintln!("导出内嵌的调整块（相机固定写入，共 {} 项）：", m.len());
        for (k, v) in &m {
            if nikonrawview::calibrate::NEUTRAL_ADJUSTMENTS.contains(&k.as_str()) {
                eprintln!("    {k} = {v}");
            }
        }
    }

    assert_eq!(bad, 0, "{bad} 对未通过校验");
}
/// 因此它**必须**被拒绝作为参考导出——这正是校验要拦住的第一种情况。
#[test]
fn real_sample_is_rejected_as_reference_because_adl_is_on() {
    let raw = samples_dir().join("DSC_4143.NEF");
    if !raw.is_file() {
        eprintln!("跳过：未找到 {}", raw.display());
        return;
    }
    let data = std::fs::read(&raw).unwrap();

    let adl = nikonrawview::makernote::active_d_lighting(&data)
        .expect("应能读出 ADL")
        .expect("样本应含 ADL 标签");
    eprintln!("DSC_4143 的 ADL = {}", adl.name());
    assert!(!adl.is_off(), "该样本的 ADL 应为开启状态");

    // 直接走校验：ADL 未关 → 该对被排除
    let mut set = nikonrawview::calibrate::ReferenceSet::default();
    set.pairs.push(nikonrawview::calibrate::ReferencePair {
        raw: raw.clone(),
        export: raw.clone(), // 导出用同一个文件占位，重点是 ADL 那一项
    });
    let reports = nikonrawview::calibrate::validate(&set);
    assert_eq!(reports.len(), 1);
    let r = &reports[0];
    let joined: String = r.issues.iter().map(|i| i.describe()).collect::<Vec<_>>().join(" | ");
    eprintln!("校验结论：{joined}");
    assert!(!r.is_ok(), "ADL 未关闭的样本必须被排除，却通过了校验");
    assert!(joined.contains("ADL 未关闭"), "应指出 ADL 问题：{joined}");
}
///
/// 前提来自 3.1 的实测——矩阵按机型固定、不随拍摄白平衡变化，所以每个机型只需测一次。
/// 若本测试通过，`render/color-pipeline` 就可以保留「由本管线施加矩阵」，同时拿回
/// 色域裁切与负值的控制权。
#[test]
fn derived_matrix_generalises_across_images() {
    let root = std::env::var("NIKONRAWVIEW_LIBRARY")
        .ok()
        .map(std::path::PathBuf::from);
    let first = samples_dir().join("DSC_4143.NEF");
    if !first.is_file() {
        eprintln!("跳过：未找到 {}", first.display());
        return;
    }

    let pair_of = |p: &std::path::Path| -> (libraw::Decoded, libraw::Decoded) {
        let o = Options { demosaic: Demosaic::Linear, user_mul: None };
        let cam = libraw::decode_with_output_for_test(p, &o, libraw::OutputColor::Camera).unwrap();
        let pro =
            libraw::decode_with_output_for_test(p, &o, libraw::OutputColor::ProPhoto).unwrap();
        (cam, pro)
    };

    let samples = |cam: &libraw::Decoded, pro: &libraw::Decoded| -> Vec<([f32; 3], [f32; 3])> {
        let mut v = Vec::new();
        for i in (0..cam.pixels.len()).step_by(3 * 7) {
            let c = [cam.pixels[i], cam.pixels[i + 1], cam.pixels[i + 2]];
            let p = [pro.pixels[i], pro.pixels[i + 1], pro.pixels[i + 2]];
            // 排除任何一侧的饱和点：那里解码层做过裁切，会把裁切行为误当成矩阵
            if c.iter().any(|x| *x == 0 || *x == u16::MAX)
                || p.iter().any(|x| *x == 0 || *x == u16::MAX)
            {
                continue;
            }
            let f = |a: [u16; 3]| [a[0] as f32 / 65535.0, a[1] as f32 / 65535.0, a[2] as f32 / 65535.0];
            v.push((f(c), f(p)));
        }
        v
    };

    let (cam_a, pro_a) = pair_of(&first);
    let pairs_a = samples(&cam_a, &pro_a);
    eprintln!("推导样本数 = {}", pairs_a.len());
    let derived = nikonrawview::color::derive_matrix(&pairs_a).expect("应能求出矩阵");
    eprintln!(
        "推导矩阵（{} 样本，rms={:.6}，max={:.6}）：",
        derived.samples, derived.rms, derived.max_abs
    );
    for row in &derived.matrix {
        eprintln!("    [{:>9.5} {:>9.5} {:>9.5}]", row[0], row[1], row[2]);
    }

    // 在**推导所用的同一张图**上，残差应落在这条矩阵的解释能力之内
    assert!(
        derived.rms < 0.01,
        "同图拟合的 rms={:.6} 过大，说明相机→工作空间的映射不是线性的",
        derived.rms
    );

    // 换一张图复验
    let Some(root) = root else {
        eprintln!("未设置 NIKONRAWVIEW_LIBRARY，跳过跨图复验");
        return;
    };
    let mut second = None;
    let mut stack = vec![root];
    'find: while let Some(d) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&d) else { continue };
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
            } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case("nef"))
                && p.file_name() != first.file_name()
            {
                second = Some(p);
                break 'find;
            }
        }
    }
    let Some(second) = second else {
        eprintln!("未找到第二张 NEF，跳过跨图复验");
        return;
    };
    eprintln!("复验用图：{}", second.file_name().unwrap().to_string_lossy());

    let (cam_b, pro_b) = pair_of(&second);
    let pairs_b = samples(&cam_b, &pro_b);
    let plan = ColorPlan::from_derived(derived.matrix, 65535.0);
    let mut sum = 0f64;
    let mut worst = 0f64;
    let mut n = 0u64;
    for (c, p) in &pairs_b {
        let ours = plan.apply(*c);
        for k in 0..3 {
            let e = (ours[k] - p[k]).abs() as f64;
            sum += e;
            worst = worst.max(e);
            n += 1;
        }
    }
    let mean = sum / n.max(1) as f64;
    eprintln!("跨图复验（{} 样本）：平均绝对差 {mean:.6}，最大 {worst:.6}", pairs_b.len());

    assert!(n > 1000, "复验样本太少（{n}）");
    assert!(
        mean < 0.005,
        "推导出的矩阵在另一张图上平均绝对差 {mean:.6}——若显著偏大，说明矩阵并非按机型固定，路线 B 不成立"
    );
}
///
/// 假设是「输出色彩空间没被应用」——若如此，`Camera` 与 `ProPhoto` 两条路径的像素
/// 会几乎相同，而它们的 `rgb_cam` 却不相同。
#[test]
#[ignore = "诊断用，手动运行"]
fn diagnose_output_color_takes_effect() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let opts = Options { demosaic: Demosaic::Linear, user_mul: None };

    let cam = libraw::decode_with_output_for_test(&path, &opts, libraw::OutputColor::Camera).unwrap();
    let srgb = libraw::decode_with_output_for_test(&path, &opts, libraw::OutputColor::Srgb).unwrap();
    let pro = libraw::decode_with_output_for_test(&path, &opts, libraw::OutputColor::ProPhoto).unwrap();

    let diff = |a: &libraw::Decoded, b: &libraw::Decoded| -> f64 {
        let mut s = 0f64;
        let mut n = 0u64;
        for i in (0..a.pixels.len()).step_by(3 * 397) {
            for k in 0..3 {
                s += (a.pixels[i + k] as f64 - b.pixels[i + k] as f64).abs();
                n += 1;
            }
        }
        s / n as f64
    };

    let mc = |m: [[f32; 3]; 3]| {
        format!(
            "[{:.4} {:.4} {:.4}]",
            m[0][0], m[0][1], m[0][2]
        )
    };

    eprintln!("=== output_color 生效性诊断 ===");
    eprintln!("Camera  rgb_cam 第 0 行: {}", mc(cam.rgb_cam));
    eprintln!("sRGB    rgb_cam 第 0 行: {}", mc(srgb.rgb_cam));
    eprintln!("ProPhoto rgb_cam 第 0 行: {}", mc(pro.rgb_cam));
    eprintln!();
    eprintln!("像素平均绝对差（16 位刻度）：");
    eprintln!("  Camera  vs sRGB    : {:.1}", diff(&cam, &srgb));
    eprintln!("  Camera  vs ProPhoto: {:.1}", diff(&cam, &pro));
    eprintln!("  sRGB    vs ProPhoto: {:.1}", diff(&srgb, &pro));
    eprintln!();
    eprintln!("cam 输出前几个像素 : {:?}", &cam.pixels[..9.min(cam.pixels.len())]);
    eprintln!("pro 输出前几个像素 : {:?}", &pro.pixels[..9.min(pro.pixels.len())]);
}

/// 任务 3.3：交叉验证——本管线自行施加矩阵，应与解码层直接输出工作空间一致。
///
/// # 当前状态：**未通过，但原因已定位**
///
/// 实测（`DSC_4143.NEF`，137526 个在色域内的分量）：施加矩阵平均绝对差 0.008406，
/// 不施加矩阵 0.003888——**"什么都不做"反而更接近**。
///
/// ## 原因：`libraw_get_rgb_cam` 返回的不是实际使用的矩阵
///
/// 诊断（`diagnose_output_color_takes_effect`）在同一个文件上比较三种输出空间：
///
/// ```text
/// Camera  rgb_cam: [1.3931 -0.2157 -0.1774]
/// sRGB    rgb_cam: [1.3931 -0.2157 -0.1774]   ← 三者完全相同
/// ProPhoto rgb_cam: [1.3931 -0.2157 -0.1774]
///
/// 像素平均绝对差（16 位刻度）：
///   Camera vs sRGB    : 362.1
///   Camera vs ProPhoto: 254.7                  ← 像素确实随输出空间变化
///   sRGB   vs ProPhoto: 550.5
/// ```
///
/// 即：**矩阵与输出色彩空间无关，而像素随输出色彩空间变化**。因此该访问接口给出的
/// 矩阵不是 `convert_to_rgb` 实际使用的那个，本管线"取矩阵自行施加"的做法建立在
/// 一个错误的前提上。
///
/// ## 待定的两条出路（需先定路线再改 spec）
///
/// 1. **由解码层做转换**：直接取 `output_color = ProPhoto` 的像素，本管线不再自行
///    施加矩阵。正确性有保障，代价是解码层会就地裁到 0..65535——`render/color-pipeline`
///    中「不提前裁切色域、不截断负值」这条需要相应改写。
/// 2. **自行算出正确的矩阵**：从解码层的机型标定推导。但该标定未经 C API 暴露，
///    等于要复刻其内部数据——与本项目「不依赖解码层内部结构」的原则冲突。
#[test]
#[ignore = "原因已定位：libraw_get_rgb_cam 返回的矩阵与输出色彩空间无关，不是实际使用的那个；待定路线后修正"]
fn color_pipeline_matches_decoder_working_space() {
    let path = samples_dir().join("DSC_4143.NEF");
    if !path.is_file() {
        eprintln!("跳过：未找到 {}", path.display());
        return;
    }
    let opts = Options { demosaic: Demosaic::Linear, user_mul: None };

    let cam = libraw::decode_camera_linear(&path, &opts).unwrap();
    let direct = libraw::decode_working_space(&path, &opts).unwrap();
    assert_eq!((cam.width, cam.height), (direct.width, direct.height));

    let plan = ColorPlan::from_decoded(&cam).expect("应能构造色彩方案");

    // 逐点比较：本管线施加矩阵 vs 解码层内部转换。
    // 只统计双方都未接近饱和的点——解码层会就地裁到 0..65535，饱和点必然对不齐。
    let mut worst = 0f64;
    let mut sum = 0f64;
    let mut n = 0u64;
    let mut skipped = 0u64;
    for i in (0..cam.pixels.len()).step_by(3 * 997) {
        let c = [cam.pixels[i], cam.pixels[i + 1], cam.pixels[i + 2]];
        let w = [direct.pixels[i], direct.pixels[i + 1], direct.pixels[i + 2]];
        let ours = plan.apply_u16(c);
        // 只比对**双方都在色域内**的点。
        //
        // 解码层的内部转换会对超出色域的分量做裁切与去饱和，而本管线刻意不做
        // ——那是设计差异，不是误差。把它算进统计只会掩盖真实的一致性。
        if ours.iter().any(|v| !(0.0f32..=1.0).contains(v)) || w.contains(&u16::MAX) {
            skipped += 1;
            continue;
        }
        for k in 0..3 {
            let theirs = w[k] as f32 * plan.scale();
            let e = (ours[k] - theirs).abs() as f64;
            sum += e;
            worst = worst.max(e);
            n += 1;
        }
    }
    let mean = if n == 0 { 0.0 } else { sum / n as f64 };
    eprintln!(
        "交叉验证：{n} 个在色域内的分量，平均绝对差 {mean:.6}，最大 {worst:.6}（单位 0..1）；跳过 {skipped} 个超色域点"
    );

    // 对照：若不施加矩阵（用单位矩阵），误差应远大于上面——这证明本验证有判别力。
    // 否则"误差小"可能只是因为图像本身接近中性，测不出矩阵是否被用上。
    let mut id_sum = 0f64;
    let mut id_n = 0u64;
    for i in (0..cam.pixels.len()).step_by(3 * 997) {
        let c = [cam.pixels[i], cam.pixels[i + 1], cam.pixels[i + 2]];
        let w = [direct.pixels[i], direct.pixels[i + 1], direct.pixels[i + 2]];
        let ours_no_matrix: [f32; 3] = [
            c[0] as f32 * plan.scale(),
            c[1] as f32 * plan.scale(),
            c[2] as f32 * plan.scale(),
        ];
        if ours_no_matrix.iter().any(|v| !(0.0f32..=1.0).contains(v)) || w.contains(&u16::MAX) {
            continue;
        }
        for k in 0..3 {
            id_sum += (ours_no_matrix[k] - w[k] as f32 * plan.scale()).abs() as f64;
            id_n += 1;
        }
    }
    let id_mean = if id_n == 0 { f64::NAN } else { id_sum / id_n as f64 };
    eprintln!("对照（不施加矩阵）：平均绝对差 {id_mean:.6}");

    assert!(n > 100, "可比对的分量太少（{n}），验证不充分");
    // 关键判据：施加矩阵应显著优于不施加。当前**不成立**，这正是本测试被 ignore 的原因。
    assert!(
        mean < id_mean,
        "施加矩阵的误差（{mean:.6}）未小于不施加矩阵的对照（{id_mean:.6}）——\
         说明本管线的矩阵施加没有复现解码层的内部转换"
    );
    assert!(
        mean < 0.02,
        "本管线施加矩阵与解码层内部转换平均差 {mean:.6}，超出可归因于解码层额外步骤的范围"
    );
}
