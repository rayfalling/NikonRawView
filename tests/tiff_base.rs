//! TIFF/IFD 解析与基址处理的回归测试。
//!
//! 这些用例**不依赖本机样张**：合成一份自洽的 TIFF，再把它包进 JPEG 的 Exif
//! APP1 段。同一份字节在两种容器里必须得到**完全相同**的解析结果，而"忘记
//! 加基址"必须失败。

mod common;

use common::*;
use nikonrawview::picture_control::{self, Adjust};
use nikonrawview::tiff::Tiff;

#[test]
fn parses_synthetic_nef() {
    let pc = synth_pc("LINKS-Nature", "NEUTRAL", 0x01, 0x01);
    let tiff = build_tiff(&pc, &MARKER_V3);
    let id = picture_control::read(&tiff).expect("NEF 应可解析");
    assert_eq!(id.name, "LINKS-Nature");
    assert_eq!(id.base, "NEUTRAL");
    assert!(id.is_custom());
}

#[test]
fn parses_synthetic_jpeg_identically() {
    let pc = synth_pc("LINKS-Nature", "NEUTRAL", 0x01, 0x01);
    let tiff = build_tiff(&pc, &MARKER_V3);
    let jpeg = wrap_jpeg(&tiff);

    let a = picture_control::read(&tiff).expect("NEF 应可解析");
    let b = picture_control::read(&jpeg).expect("JPEG 应可解析");
    assert_eq!(a.raw, b.raw, "两种容器的载荷必须逐字节相同");
}

/// 回归：把 JPEG 的字节当成"TIFF 头在偏移 0"来解析必须失败。
///
/// 这正是早前 JPEG 侧 5678 个文件全部解析失败的原因——偏移是相对 APP1 内
/// TIFF 头的，不加基址就会在错误的位置找 IFD。
#[test]
fn forgetting_the_jpeg_base_offset_fails() {
    let pc = synth_pc("X", "NEUTRAL", 0x80, 0x80);
    let jpeg = wrap_jpeg(&build_tiff(&pc, &MARKER_V3));

    assert!(
        Tiff::at(&jpeg, 0).is_err(),
        "以偏移 0 为基址必须失败（JPEG 开头是 FFD8，不是 TIFF 头）"
    );
    let t = Tiff::locate(&jpeg).expect("locate 应找到 Exif APP1 内的 TIFF 头");
    assert!(t.base() > 0, "JPEG 的 TIFF 头基址必须大于 0");
    assert_eq!(&jpeg[t.base()..t.base() + 2], b"II");
}

/// 规格要求：前导字节中那个 u16 会变化，实现 MUST NOT 依赖它。
#[test]
fn variable_leading_bytes_do_not_affect_result() {
    let pc = synth_pc("LINKS-Nature", "NEUTRAL", 0x01, 0x01);
    let a = picture_control::read(&build_tiff(&pc, &MARKER_V1)).unwrap();
    let b = picture_control::read(&build_tiff(&pc, &MARKER_V3)).unwrap();
    assert_eq!(a.name, b.name);
    assert_eq!(a.base, b.base);
    assert_eq!(a.raw, b.raw, "前导字节不同不应改变解析出的载荷");
}

#[test]
fn missing_picture_control_tag_is_explicit() {
    // 构造一个没有 tag 0x0023 的 MakerNote：把条目 tag 改成别的
    let pc = synth_pc("X", "NEUTRAL", 0x80, 0x80);
    let mut tiff = build_tiff(&pc, &MARKER_V3);
    // MakerNote IFD 里唯一的 tag 位于 ExifIFD 之后
    let pos = tiff
        .windows(2)
        .position(|w| w == 0x0023u16.to_le_bytes())
        .expect("合成文件里应有 tag 0x0023");
    tiff[pos..pos + 2].copy_from_slice(&0x00FFu16.to_le_bytes());
    let e = picture_control::read(&tiff).unwrap_err();
    assert!(matches!(e, nikonrawview::Error::NoPictureControlTag), "得到 {e:?}");
}

#[test]
fn non_tiff_input_reports_no_exif() {
    let junk = vec![0xAAu8; 128];
    let e = picture_control::read(&junk).unwrap_err();
    assert!(matches!(e, nikonrawview::Error::NoExifSegment), "得到 {e:?}");
}

#[test]
fn decodes_adjust_mode_and_scalars() {
    let pc = synth_pc("LINKS-Nature", "NEUTRAL", 0x01, 0x01);
    let id = picture_control::read(&build_tiff(&pc, &MARKER_V3)).unwrap();
    assert_eq!(id.adjust, Adjust::QuickAdjust);
    assert_eq!(id.sharpness, 0.5);
    assert_eq!(id.clarity, 0.5);
    assert_eq!(id.saturation, 5);
    assert_eq!(id.hue, -8);
    assert!(id.curve_sentinel());
}

/// 规格要求：哨兵是**提示**，不得构成"是否需要解析配方库"的短路。
///
/// 实测存在反例——某配方确含 257 级非恒等曲线，其哨兵却为假。因此无论哨兵
/// 取值如何，身份里的配方名都必须被用于解析配方库。
#[test]
fn sentinel_never_short_circuits_library_lookup() {
    let lib = {
        let mut l = nikonrawview::library::Library::new();
        // 放一份与哨兵无关的配方，模拟"配方库里确实存在"
        l.push(nikonrawview::library::Recipe {
            name: "LINKS-Nature".into(),
            source: nikonrawview::library::Source::File { path: "x".into() },
            container: nikonrawview::np3::parse(&{
                let mut p = Vec::new();
                p.extend_from_slice(b"NCP\0");
                p.extend_from_slice(&1u32.to_be_bytes());
                p.extend_from_slice(&4u32.to_be_bytes());
                p.extend_from_slice(b"0310");
                p.extend_from_slice(&0u32.to_be_bytes());
                p
            })
            .unwrap(),
        });
        l
    };

    for (contrast, brightness, expect_sentinel) in
        [(0x01u8, 0x01u8, true), (0x80, 0x80, false), (0x01, 0x80, false)]
    {
        let pc = synth_pc("LINKS-Nature", "NEUTRAL", contrast, brightness);
        let id = picture_control::read(&build_tiff(&pc, &MARKER_V3)).unwrap();
        assert_eq!(id.curve_sentinel(), expect_sentinel);
        // 关键断言：与哨兵取值无关，配方库解析照常进行并命中
        assert_eq!(
            lib.find(&id.name).len(),
            1,
            "哨兵={expect_sentinel} 时仍必须按名解析配方库"
        );
    }
}
