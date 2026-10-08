//! 针对本机真实样张的测试。
//!
//! 样张不入库（见 `.gitignore`），因此这些用例在样本缺席时打印提示并跳过；
//! 期望值本身作为夹具提交在 `tests/common/mod.rs` 中。

mod common;

use common::*;
use nikonrawview::library;
use nikonrawview::{np3, picture_control};

#[test]
fn dsc4143_nef_identity() {
    let Some(data) = sample(DSC4143.file) else {
        return skip("dsc4143_nef_identity", DSC4143.file);
    };
    let id = picture_control::read(&data).expect("应能解析");
    assert_eq!(id.name, DSC4143.name);
    assert_eq!(id.base, DSC4143.base);
    assert_eq!(id.raw.len(), DSC4143.payload_len);
    assert_eq!(id.curve_sentinel(), DSC4143.curve_sentinel);
    assert_eq!(id.sharpness, DSC4143.sharpness);
    assert_eq!(id.clarity, DSC4143.clarity);
    assert_eq!(id.saturation, DSC4143.saturation);
    assert_eq!(id.hue, DSC4143.hue);
}

/// 配对 NEF 与 JPEG 的载荷必须逐字节相同——可用于交叉校验两侧解析器。
#[test]
fn dsc4143_nef_and_jpeg_payloads_match() {
    let (Some(nef), Some(jpg)) = (sample(DSC4143.file), sample("DSC_4143.JPG")) else {
        return skip("dsc4143_nef_and_jpeg_payloads_match", "DSC_4143.NEF/.JPG");
    };
    let a = picture_control::read(&nef).expect("NEF 应能解析");
    let b = picture_control::read(&jpg).expect("JPEG 应能解析");
    assert_eq!(a, b, "两侧载荷应逐字节相同");
}

/// 本机配方库：注册表 + 边车，合并后应能按名解析出配方本体。
#[test]
fn local_library_resolves_the_sample_recipe() {
    let Some(nef) = sample(DSC4143.file) else {
        return skip("local_library_resolves_the_sample_recipe", DSC4143.file);
    };
    let id = picture_control::read(&nef).unwrap();

    let roots: Vec<std::path::PathBuf> = std::env::var_os("NIKONRAWVIEW_SIDECARS")
        .map(|v| v.to_string_lossy().split(';').map(std::path::PathBuf::from).collect())
        .unwrap_or_default();
    let lib = library::scan_all(&roots);
    if lib.recipes.is_empty() {
        eprintln!("跳过：本机无可读配方来源（注册表 / 边车）");
        return;
    }
    eprintln!("本机配方：{:?}", lib.names());
    let hit = lib.find(&id.name);
    assert!(
        !hit.is_empty(),
        "样张用了 {}，但本机配方库里没有同名配方：{:?}",
        id.name,
        lib.names()
    );
    for r in hit {
        eprintln!(
            "  {} <- {}  代次={}  含曲线={}",
            r.name,
            r.source.label(),
            r.family(),
            r.container.has_custom_curve()
        );
    }
}

/// 读文件前缀即可拿到配方身份；前缀不足时退回整份读取。
///
/// MakerNote 在第 35 KB 附近、载荷紧随其后，因此 512 KB 前缀对 NEF 与 JPEG 都
/// 绰绰有余。整份读入一个 50 MB 的 NEF 只为取 108 字节，会让全库校验产生
/// 300 GB 以上的无谓 I/O。
fn parse_file(p: &std::path::Path) -> Result<nikonrawview::picture_control::Identity, String> {
    use std::io::Read;
    const PREFIX: u64 = 512 * 1024;

    let mut f = std::fs::File::open(p).map_err(|e| e.to_string())?;
    let mut buf = Vec::new();
    f.by_ref()
        .take(PREFIX)
        .read_to_end(&mut buf)
        .map_err(|e| e.to_string())?;

    match picture_control::read(&buf) {
        Ok(id) => Ok(id),
        Err(first) => {
            let full = std::fs::read(p).map_err(|e| e.to_string())?;
            picture_control::read(&full).map_err(|_| first.to_string())
        }
    }
}

/// 全库批量校验。默认忽略，因为要读上万个文件：
///
/// ```text
/// $env:NIKONRAWVIEW_LIBRARY="E:\Nikon\Z8"
/// cargo test --test local_samples -- --ignored --nocapture
/// ```
#[test]
#[ignore = "需指定 NIKONRAWVIEW_LIBRARY，且会读取上万个文件"]
fn batch_parses_whole_library() {
    let root = std::env::var("NIKONRAWVIEW_LIBRARY")
        .expect("请设置 NIKONRAWVIEW_LIBRARY 指向包含 NEF/JPG 的目录");
    let mut total = 0usize;
    let mut ok = 0usize;
    let mut with_identity = std::collections::BTreeMap::<String, usize>::new();
    let mut failures: Vec<String> = Vec::new();
    // stem -> (扩展名, 载荷)：用于配对复验
    let mut by_stem = std::collections::BTreeMap::<String, Vec<(String, Vec<u8>)>>::new();

    let mut stack = vec![std::path::PathBuf::from(&root)];
    while let Some(dir) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&dir) else { continue };
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
                continue;
            }
            let ext = p.extension().map(|x| x.to_string_lossy().to_uppercase()).unwrap_or_default();
            if ext != "NEF" && ext != "JPG" && ext != "JPEG" {
                continue;
            }
            total += 1;
            match parse_file(&p) {
                Ok(id) => {
                    ok += 1;
                    *with_identity
                        .entry(format!("{} / {}", id.name, id.base))
                        .or_default() += 1;
                    if let Some(stem) = p.file_stem().map(|s| s.to_string_lossy().to_string()) {
                        by_stem.entry(stem).or_default().push((ext.clone(), id.raw.clone()));
                    }
                }
                Err(e) => {
                    if failures.len() < 20 {
                        failures.push(format!("{}: {e}", p.display()));
                    }
                }
            }
        }
    }

    eprintln!("共 {total} 个文件，解析成功 {ok}");
    for (k, v) in &with_identity {
        eprintln!("  {v:6}  {k}");
    }
    for f in &failures {
        eprintln!("  失败 {f}");
    }

    // 配对复验：同名 NEF/JPG 的载荷必须逐字节相同（跨全部行程，而非单个样本）
    let mut pairs = 0usize;
    let mut mismatched = Vec::new();
    for (stem, entries) in &by_stem {
        if entries.len() < 2 {
            continue;
        }
        let first = &entries[0].1;
        for (ext, raw) in &entries[1..] {
            pairs += 1;
            if raw != first && mismatched.len() < 10 {
                mismatched.push(format!("{stem}：{} 与 {} 载荷不一致", entries[0].0, ext));
            }
        }
    }
    eprintln!("配对复验：{pairs} 对同名词条，载荷不一致 {} 对", mismatched.len());
    for m in &mismatched {
        eprintln!("  {m}");
    }

    assert_eq!(ok, total, "存在解析失败的文件");
    assert!(total > 0, "未找到任何 NEF/JPG");
    assert!(pairs > 0, "未发现任何可配对的同名文件");
    assert!(mismatched.is_empty(), "配对载荷不一致");
}

/// 注册表四份配方的名称与 base 编码必须与实测一致，且 `0310` 代含曲线者
/// 应解析出 5 个控制点与单调 LUT。
#[test]
fn registry_profiles_match_measured_values() {
    let lib = library::Library {
        recipes: library::scan_registry(),
    };
    if lib.recipes.is_empty() {
        eprintln!("跳过：注册表不可读（非 Windows 或键不存在）");
        return;
    }
    // 实测：名称 -> base 编码
    let expected: &[(&str, u16)] = &[
        ("Kodak Ektar Green", 0x00C3),
        ("Kodak-Sun-Nature", 0x03C2),
        ("LINKS-Nature", 0x03C2),
        ("MSLT-XTRA400-V0.5", 0x0020),
    ];
    eprintln!("注册表 {} 份：{:?}", lib.recipes.len(), lib.names());
    for (name, code) in expected {
        let found = lib.find(name);
        assert!(!found.is_empty(), "注册表应含 {name}，实际 {:?}", lib.names());
        let r = found[0];
        assert_eq!(
            r.container.base_code,
            Some(*code),
            "{name} 的 base 编码应为 0x{code:04X}"
        );
        eprintln!(
            "  {name:<22} base=0x{:04X}  代次={}  含曲线={}",
            r.container.base_code.unwrap_or(0),
            r.family(),
            r.container.has_custom_curve()
        );
    }

    // 曲线：`0310` 代的 MSLT 应给出 5 个控制点与单调 LUT
    let mslt = lib.find("MSLT-XTRA400-V0.5")[0];
    assert_eq!(mslt.family(), "0310");
    let curve = mslt.container.curve.as_ref().expect("MSLT 应含曲线块");
    assert_eq!(curve.point_count, 5);
    assert_eq!(curve.points.len(), 5);
    assert_eq!(curve.lut.len(), 257);
    assert!(curve.is_monotonic());
    assert!(mslt.container.has_custom_curve(), "MSLT 的曲线应非恒等");
    eprintln!("  MSLT 曲线：LUT 偏移 {}，首={}，末={}",
        curve.lut_offset, curve.lut[0], curve.lut[256]);
}

/// 边车来源：载荷应为 978 字节的 `0310` 代，含 5 个控制点。
#[test]
fn sidecar_payload_is_a_valid_container() {
    let roots: Vec<std::path::PathBuf> = std::env::var_os("NIKONRAWVIEW_SIDECARS")
        .map(|v| v.to_string_lossy().split(';').map(std::path::PathBuf::from).collect())
        .unwrap_or_default();
    if roots.is_empty() {
        eprintln!("跳过：未设置 NIKONRAWVIEW_SIDECARS");
        return;
    }
    let mut n = 0;
    for r in &roots {
        for rec in library::scan_sidecars(r) {
            n += 1;
            assert!(!rec.name.is_empty(), "边车配方应有名称");
            assert!(
                rec.container.family == "0310" || rec.container.family == "0300",
                "家族码应为 0310 或 0300，得到 {}",
                rec.container.family
            );
            if let Some(c) = &rec.container.curve {
                assert_eq!(c.lut.len(), 257);
                assert!(c.is_monotonic(), "曲线 LUT 应单调");
                let _ = np3::parse;
            }
        }
    }
    eprintln!("边车解析出 {n} 份配方");
}
