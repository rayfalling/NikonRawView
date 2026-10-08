//! `nikon-pc` —— 打印一张照片所用的 Picture Control 配方。
//!
//! ```text
//! nikon-pc <文件>...            打印配方身份；若能解析出配方库则一并打印本体
//! nikon-pc --json <文件>...     机器可读输出
//! nikon-pc --library            列出本机可用配方
//! nikon-pc --sidecars <目录>    追加边车目录（可重复）
//! ```

use std::path::PathBuf;
use std::process::ExitCode;

use nikonrawview::library::{self, Library, Recipe};
use nikonrawview::np3::{self, Scale};
use nikonrawview::picture_control::{self, Adjust, Identity};

struct Opts {
    json: bool,
    library_only: bool,
    sidecars: Vec<PathBuf>,
    files: Vec<PathBuf>,
}

fn usage() -> &'static str {
    "用法：nikon-pc [--json] [--library] [--sidecars <目录>]... <文件>...\n\
     退出码：0 全部成功；1 有文件解析失败；2 参数错误"
}

fn parse_args() -> Result<Opts, String> {
    let mut o = Opts { json: false, library_only: false, sidecars: Vec::new(), files: Vec::new() };
    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--json" => o.json = true,
            "--library" => o.library_only = true,
            "--sidecars" => {
                let v = it.next().ok_or("--sidecars 需要一个目录参数")?;
                o.sidecars.push(PathBuf::from(v));
            }
            "-h" | "--help" => return Err(String::new()),
            s if s.starts_with("--") => return Err(format!("未知选项：{s}")),
            _ => o.files.push(PathBuf::from(a)),
        }
    }
    if o.files.is_empty() && !o.library_only {
        return Err("至少需要一个文件，或使用 --library".into());
    }
    Ok(o)
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn adjust_label(a: Adjust) -> String {
    match a {
        Adjust::DefaultSettings => "Default Settings".into(),
        Adjust::QuickAdjust => "Quick Adjust".into(),
        Adjust::FullControl => "Full Control".into(),
        Adjust::Other(v) => format!("其他({v})"),
    }
}

fn print_identity_text(path: &std::path::Path, id: &Identity, hits: &[&Recipe]) {
    println!("{}", path.display());
    println!(
        "  配方      : {}   （基准：{}{}）",
        id.name,
        id.base,
        if id.is_custom() { "，自定义" } else { "，内置" }
    );
    println!(
        "  载荷      : {} 字节，版本 {}，调整模式 {}",
        id.raw.len(),
        id.version,
        adjust_label(id.adjust)
    );
    // 哨兵值 0x01 的语义是"该参数不适用"，不是 −127；直接打印数值会误导。
    let (contrast, brightness) = if id.curve_sentinel() {
        ("不适用(哨兵)".to_string(), "不适用(哨兵)".to_string())
    } else {
        (format!("{:+}", id.contrast), format!("{:+}", id.brightness))
    };
    println!(
        "  标量      : 锐化 {:+.2}  清晰度 {:+.2}  对比度 {}  亮度 {}  饱和度 {:+}  色相 {:+}",
        id.sharpness, id.clarity, contrast, brightness, id.saturation, id.hue
    );
    println!(
        "  曲线哨兵  : {}（提示，不构成是否需要解析配方库的判据）",
        if id.curve_sentinel() { "存在" } else { "无" }
    );
    if hits.is_empty() {
        println!("  配方本体  : 未在本机配方库中找到（用 --sidecars 指定边车目录，或核对名称）");
    } else {
        for r in hits {
            let c = &r.container;
            println!(
                "  配方本体  : {} <- {}  代次={}  含曲线={}",
                r.name,
                r.source.label(),
                r.family(),
                c.has_custom_curve()
            );
            if let Some(curve) = &c.curve {
                println!(
                    "              曲线：{} 控制点，257 级 LUT（偏移 {}，单调={}）",
                    curve.point_count,
                    curve.lut_offset,
                    curve.is_monotonic()
                );
            }
            let fmt = |v: Option<f32>| v.map(|x| format!("{x:+.2}")).unwrap_or_else(|| "-".into());
            println!(
                "              标量：锐化 {}  清晰度 {}  中频锐化 {}  对比度 {}  饱和度 {}",
                fmt(c.sharpening()),
                fmt(c.clarity()),
                fmt(c.mid_range_sharpening()),
                fmt(c.contrast()),
                fmt(c.saturation())
            );
            if let Some(b) = &c.blender {
                let segs: Vec<String> = np3::BAND_ORDER
                    .iter()
                    .zip(b.iter())
                    .map(|(n, s)| {
                        format!("{} {:+}/{:+}/{:+}", &n[..2], s.hue, s.chroma, s.brightness)
                    })
                    .collect();
                println!("              混合器：{}", segs.join("  "));
            }
        }
    }
}

fn print_identity_json(path: &std::path::Path, id: &Identity, hits: &[&Recipe]) {
    let mut s = String::new();
    s.push('{');
    s.push_str(&format!("\"file\":\"{}\"", json_escape(&path.to_string_lossy())));
    s.push_str(&format!(",\"name\":\"{}\"", json_escape(&id.name)));
    s.push_str(&format!(",\"base\":\"{}\"", json_escape(&id.base)));
    s.push_str(&format!(",\"custom\":{}", id.is_custom()));
    s.push_str(&format!(",\"payload_len\":{}", id.raw.len()));
    s.push_str(&format!(",\"version\":\"{}\"", json_escape(&id.version)));
    s.push_str(&format!(",\"adjust\":\"{}\"", adjust_label(id.adjust)));
    s.push_str(&format!(",\"curve_sentinel\":{}", id.curve_sentinel()));
    s.push_str(&format!(",\"sharpness\":{}", id.sharpness));
    s.push_str(&format!(",\"clarity\":{}", id.clarity));
    s.push_str(&format!(",\"contrast\":{}", id.contrast));
    s.push_str(&format!(",\"brightness\":{}", id.brightness));
    s.push_str(&format!(",\"saturation\":{}", id.saturation));
    s.push_str(&format!(",\"hue\":{}", id.hue));
    s.push_str(",\"library\":[");
    for (i, r) in hits.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!(
            "{{\"source\":\"{}\",\"family\":\"{}\",\"has_custom_curve\":{}}}",
            json_escape(&r.source.label()),
            json_escape(r.family()),
            r.container.has_custom_curve()
        ));
    }
    s.push_str("]}");
    println!("{s}");
}

fn print_library(lib: &Library, json: bool) {
    if json {
        let items: Vec<String> = lib
            .recipes
            .iter()
            .map(|r| {
                format!(
                    "{{\"name\":\"{}\",\"source\":\"{}\",\"family\":\"{}\"}}",
                    json_escape(&r.name),
                    json_escape(&r.source.label()),
                    json_escape(r.family())
                )
            })
            .collect();
        println!("{{\"recipes\":[{}]}}", items.join(","));
        return;
    }
    if lib.recipes.is_empty() {
        println!("本机未发现任何配方来源（注册表不可读，且未给出边车目录）");
        return;
    }
    println!("本机配方（去重后 {} 个，共 {} 条来源）：", lib.names().len(), lib.recipes.len());
    for n in lib.names() {
        println!("  {n}");
        for r in lib.find(&n) {
            println!(
                "      <- {:<52} 代次={} 含曲线={} 混合器={}",
                r.source.label(),
                r.family(),
                r.container.has_custom_curve(),
                r.container.blender.is_some()
            );
        }
    }
}

fn main() -> ExitCode {
    let opts = match parse_args() {
        Ok(o) => o,
        Err(e) => {
            if !e.is_empty() {
                eprintln!("{e}");
            }
            eprintln!("{}", usage());
            return if e.is_empty() { ExitCode::SUCCESS } else { ExitCode::from(2) };
        }
    };

    let lib = library::scan_all(&opts.sidecars);

    if opts.library_only {
        print_library(&lib, opts.json);
        return ExitCode::SUCCESS;
    }

    let mut failed = 0usize;
    for path in &opts.files {
        let data = match std::fs::read(path) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("{}：无法读取（{e}）", path.display());
                failed += 1;
                continue;
            }
        };
        match picture_control::read(&data) {
            Ok(id) => {
                let hits = lib.find(&id.name);
                if opts.json {
                    print_identity_json(path, &id, &hits);
                } else {
                    print_identity_text(path, &id, &hits);
                }
            }
            Err(e) => {
                eprintln!("{}：解析失败（{e}）", path.display());
                failed += 1;
            }
        }
    }

    // 顺带暴露一个便于调试的标量解码入口，避免未使用告警
    let _ = Scale::Quarter;

    if failed == 0 {
        ExitCode::SUCCESS
    } else {
        eprintln!("{failed} 个文件解析失败");
        ExitCode::FAILURE
    }
}
