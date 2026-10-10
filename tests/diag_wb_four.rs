//! 四元白平衡：**`G1` 与 `G2` 是否相等，是尼康调色偏移（A/M 微调）的直接编码处**。
//!
//! # 为什么看这个
//!
//! 用户指出色温扫描组带 **A2.0 / M2.25** 的调色偏移（从 `DSC_0563` 起全程固定）。
//! 而 [`nikonrawview::libraw::read_wb`] 只返回四元组的前三项，把 `G2` 丢了。
//!
//! 尼康的白平衡微调有两条轴：
//!
//! ```text
//! A-B 轴（琥珀-蓝）：改变 R 与 B 的相对关系
//! G-M 轴（绿-品红）：改变 G 与 (R,B) 的相对关系  ← 这一维只体现在 G 上
//! ```
//!
//! 所以 `G1 != G2` 时，那个差值就是调色偏移的读数；`G1 == G2` 时没有可读的偏移。
//! 本测试把四元组打出来，看同一组内是否恒定、跨组是否不同。

use std::path::Path;

fn show(dir: &Path, stem: &str) {
    let p = dir.join(format!("{stem}.NEF"));
    if !p.is_file() {
        eprintln!("  {stem:<12} （缺文件）");
        return;
    }
    let k = nikonrawview::makernote::color_temperature(&std::fs::read(&p).unwrap())
        .ok()
        .flatten();
    match nikonrawview::libraw::read_wb4(&p) {
        Ok(w) => {
            let g = if w[3].abs() > 1e-9 { (w[1] / w[3]) as f64 } else { f64::NAN };
            let bg = (w[2] / w[1]) as f64;
            let rg = (w[0] / w[1]) as f64;
            eprintln!(
                "  {stem:<12} K={:<6} R/G={rg:.5} B/G={bg:.5}  G1={:.5} G2={:.5}  G1/G2={g:.6}{}",
                k.map(|v| v.to_string()).unwrap_or_else(|| "?".into()),
                w[1],
                w[3],
                if (w[1] - w[3]).abs() < 1e-7 { "（相等）" } else { "  <-- 不等" }
            );
        }
        Err(e) => eprintln!("  {stem:<12} 读取失败：{e}"),
    }
}

#[test]
fn four_element_white_balance_shows_the_tint() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("simple");

    eprintln!("=== 色温扫描组（用户说 A2.0 / M2.25 全程固定，从 563 起）===");
    for i in 561..=570 {
        show(&dir, &format!("DSC_{i:04}"));
    }

    eprintln!("\n=== 其它参考（调色偏移应当不同）===");
    for s in ["DSC_0001", "DSC_0010", "DSC_0141", "DSC_0314", "DSC_0345", "DSC_0507", "DSC_0540", "DSC_8562"] {
        show(&dir, s);
    }

    eprintln!("\n判读：");
    eprintln!("  G1/G2 在扫描组内应恒定（用户说偏移全程固定），而与其它组不同；");
    eprintln!("  若 G1/G2 恒为 1，则尼康把该偏移存在别处（可能是加密块），");
    eprintln!("  此时只能用 R/G 作为调色偏移的代理变量。");
}
