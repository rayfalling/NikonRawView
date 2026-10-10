//! 白平衡来源核对：**我们读到的 `cam_mul` 是不是尼康记录的那组系数？**
//!
//! # 为什么这条现在最该做
//!
//! 诊断链条已经把问题收敛到唯一剩下的一环：
//!
//! | 环节 | 状态 | 依据 |
//! |---|---|---|
//! | 相机矩阵本身 | 没问题 | 跨图残差 1.42×，浮点精度级 |
//! | 参考样本的白平衡 | 没被改过 | NX Studio 显示「原始值」；`.nksc` 里无任何 WB 字段 |
//! | 差异形态 | **逐通道增益，加在矩阵之前** | `M0·diag(w)` 吃掉色度落差的 106% |
//! | 只与白平衡档有关 | 是 | 跨图迁移：同档 ≈ 自用，跨档 12.9× |
//!
//! `diag(w)` 加在矩阵**之前**的位置，正是 LibRaw 里 `scale_colors` / `pre_mul` 所在
//! 的那一步。而**白平衡的读取是 LibRaw 自己实现的，不是从尼康那里抄的**——所以最后
//! 剩下的可能就是它。
//!
//! 本测试把两个来源摆在一起：NEF 的 MakerNote `WB_RBLevels`（tag `0x000C`，尼康自己
//! 写的）与 `libraw_get_cam_mul`（LibRaw 报的）。
//!
//! - 两者一致 → LibRaw 忠实读了，差异在它对 `pre_mul` 的后续处理
//! - 两者不同 → **LibRaw 在读取或换算白平衡时就改了值**，直接命中

mod common;

use common::samples_dir;
use nikonrawview::makernote;

/// MakerNote 里白平衡相关的标签。
const TAG_WB_RBLEVELS: u16 = 0x000C;
const TAG_WB_RBGGLEVELS: u16 = 0x0097;
/// 白平衡版本 / 色彩平衡版本，不同版本的多重定义不同。
const TAG_COLOR_BALANCE_VERSION: u16 = 0x0096;

/// 把标签原始字节按 SHORT 解析成一组数。
fn as_shorts(b: &[u8]) -> Vec<f64> {
    b.chunks_exact(2)
        .map(|c| u16::from_le_bytes([c[0], c[1]]) as f64)
        .collect()
}

/// 把标签原始字节按 RATIONAL（u32 分子 / u32 分母，大端）解析。
fn as_rationals(b: &[u8]) -> Vec<f64> {
    b.chunks_exact(8)
        .filter_map(|c| {
            let n = u32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64;
            let d = u32::from_le_bytes([c[4], c[5], c[6], c[7]]) as f64;
            (d != 0.0).then_some(n / d)
        })
        .collect()
}

/// `WB_RBLevels` 的前两项**本身就是比值** `[R/G, B/G, 1, 1]`（后两项恒为 1）。
///
/// 这一点踩过坑：原先按「R, G, B 三通道电平」解释再除以第二项，得到的
/// `1073/599 = 1.791` 之类的数字与 LibRaw 对不上；实际上 `1073/512 = 2.0957`
/// 才是 R/G，`599/512 = 1.1699` 才是 B/G——与 `cam_mul` 逐位吻合。
fn norm_g(v: &[f64]) -> Option<[f64; 3]> {
    if v.len() < 2 {
        return None;
    }
    Some([v[0], 1.0, v[1]])
}

#[test]
fn cam_mul_matches_nikons_own_wb_levels() {
    let dir = samples_dir();
    let mut any = false;
    eprintln!("=== 白平衡来源核对：MakerNote WB_RBLevels vs libraw cam_mul ===");
    eprintln!();
    eprintln!(
        "  {:<10} {:>22} {:>22} {:>10}",
        "图", "MakerNote 归一化 R/G,B/G", "LibRaw cam_mul R/G,B/G", "判定"
    );

    let mut stems: Vec<String> = std::fs::read_dir(&dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter_map(|e| {
            let p = e.path();
            let stem = p.file_name()?.to_string_lossy().into_owned();
            let stem = stem.strip_suffix(".NEF")?.to_string();
            // 只统计有配对参考导出或样张的
            let has_tif = p.with_extension("TIF").is_file();
            (has_tif || stem == "DSC_4143").then_some(stem)
        })
        .collect();
    stems.sort();
    for stem in stems {
        let stem = stem.as_str();
        let nef = dir.join(format!("{stem}.NEF"));
        if !nef.is_file() {
            continue;
        }
        let Ok(data) = std::fs::read(&nef) else { continue };
        any = true;

        // 尼康自己写的
        let raw = makernote::tag_bytes(&data, TAG_WB_RBLEVELS).ok().flatten();
        let mk_short = raw.as_deref().map(as_shorts).and_then(|v| norm_g(&v));
        let mk_rat = raw.as_deref().map(as_rationals).and_then(|v| norm_g(&v));
        // 优先 RATIONAL：WB_RBLevels 是 4 个 RATIONAL，按 SHORT 读会得到一串电平而非比值
        let mk = mk_rat.or(mk_short);

        // LibRaw 报的
        let lr = nikonrawview::libraw::read_wb(&nef).ok().map(|w| {
            let g = w[1] as f64;
            [w[0] as f64 / g, 1.0, w[2] as f64 / g]
        });

        let raw_len = raw.as_ref().map(|b| b.len()).unwrap_or(0);
        let (mks, lrs, verdict) = match (&mk, &lr) {
            (Some(a), Some(b)) => {
                let d = (a[0] - b[0]).abs() + (a[2] - b[2]).abs();
                let v = if d < 0.02 {
                    "一致".to_string()
                } else {
                    format!("**不一致（差 {d:.4}）**")
                };
                (
                    format!("{:.5} {:.5}", a[0], a[2]),
                    format!("{:.5} {:.5}", b[0], b[2]),
                    v,
                )
            }
            (None, Some(b)) => (
                format!("（无标签，{raw_len} 字节）"),
                format!("{:.5} {:.5}", b[0], b[2]),
                "无法比对".to_string(),
            ),
            (Some(a), None) => (
                format!("{:.5} {:.5}", a[0], a[2]),
                "（读取失败）".to_string(),
                "无法比对".to_string(),
            ),
            (None, None) => ("—".into(), "—".into(), "两者都取不到".into()),
        };
        eprintln!("  {:<10} {:>22} {:>22} {:>10}", stem, mks, lrs, verdict);
        if let Some(r) = raw.as_ref() {
            eprintln!("       WB_RBLevels 原始 {} 字节：{:02X?}", r.len(), &r[..r.len().min(24)]);
        }

        // 黑电平：LibRaw 的唯一来源就是 MakerNote tag 0x003d，用同一个解析器读它。
        // 诊断 L 已从源码确证 LibRaw 读的就是这里，且逐通道扣、只扣一次、整数、钳到 0。
        match makernote::black_level(&data) {
            Ok(Some(bl)) => eprintln!("       黑电平：{}", bl.describe()),
            Ok(None) => eprintln!("       黑电平：**读不到 tag 0x003d**（不要当成 0）"),
            Err(e) => eprintln!("       黑电平：解析失败 {e}"),
        }

        // 顺带看多重定义的变体与版本
        if let Ok(Some(b)) = makernote::tag_bytes(&data, TAG_WB_RBGGLEVELS) {
            let v = as_rationals(&b);
            eprintln!(
                "       WB_RBGGLevels（{} 项）：{}",
                v.len(),
                v.iter().take(5).map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(" ")
            );
        }
        if let Ok(Some(b)) = makernote::tag_bytes(&data, TAG_COLOR_BALANCE_VERSION) {
            eprintln!("       ColorBalanceVersion 原始：{:02X?}", &b[..b.len().min(8)]);
        }
        eprintln!();
    }

    if !any {
        eprintln!("跳过：没有可用 NEF");
    }

    eprintln!("判读：");
    eprintln!("  两者一致 → LibRaw 忠实读了尼康的系数，差异在它对 pre_mul 的后续处理；");
    eprintln!("  两者不同 → LibRaw 在读取或换算白平衡时就改了值，这就是命中点。");
}
