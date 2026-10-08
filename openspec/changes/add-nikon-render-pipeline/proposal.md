## Why

两个痛点中的第二个：**Lightroom Classic 无法呈现相机上载入的第三方 Picture Control 色彩**。

解析层已经完成——能读出"这张照片用了哪个配方"，也能从 NP3/NCP 容器取出配方的完整参数（257 级曲线、8 段色相混合器、3 区校色、标量）。但**还没有任何东西把这些参数变成像素**。

本变更补上这一段，它同时是后续 `add-adl-extraction`（ADL 只能烘焙、需要先有渲染）与 `add-tiff-export`（交付给 LrC 的 TIFF）的共同前置。

## What Changes

- **NEF → 线性 RGB**：引入 `rawler`（纯 Rust）做 RAW 解码，自行实现去马赛克
- **相机线性 RGB → 工作空间线性 RGB**：白平衡、相机色彩矩阵（A / D65 双光源插值）、转入 ProPhoto 线性
- **基准配方渲染与配方施加**：把"基准配方（NEUTRAL 等）+ 配方增量"施加到线性 RGB；含 257 级曲线、8 段混合器、3 区校色与标量
- **基准配方标定**：基准配方的渲染变换在固件里，任何文件都取不到（已实测确认：NX Studio 的 `PicCon.bin` / `PicCon21.bin` 是携带配方记录的模板 TIFF，**不含 257 级曲线**）。因此从 NX Studio 的参考导出拟合，并提供导出清单与拟合工具
- **记录一条运行约束**：本仓库首次引入外部依赖，而沙箱内 `~/.cargo` 不可写、crates.io 拉取被 schannel 拦截，**构建必须提权执行**

## Capabilities

### New Capabilities

- `render/raw-decode`: 把 NEF 解成线性 RGB。覆盖 CFA 解析、去马赛克，以及白平衡系数、相机色彩矩阵、黑白电平、裁切区等标定元数据的取出
- `render/color-pipeline`: 相机线性 RGB → 工作空间线性 RGB。覆盖白平衡、按光源插值的相机矩阵、转入 ProPhoto 线性
- `picture-control/render`: 把基准配方与配方施加到线性 RGB。覆盖曲线、色相混合器、校色、标量的语义
- `picture-control/calibration`: 从参考导出拟合基准配方的变换，并报告拟合质量

### Modified Capabilities

（无——本变更新增能力，不改变既有 `nef/picture-control`、`picture-control/np3`、`picture-control/library` 的规范行为）

## Impact

- **新增依赖**：`rawler`。这是本仓库第一次引入外部 crate——解析层此前刻意保持零依赖。已验证 `rawler` 能解 Z8 的 NEF（277 ms），且白平衡与尺寸结论与 LibRaw 交叉验证一致
- **构建约束**：沙箱内 `~/.cargo` 不可写、crates.io 被 schannel 拦截，`cargo build` 必须提权执行。这会影响后续所有构建步骤
- **需要你的手工投入**：用 NX Studio 重导一批参考样张（同一 NEF、PC=NEUTRAL、ADL 关闭、无其他调整），用于标定基准配方
- **性能**：`rawler` 解 CFA 约 277 ms；去马赛克与后续色彩变换的开销另计。本变更面向**导出**而非浏览——浏览的答案是抽内嵌预览（约 15 ms），属另一个变更
- **许可**：`rawler` 为 MIT/Apache-2.0；去马赛克自行实现，不移植 GPL 代码
