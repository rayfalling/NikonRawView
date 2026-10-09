## Why

两个痛点中的第二个：**Lightroom Classic 无法呈现相机上载入的第三方 Picture Control 色彩**。

解析层已经完成——能读出"这张照片用了哪个配方"，也能从 NP3/NCP 容器取出配方的完整参数（257 级曲线、8 段色相混合器、3 区校色、标量）。但**还没有任何东西把这些参数变成像素**。

本变更补上这一段，它同时是后续 `add-adl-extraction`（ADL 只能烘焙、需要先有渲染）与 `add-tiff-export`（交付给 LrC 的 TIFF）的共同前置。

## What Changes

- **NEF → 线性相机 RGB**：用 LibRaw 解码（去马赛克亦由 LibRaw 提供），输出线性相机空间像素并取出标定元数据
- **相机线性 RGB → 工作空间线性 RGB**：白平衡、相机 → 输出空间矩阵、由本管线自行施加并转入 ProPhoto 线性
- **基准配方渲染与配方施加**：把"基准配方（NEUTRAL 等）+ 配方增量"施加到线性 RGB；含 257 级曲线、8 段混合器、3 区校色与标量
- **基准配方标定**：基准配方的渲染变换在固件里，任何文件都取不到（已实测确认：NX Studio 的 `PicCon.bin` / `PicCon21.bin` 是携带配方记录的模板 TIFF，**不含 257 级曲线**）。因此从 NX Studio 的参考导出拟合，并提供导出清单与拟合工具
- **许可合规**：LibRaw 走 **CDDL-1.0**，源码随仓库分发，并补齐仓库自身缺失的 LICENSE 文件与第三方许可清单
- **与尼康工坊交叉验证**：渲染结果与 NX Studio 对同一 NEF 的导出比对，ΔE00 数值记入仓库作为回归基线

## Capabilities

### New Capabilities

- `render/raw-decode`: 把 NEF 解成线性相机 RGB。覆盖标定元数据的取出、线性输出的参数设置、去马赛克算法选择与可复现性
- `render/color-pipeline`: 相机线性 RGB → 工作空间线性 RGB。覆盖白平衡、相机 → 输出空间矩阵的取得与自行施加、转入 ProPhoto 线性
- `picture-control/render`: 把基准配方与配方施加到线性 RGB。覆盖曲线、色相混合器、校色、标量的语义
- `picture-control/calibration`: 从参考导出拟合基准配方的变换，并报告拟合质量

### Modified Capabilities

（无——本变更新增能力，不改变既有 `nef/picture-control`、`picture-control/np3`、`picture-control/library` 的规范行为）

## Impact

- **新增依赖**：LibRaw 0.22（C 库，源码 vendor 在 `third_party/libraw/`）。**crate 依赖仍为零**——以 `build.rs` + 手写 FFI 接入。曾评估纯 Rust 的 `rawler`，因其为 LGPL-2.1 而改用 CDDL 的 LibRaw，顺带把依赖从 129 个 crate 降到 1 个 C 库
- **构建方式**：`nmake /f Makefile.msvc` + `vcvars64`。**不能用 MSBuild**——其 vcxproj 钉死 Windows SDK `10.0.18362.0` 与工具集 `v142`，且在沙箱内会被「同名不同大小写的代理环境变量」触发的 .NET 字典异常打断
- **需要你的手工投入**：用 NX Studio 重导一批参考样张（同一 NEF、PC=NEUTRAL、ADL 关闭、无其他调整），用于标定基准配方
- **性能**：LibRaw 解一张 45 MP 的 NEF 并完成去马赛克约 3.7 s。本变更面向**导出**而非浏览——浏览的答案是抽内嵌预览（约 15 ms），属另一个变更
- **许可**：本仓库自身为 MIT OR Apache-2.0；LibRaw 为 **CDDL-1.0**（文件级 copyleft）。CDDL 允许链入任意许可的整体作品，不要求可重新链接，但分发二进制时须随附 LibRaw 的源码——`third_party/libraw/` 即满足
