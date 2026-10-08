## Why

解决两个痛点：**尼康工坊（NX Studio）浏览 NEF 太慢**；**Lightroom Classic 无法还原相机上载入的第三方 Picture Control 色彩**（本项目的配方均为第三方下载的 NCP/NP3 文件）。

本次变更只做**解析层**。它是后续所有能力的前置：不先知道"这张照片用了哪个配方"，就无从渲染、选图与导出。

实测结论决定了这一层的形态：**NEF 装不下完整配方**。相机只写入配方身份（约 108 字节）+ 标量滑块 + 一个曲线标志位，而完整配方（含 257 级色调曲线）在 NP3/NCP 文件里。因此解析必须做成两段式——**从文件读身份，按名字去配方库取本体**。

## What Changes

- 新增 **NEF/JPEG 配方身份读取**：走 Nikon MakerNote tag `0x0023`（`PictureControl3` 布局），取出自定义配方名、基准色彩名、调整模式、标量值、曲线哨兵
- 新增 **NP3/NCP 配方容器解析**，**同时支持两代模板**（`0310` 与 `0300`）
- 新增 **本机配方库发现**：注册表 `HKCU\Software\Nikon\Common\PictureControl\NP3\CustomCurves` 与 NX Studio `.nksc` 边车
- 新增 **CLI**：给定一个 NEF/JPEG，打印所用配方身份、解析出的配方本体参数、以及解析来源
- 将两条**实测踩坑**固化为规范约束：禁止以字节模式匹配定位 Nikon 元数据；禁止依赖前导头中那个可变 u16

## Capabilities

### New Capabilities

- `nef/picture-control`: 从 NEF/JPEG 的 MakerNote 读出配方身份。覆盖 TIFF/IFD 遍历、NEF 与 JPEG 两种基址规则、`PictureControl3` 字段布局、曲线哨兵语义
- `picture-control/np3`: 解析 NP3/NCP 配方容器。覆盖 chunk 遍历（禁用固定偏移）、两代模板、曲线块、8 段色相混合器、3 区校色、base 编码
- `picture-control/library`: 发现本机可用配方并建立索引，按配方名解析出配方本体，并报告来源与代次

### Modified Capabilities

（无——本仓库尚无既有 specs）

## Impact

- **范围**：纯解析。不含渲染、GPU、色彩拟合、LrC 交互
- **语言**：Rust（后续渲染与选图同样使用 Rust）
- **外部依赖**：无。TIFF/IFD 与 chunk 解析均自行实现
- **许可**：不引入受限第三方代码。ExifTool（GPL/Artistic 双许可）仅作为字段定义的参考资料，不拷贝其源码
- **测试数据**：`simple/DSC_4143.NEF` + `.JPG`（配对，108 字节载荷逐字节一致）作为验收向量
- **风险**：NP3 容器的 `0300` 代次缺少公开布局文档，其 chunk 语义需靠本机 3 个真实样本反推
