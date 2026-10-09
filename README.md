# NikonRawView

解决两件事：**尼康工坊（NX Studio）看图太慢**，以及 **Lightroom Classic 无法还原相机上载入的自定义 Picture Control 色彩**。

核心思路：**从 NEF 读出「用了哪个配方」，从本机配方库取出配方本体，用它渲染——同时用 NEF 内嵌的相机预览实现秒开浏览。**

> 项目刚初始化，目前只有这份说明与 `.gitignore`，**尚未写任何代码**。下文描述的是经实测确定的设计意图。

---

## 问题

### 1. 尼康工坊看图慢

尼康工坊（NX Studio）浏览、筛选大量 NEF 时响应很慢，逐张打开确认对焦与构图的开销尤其明显，不适合做快速初筛。

根因是它每次都在**解 RAW**：本机实测一张 45 MP NEF 的完整去马赛克+后处理需要 **约 1000 ms**。

### 2. Lightroom Classic 无法还原自定义 Picture Control

相机上的色彩（本项目中指**第三方下载并载入相机的 NCP/NP3 配方**，例如 LINKS-Nature、MSLT 系列）在拍摄时由相机 ISP **烘焙进机内 JPEG**。而：

- 机内 JPEG 里看到的是**已应用该配方**的结果；
- Lightroom Classic 只会用自己那套相机配置文件解释 RAW。Adobe 的「Camera Matching」只覆盖**内置**色彩（Standard / Neutral / Vivid / Portrait / Landscape / Flat / Monochrome 等），**不包含任何自定义或第三方配方**；
- 结果就是**自定义配方的色彩在 LrC 里丢失或明显偏色**，NEF 与机内 JPEG 观感对不上。

## 一个决定架构的实测结论

**NEF 里装不下完整配方。** 相机只把「配方身份 + 滑块参数 + 有无自定义曲线的标记」写进 MakerNote：

| 来源 | 含有什么 | 完整配方？ |
|---|---|---|
| NEF / JPEG 的 MakerNote | 配方名、基准色彩名、标量滑块、曲线哨兵 | **否** |
| NCP / NP3 配方文件 | 完整参数 + 257 级 15-bit 色调曲线 + 8 段色相混合器 + 3 区校色 | **是** |
| 相机固件 | 基准色彩自身的曲线与混合器权重 | 取不到 |

所以任何「只靠解析 NEF 就能还原色彩」的方案都不成立，架构必须是混合式。

## 思路

```
NEF / JPEG ──→ 身份：配方名 + 基准色彩 + 是否曲线型
                    │
                    ▼  按名字查本机配方库
NCP / NP3  ──→ 配方本体：完整参数 + 色调曲线 LUT + 混合器 + 校色
                    │
        ┌───────────┴───────────┐
        ▼                       ▼
   快速看图                  色彩还原
 （内嵌预览直出）        （DCP 配置文件 / GPU 渲染）
```

配方库的来源：相机随附的配方文件、NX Studio 的导入记录、以及 NX Studio 保存的编辑边车文件。

### 快速看图为什么不需要 GPU

NEF 内嵌了多条 JPEG 流，其中带**屏幕尺寸（约 1620×1080）**和**全分辨率**两档，都在文件结构中正规寻址。这些预览**已经烘焙了相机当前的配方**。

| 路径 | 实测耗时 |
|---|---|
| 抽取内嵌预览 | **约 15 ms** |
| 完整 45 MP 去马赛克 | 约 1000 ms |

即使 GPU 把去马赛克压到 50 ms，仍比直接抽预览慢数倍，而色彩还不如预览准确。**GPU 只在编辑路径（改曝光、高光恢复）才有意义**，那是另一个目标。

### LrC 的交付物是 DCP，不是 LUT

Lightroom Classic **不支持 `.cube` 3D LUT**。它能装的是自定义 **DCP 相机配置文件**，直接作用于 NEF，无需转 DNG。

因此交付物为：**沿用 Adobe 该机型的色彩矩阵，替换其中的观感表**——即在 Adobe 已经正确的基准上叠加一个真正的 delta，并可选附一份 XMP 预设供一键套用。

## 计划中的流程

1. **配方解析** —— 解析 NCP/NP3 容器，取出完整参数与色调曲线；从 NEF 读出配方身份。
2. **配方库** —— 汇集本机可用的配方，按名字与基准色彩索引。
3. **快速看图** —— 抽取内嵌预览，实现秒开浏览与初筛。
4. **标定** —— 以成对的 RAW+JPEG 为样本，在干净目标（关闭相机其它处理）上求解 delta。
5. **评估** —— 用 ΔE 量化误差，并做可视化对比。
6. **导出** —— 生成 LrC 可用的 DCP 相机配置文件（可选附预设）。

## 规划中的目录（暂未创建）

```
docs/      设计说明与实测结论
src/       配方解析、预览抽取、标定
samples/   配对样张（不入库）
out/       生成的 DCP / 预设
```

## 现状

- [x] 仓库初始化（README + `.gitignore`）
- [x] **技术路径验证**：NEF 配方取证、NCP/NP3 容器结构、RAW 解码库可用性与性能实测
- [ ] 配方解析器（NCP/NP3 读取 + NEF 配方身份读取）
- [ ] 快速看图器
- [ ] 标定与 DCP 导出

## 构建

```bash
cargo build --release
```

### 前置条件

RAW 解码使用 [LibRaw](https://www.libraw.org/)，其源码随仓库分发在 `third_party/libraw/`，由 `build.rs` 自动构建。因此构建机需要：

- **Visual Studio 的「使用 C++ 的桌面开发」工作负载**（提供 `cl.exe`、`nmake` 与 `vcvars64.bat`）
- Rust 工具链

`build.rs` 通过 `vswhere` 定位 `vcvars64.bat`；若定位失败，可设环境变量 `VCVARS64` 指向它的完整路径。

首次构建会花约 2 分钟编译 LibRaw（80 个源文件）；静态库已存在时跳过。

### 为什么用 nmake 而不是 MSBuild

LibRaw 0.22 **没有 CMakeLists.txt**，官方构建路径之一是 `Makefile.msvc`。选它而非 `LibRaw.sln`／MSBuild 有两个具体原因：

1. `buildfiles/libraw.vcxproj` 钉死 **Windows SDK `10.0.18362.0`** 与平台工具集 **`v142`**（VS 2019）。用较新的 Visual Studio 构建会直接报 `MSB8036`／`MSB8020`。
2. 即使覆盖了上述两项，仍会失败于 `MSB6001`：当进程环境块里**同时存在大小写两份代理变量**（如 `HTTP_PROXY` 与 `http_proxy`）时，MSBuild 用于构造子进程环境的 .NET 字典不区分大小写，遇到重复键直接抛异常。PowerShell 的 `Env:` 提供程序会把两者规范成一条，因此**从 PowerShell 里删不掉**。

`nmake` 是纯 Win32 工具，不走 .NET，不受第 2 条影响；改用 `Makefile.msvc` 也完全绕开了第 1 条。

## 开发流程

本项目使用 [OpenSpec](https://github.com/Fission-AI/OpenSpec) 做规范驱动开发（`schema: spec-driven`）：

- 规划产物位于 `openspec/changes/<change-name>/`，包含 `proposal.md`、`design.md`、`specs/`、`tasks.md`
- 每个变更对应一个**同名 git 分支**（kebab-case，例如 `add-picture-control-parsing`）
- **变更完成后合入 `main`**，并用 `openspec archive` 把 spec delta 并入 `openspec/specs/`
- `main` 上只保留已归档的 specs 与规划产物，实现工作都在变更分支上进行

```bash
openspec list                          # 查看进行中的变更
openspec validate <change-name>        # 校验规划产物
openspec status --change <change-name> # 查看产物完成度
```

## 精度预期

如实说明：把相机观感搬到第三方 RAW 处理流程中**无法做到逐像素一致**——相机的动态范围处理、降噪、镜头校正等是独立于色彩配方的模块，其中一部分（尤其是**局部的**明暗处理）在数学上无法被一个逐像素的色彩变换表达。

可验证的目标是**视觉一致**，并用 ΔE 量化报告误差。参考量级：ΔE00 约 2.3 为可察觉阈值。

## 免责声明

非尼康官方项目，与尼康公司无任何关联。Nikon、NEF 等为尼康公司商标。

本项目处理的**第三方 Picture Control 配方是其各自作者的创作成果**。本仓库不分发配方文件本身；从配方推导出的色彩描述是否可再分发，取决于原配方的授权条款，使用者需自行确认。

## License

本仓库自身为 **MIT OR Apache-2.0** 双许可，可任选其一：

- [LICENSE-MIT](LICENSE-MIT)
- [LICENSE-APACHE](LICENSE-APACHE)

### 第三方组件

**构建出的二进制包含第三方代码，其许可与本仓库自身不同。**

| 组件 | 许可 | 源码位置 |
|---|---|---|
| LibRaw 0.22 | **CDDL-1.0** | [`third_party/libraw/`](third_party/libraw/)（随仓库分发） |

CDDL-1.0 是**文件级** copyleft：把它链接进本仓库不改变本仓库自身代码的许可，也不要求"可重新链接"；但**分发二进制时须随附 LibRaw 的源码**。本仓库已将其源码入库（含 `LICENSE.CDDL` 与 `COPYRIGHT`），因此由本仓库构建并分发的二进制满足该要求。

完整的组件清单、许可全文位置与合规说明见 [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md)。
