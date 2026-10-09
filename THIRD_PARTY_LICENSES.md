# 第三方组件许可

本文件列出**构建产物中包含**的第三方组件。本仓库自身的许可（MIT OR Apache-2.0）**不覆盖**这些组件。

> 本仓库的 crate 依赖为空——RAW 解码以外的部分全部自行实现。因此本清单当前只有一个组件。

---

## LibRaw 0.22.0

| 项 | 值 |
|---|---|
| 用途 | RAW 解码（NEF 解析、去马赛克） |
| 许可 | **CDDL-1.0**（Common Development and Distribution License） |
| 版权 | 见 [`third_party/libraw/COPYRIGHT`](third_party/libraw/COPYRIGHT) |
| 许可全文 | [`third_party/libraw/LICENSE.CDDL`](third_party/libraw/LICENSE.CDDL) |
| 上游 | <https://www.libraw.org/> ・ <https://github.com/LibRaw/LibRaw> |
| 源码位置 | [`third_party/libraw/`](third_party/libraw/)（随本仓库分发，未做修改） |
| 链接方式 | 静态链接（`libraw_static.lib`），由 `build.rs` 从上述源码构建 |

### 为什么选 CDDL 而不是 LGPL

LibRaw 是**三重许可**：LGPL-2.1 / CDDL-1.0 / 商业许可，使用者可任选其一。本项目选 **CDDL-1.0**，因为它对本项目的约束最小：

- CDDL 是**文件级** copyleft——它约束的是 LibRaw 自身的文件，不传染到链接它的本项目代码
- CDDL **没有** LGPL 那样的"必须允许用户重新链接"要求
- 义务集中在一点：分发包含 LibRaw 的二进制时，须提供 LibRaw 的**对应源码**

### 选型时评估过的替代方案

记录在此以便日后复核：

| 方案 | 未采用的原因 |
|---|---|
| `rawler` / `rawloader` / `quickraw` | 均为 **LGPL-2.1**；且 `rawler` 强制依赖 `tokio`(full)、`image`、`exr`，为读拜耳阵列拖入 129 个 crate |
| `imagepipe` | LGPL-3.0-only |
| `libraw-sys` / `libraw-rs` | 绑定本身是 MIT/Apache，但经查证 `libraw-sys` 不打包 LibRaw 源码，依赖 `pkg-config` 找系统库，无法直接使用 |
| 自行实现尼康解码 | 工作量与风险远超其价值；且相机色彩矩阵不在文件里，自行实现同样需要外部机型标定数据 |

---

## 分发二进制时的合规要点

由本仓库源码构建并分发的二进制包含 LibRaw。分发者需：

1. **随附 LibRaw 的源码**——即 `third_party/libraw/` 的全部内容。只要随发布产物一并提供（或提供等效的获取途径），即满足 CDDL 的要求。
2. **保留版权与许可声明**——即 `third_party/libraw/COPYRIGHT` 与 `LICENSE.CDDL`。
3. **若要修改 LibRaw**，需在修改的文件中注明改动，且这些文件仍受 CDDL 约束。

CDDL **不要求**本项目的代码改用 CDDL，也**不要求**提供重新链接的能力。

---

## 维护约定

- **新增任何第三方组件时，必须同步更新本文件**，至少写明：名称与版本、许可、许可全文的位置、源码位置、链接方式。
- 若引入的是 crate 依赖，还必须说明为何值得引入——本项目此前刻意保持零 crate 依赖。
- `openspec/changes/add-nikon-render-pipeline/tasks.md` 的任务 1.4 是本次补全这些文件的来源。
