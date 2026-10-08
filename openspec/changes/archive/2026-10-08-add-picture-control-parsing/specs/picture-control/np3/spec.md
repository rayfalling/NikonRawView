## Purpose

解析 NP3/NCP 配方容器，取出完整配方参数（含 257 级色调曲线、8 段色相混合器、3 区校色），使上层能在不依赖尼康工坊的前提下复现该配方的观感。

## ADDED Requirements

### Requirement: 容器识别

系统 SHALL 以 ASCII magic `NCP\0`（`4E 43 50 00`）识别配方容器，并读出家族码与版本。

#### Scenario: 识别配方文件

- **WHEN** 输入以 `NCP\0` 开头的文件
- **THEN** 系统报告家族码（`0310` 或 `0300`）与版本字段

#### Scenario: 非配方文件

- **WHEN** 输入不以 `NCP\0` 开头的文件
- **THEN** 系统返回明确的"非配方容器"结果，而非错误

### Requirement: 按 chunk 遍历解析

系统 SHALL 自偏移 `0x10` 起按 `[u32BE tag][u32BE len][payload]` 逐块遍历，遇到 `tag = 0` 且 `len = 0` 结束。系统 MUST NOT 以固定偏移定位任何字段。

#### Scenario: 两代模板均可解析

- **WHEN** 分别解析 `0310` 代（32 个 chunk，约 978 字节）与 `0300` 代（22–23 个 chunk，约 834/904 字节）的文件
- **THEN** 两者均解析成功，且各自的曲线块偏移被正确识别

#### Scenario: 固定偏移被拒绝

- **WHEN** 实现试图以固定偏移 `0x1CC` 读取曲线 LUT
- **THEN** 该实现在 `0300` 代上必然失败（其曲线块位于 `0xF4` 或 `0x13A`），因此该做法被规范禁止

### Requirement: 解析配方参数

系统 SHALL 取出配方名称、base 编码、标量参数、色相混合器与校色。

#### Scenario: 名称

- **WHEN** 解析 `0x0200` 块
- **THEN** 得到 NUL 结尾的配方的显示名

#### Scenario: 色相混合器

- **WHEN** 解析 `0x1F00` 块（长度 28）
- **THEN** 得到恰好 8 段，顺序为 R、O、Y、G、C、B、P、M，每段含 hue / chroma / brightness 三个值

#### Scenario: 校色

- **WHEN** 解析 `0x2000` 块（长度 20）
- **THEN** 得到 highlights / mid-tone / shadows 三区，以及 blending 与 balance 两个公共值

#### Scenario: base 编码

- **WHEN** 解析 `0x0300` 块的 2 字节载荷
- **THEN** 得到基准色彩编码，且系统能将其与 NEF 侧的基准名字对应（`0x03C2` ↔ `NEUTRAL`，`0x0020` ↔ `FLEXIBLE COLOR`）

### Requirement: 解析色调曲线

系统 SHALL 识别曲线块（`tag = 0x00000002`，`len = 578`），读出控制点数量与 257 级 16 位大端 LUT。

#### Scenario: 含曲线的配方

- **WHEN** 配方含曲线块
- **THEN** 曲线块载荷以 `BI0\0` 标记开头，系统取出控制点数量与 257 级 LUT，并报告 LUT 是否单调

#### Scenario: 不含曲线的配方

- **WHEN** 配方不含曲线块
- **THEN** 系统报告"无自定义曲线"，且不将其视为解析失败

### Requirement: 不以字节同一性判定配方等价

系统 SHALL 以解析出的参数判断两份配方是否等价。系统 MUST NOT 以文件哈希或字节长度判等。

#### Scenario: 同一配方的两代存储

- **WHEN** 同一配方分别以 `0300` 代（904 字节）与 `0310` 代（978 字节）存在
- **THEN** 系统识别为同名配方，并分别报告各自代次与参数差异

### Requirement: 报告未知块

系统 SHALL 保留未识别 tag 的原始载荷并在结果中报告，使后续可实现无需改动解析器主体即可读取新字段。

#### Scenario: 遇到未知 tag

- **WHEN** 遍历到规范未定义的 tag
- **THEN** 该块的 tag、长度与原始载荷被保留并报告，解析继续
