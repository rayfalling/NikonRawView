# nef/picture-control Specification

## Purpose
从 NEF 与 JPEG 的 Nikon MakerNote 中读出该照片所用 Picture Control 的身份信息，使上层无需依赖尼康工坊即可知道"这张照片用了哪个配方、基于哪个基准色彩"。

## Requirements

### Requirement: 定位配方数据

系统 SHALL 沿 TIFF/IFD 链（IFD0 → ExifIFD → MakerNote → tag `0x0023`）定位配方载荷。系统 MUST NOT 以字节模式匹配来定位 Nikon 元数据。

#### Scenario: NEF 文件

- **WHEN** 输入一个 Nikon Z8 生成的 NEF
- **THEN** 系统返回 tag `0x0023` 的载荷及其长度

#### Scenario: JPEG 文件

- **WHEN** 输入由同一相机生成的配对 JPEG
- **THEN** 系统以 APP1 段内 TIFF 头为基址解析所有偏移，返回与配对 NEF **逐字节相同**的载荷

#### Scenario: 缺少配方数据

- **WHEN** 目标文件不含 tag `0x0023`
- **THEN** 系统返回明确的"无配方数据"结果，而非错误或空值

### Requirement: 解析配方身份字段

系统 SHALL 按 PictureControl3 布局解析载荷，取出版本、自定义配方名、基准色彩名、调整模式、标量参数与滤镜/色调效果。

#### Scenario: 自定义配方

- **WHEN** 照片以自定义配方拍摄
- **THEN** 自定义名与基准名均被取出且互不相同（示例：`LINKS-Nature` / `NEUTRAL`）

#### Scenario: 内置配方

- **WHEN** 照片以内置配方拍摄
- **THEN** 自定义名与基准名相同（示例：`PORTRAIT` / `PORTRAIT`）

#### Scenario: 标量解码

- **WHEN** 解析锐化、中频锐化、清晰度
- **THEN** 按 `(b − 0x80) / 4` 解码

#### Scenario: 标量解码（线性字段）

- **WHEN** 解析对比度、亮度、饱和度、色相
- **THEN** 按 `b − 0x80` 解码

### Requirement: 报告曲线哨兵

系统 SHALL 报告对比度与亮度是否同时等于哨兵值 `0x01`，作为"该配方可能含自定义色调曲线"的提示。该标志 MUST NOT 被用作是否需要解析配方库的判据。

#### Scenario: 哨兵为真

- **WHEN** 对比度与亮度均为 `0x01`
- **THEN** 系统标记"可能存在自定义曲线"，并**继续**按配方名解析配方库

#### Scenario: 哨兵不构成短路

- **WHEN** 哨兵为假，或某配方的哨兵与其 NP3 文件中的曲线标志不一致
- **THEN** 系统 MUST NOT 据此跳过配方库解析，因为实测存在哨兵为假但配方确含曲线的反例

### Requirement: 不依赖前导头的可变字段

配方载荷前存在 8 字节前导头，其中间一个 u16 会随文件变化（实测取值 `0x01` 与 `0x03`）。系统 MUST NOT 依赖该字段识别配方或定位载荷。

#### Scenario: 两种前导取值

- **WHEN** 分别解析前导取值为 `0x01` 与 `0x03` 的两个文件
- **THEN** 两者的配方身份被正确解析，且结果不因该字段而不同
