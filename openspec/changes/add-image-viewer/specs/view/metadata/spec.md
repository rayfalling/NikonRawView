## Purpose

提取并展示相机型号与拍摄参数，使浏览时能直接看到"这是哪台机器、什么镜头、什么曝光拍的"，无需另开工具。

## ADDED Requirements

### Requirement: 显示机型与镜头

系统 SHALL 显示相机厂商与型号、以及镜头标识。

#### Scenario: 机型

- **WHEN** 打开一张 Nikon Z8 拍摄的 NEF
- **THEN** 界面显示可读的机型名（如 `Nikon Z 8`），而非仅厂商代码

#### Scenario: 镜头

- **WHEN** 文件中含镜头信息
- **THEN** 显示镜头标识；若文件只含镜头 ID，则显示该 ID 并标明未能解析为名称

### Requirement: 显示曝光参数

系统 SHALL 显示快门、光圈、ISO 与焦距。

#### Scenario: 完整参数

- **WHEN** 文件含完整的曝光信息
- **THEN** 显示快门（如 `1/500`）、光圈（如 `f/5.6`）、ISO 与焦距（如 `400mm`）

#### Scenario: 缺失字段

- **WHEN** 某一项在文件中不存在
- **THEN** 该项显示为明确的"未知"，MUST NOT 用 0 或空字符串冒充

#### Scenario: 数值格式化

- **WHEN** 展示快门速度
- **THEN** 慢于 1 秒的显示为分数形式，不显示为 `0.002` 一类难以速读的形式

### Requirement: 显示拍摄时间

系统 SHALL 显示拍摄时间，并 SHALL 处理文件中没有时区信息的情况。

#### Scenario: 含时间

- **WHEN** 文件含拍摄时间
- **THEN** 显示该时间

#### Scenario: 无时区信息

- **WHEN** 文件的时间字段不带时区
- **THEN** 按文件的本地时间原样显示并标明"无时区信息"，MUST NOT 擅自换算为浏览者所在时区

### Requirement: 取值为只读展示

系统 SHALL 保证这些字段仅用于展示，不参与且不影响渲染结果。

#### Scenario: 展示不改变像素

- **WHEN** 查看某张图的信息
- **THEN** 渲染结果与未查看时完全一致
