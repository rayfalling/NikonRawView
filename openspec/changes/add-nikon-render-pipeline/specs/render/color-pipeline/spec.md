## Purpose

把相机线性 RGB 转入工作空间线性 RGB，使后续的配方施加发生在正确的色彩空间里，而不是在相机原生空间或显示空间里凑合。

## ADDED Requirements

### Requirement: 按光源插值相机矩阵

系统 SHALL 依据照片的白平衡在两个已知光源的相机色彩矩阵之间插值，MUST NOT 固定使用单一光源的矩阵。

#### Scenario: 日光条件

- **WHEN** 照片的白平衡接近日光
- **THEN** 采用的矩阵接近 D65 光源下的矩阵

#### Scenario: 白炽灯条件

- **WHEN** 照片的白平衡明显偏暖
- **THEN** 采用的矩阵向 A 光源下的矩阵偏移，而非仍用 D65 矩阵

#### Scenario: 插值可解释

- **WHEN** 报告实际采用的矩阵
- **THEN** 同时报告两个光源的原始矩阵与插值权重，使结果可复核

### Requirement: 白平衡来自文件

系统 SHALL 使用 NEF 中记录的拍摄白平衡系数。系统 MUST NOT 假定日光白平衡。

#### Scenario: 非日光白平衡

- **WHEN** 照片以自定义或非日光白平衡拍摄
- **THEN** 输出反映该白平衡，而非被强制回日光

#### Scenario: 可覆盖

- **WHEN** 调用方显式指定白平衡
- **THEN** 该指定覆盖文件中的值，且输出中标明使用的是指定值

### Requirement: 转入 ProPhoto 线性工作空间

系统 SHALL 把结果转入 ProPhoto 线性，与 Lightroom Classic / Camera Raw 的内部工作空间一致，以便后续交付无需二次转换。

#### Scenario: 输出空间

- **WHEN** 管线完成
- **THEN** 输出为 ProPhoto 线性，且该事实可被调用方查询

#### Scenario: 与交付环节一致

- **WHEN** 后续导出 TIFF 时嵌入 ProPhoto ICC
- **THEN** 像素值无需再转换——管线的输出空间即交付空间

### Requirement: 不提前裁切色域

系统 SHALL NOT 在工作空间阶段裁切超出色域的成分。色域映射 MUST 推迟到最终输出编码时进行。

#### Scenario: 高饱和颜色

- **WHEN** 图像含超出 sRGB 且部分超出 ProPhoto 的高饱和颜色
- **THEN** 管线各阶段不得提前裁切，使其在最终输出阶段仍可被正确映射

#### Scenario: 负值保留

- **WHEN** 色彩变换产生负值分量
- **THEN** 负值在中途 MUST NOT 被截断为 0，除非是最终编码阶段
