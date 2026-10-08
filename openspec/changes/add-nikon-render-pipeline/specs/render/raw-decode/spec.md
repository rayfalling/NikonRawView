## Purpose

把 NEF 解成线性 RGB，并取出后续色彩变换所需的全部标定元数据，使上层不必依赖尼康工坊即可拿到可运算的像素。

## ADDED Requirements

### Requirement: 取出标定元数据

系统 SHALL 从 NEF 中取出白平衡系数、相机色彩矩阵（按光源）、白电平、黑电平、CFA 排布与有效像素区。相机色彩矩阵 MUST 按光源分别取出，而非只取一个。

#### Scenario: Z8 的 NEF

- **WHEN** 解析一个 Nikon Z8 生成的 NEF
- **THEN** 系统报告白平衡系数、至少两个光源下的相机色彩矩阵、白电平、黑电平、CFA 排布与裁切区

#### Scenario: 数值与独立实现一致

- **WHEN** 取出的白平衡与尺寸与另一套独立实现（如 LibRaw）比对
- **THEN** 两者一致；这是解码正确性的交叉验证方式

### Requirement: 输出线性 RGB

系统 SHALL 在施加白平衡与扣除黑电平后做去马赛克，输出线性 RGB。输出 MUST NOT 含 gamma 曲线、色调曲线或任何显示编码。

#### Scenario: 全画幅输出尺寸

- **WHEN** 解码一个 45 MP 的 NEF
- **THEN** 输出为 `8256×5504` 的线性 RGB，与相机标称有效像素一致

#### Scenario: 线性性

- **WHEN** 检查输出的数值分布
- **THEN** 输出为线性光，不含 gamma 或色调映射

#### Scenario: 不施加相机内观感

- **WHEN** 解码一张以自定义配方拍摄的 NEF
- **THEN** 输出 MUST NOT 含该配方的任何影响——配方的施加是渲染阶段的职责

### Requirement: 去马赛克质量

系统 SHALL 抑制去马赛克的典型伪像。系统 MUST NOT 移植 GPL 许可的实现。

#### Scenario: 高频细节

- **WHEN** 对含高频细节（如细密织物、远处枝叶）的区域去马赛克
- **THEN** 不出现明显的拉链效应，且不出现成片伪色

#### Scenario: 许可合规

- **WHEN** 实现去马赛克
- **THEN** 采用的算法来源为许可兼容者，且来源在代码中标注

### Requirement: 不支持的相机须明确失败

系统 SHALL 在不认识的机型上返回明确错误。系统 MUST NOT 猜测参数后继续。

#### Scenario: 未知机型

- **WHEN** 输入一个解码器不支持的 RAW
- **THEN** 返回指明机型与原因的失败，而非产出错误图像
