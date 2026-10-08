## Purpose

把「基准配方 + 配方」施加到线性 RGB，产出与相机内观感一致的图像；这是解析层取出的配方参数真正变成像素的地方。

## ADDED Requirements

### Requirement: 先基准、后配方

系统 SHALL 先施加基准配方的渲染变换，再施加配方相对基准的增量。系统 MUST NOT 把配方当作独立于基准的完整变换。

#### Scenario: 自定义配方

- **WHEN** 渲染一张身份为 `LINKS-Nature` / 基准 `NEUTRAL` 的照片
- **THEN** 先施加 `NEUTRAL` 的基准变换，再施加 `LINKS-Nature` 的增量

#### Scenario: 内置配方

- **WHEN** 渲染一张身份为 `PORTRAIT` / 基准 `PORTRAIT` 的照片
- **THEN** 基准与配方同名，增量部分为空，输出等于基准变换的结果

#### Scenario: 基准名与配方名不同却指向同一基准

- **WHEN** 两个配方的基准名不同但映射到同一编码
- **THEN** 系统按基准**编码**而非名称字符串选择基准变换

### Requirement: 施加色调曲线

系统 SHALL 施加配方中的 257 级色调曲线。

#### Scenario: 含曲线的配方

- **WHEN** 配方含非恒等曲线
- **THEN** 曲线被施加，且可用一条已知输入验证其映射结果

#### Scenario: 不含曲线

- **WHEN** 配方不含曲线块，或曲线为恒等
- **THEN** 该步骤被跳过且不报错

#### Scenario: 曲线单调性

- **WHEN** 施加一条已解析出的曲线
- **THEN** 系统不得假设曲线单调；非单调曲线必须被原样施加

### Requirement: 施加色相混合器

系统 SHALL 按 8 段施加色相混合器，段序为 Red、Orange、Yellow、Green、Cyan、Blue、Purple、Magenta，各段含色相、彩度、明度三个量。

#### Scenario: 单段调整

- **WHEN** 配方只调整了某一段
- **THEN** 只有该段对应的色相范围被改变，其余段基本不受影响

#### Scenario: 中性混合器

- **WHEN** 8 段全为中性值
- **THEN** 输出与未施加混合器时一致

### Requirement: 施加校色

系统 SHALL 按 highlights、mid-tone、shadows 三区施加校色，并应用其 blending 与 balance 两个公共参数。

#### Scenario: 三区分别生效

- **WHEN** 配方对三个区设置了不同的量
- **THEN** 亮部、中间调、暗部分别按其量被调整

#### Scenario: 中性校色

- **WHEN** 三区全为中性值
- **THEN** 输出与未施加校色时一致

### Requirement: 施加标量参数

系统 SHALL 按各标量的语义施加锐化、中频锐化、清晰度、对比度、亮度、饱和度与色相。

#### Scenario: 哨兵值必须被识别为不适用

- **WHEN** 对比度或亮度取哨兵值
- **THEN** 系统 MUST 将其视为「该参数不适用」并跳过，MUST NOT 当作 `−127` 施加

#### Scenario: 四分之一步进

- **WHEN** 施加锐化、中频锐化或清晰度
- **THEN** 按四分之一步进解释其数值

### Requirement: 渲染是纯函数

系统 SHALL 保证渲染不修改输入数据，且相同输入产生相同输出。

#### Scenario: 可重复

- **WHEN** 对同一输入渲染两次
- **THEN** 两次输出逐位相同

#### Scenario: 输入不被修改

- **WHEN** 渲染完成后检查输入缓冲
- **THEN** 输入未被改动

#### Scenario: 管线顺序固定

- **WHEN** 查询渲染的步骤顺序
- **THEN** 系统报告确定的、可复核的处理顺序，而非隐式依赖调用次序
