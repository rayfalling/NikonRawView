# picture-control/library Specification

## Purpose
发现本机可用的 Picture Control 配方并建立索引，使上层能由照片中的配方名解析出配方本体，并知道它来自哪里、属于哪一代格式。

## Requirements

### Requirement: 发现配方来源

系统 SHALL 至少覆盖两个本机来源：注册表 `HKCU\Software\Nikon\Common\PictureControl\NP3\CustomCurves`，以及 NX Studio 的 `.nksc` 边车文件。

#### Scenario: 从注册表读取

- **WHEN** 系统扫描注册表该键
- **THEN** 每个 `REG_BINARY` 值被作为一份配方载荷解析，配方名从载荷内部读取（而非取自注册表值名）

#### Scenario: 从边车读取

- **WHEN** 目录中存在 `.nksc` 文件
- **THEN** 系统从 XMP 中 `id="nikon::PictureControl"` 的滤镜里取出 Base64 编码的配方载荷并解析

#### Scenario: 边车内的多份载荷

- **WHEN** 边车同时含 `ExportData` 与 `CustomData`
- **THEN** 两份都被取出，并各自报告其在边车中的字段来源

#### Scenario: 来源不可用

- **WHEN** 注册表键不存在或目录无 `.nksc` 文件
- **THEN** 系统正常返回其余来源的结果，不报错

### Requirement: 按配方名解析配方本体

系统 SHALL 支持"给定配方名 → 返回配方本体"的解析。

#### Scenario: 命中

- **WHEN** 请求一个本机存在的配方名
- **THEN** 返回该配方的完整参数（含标量、色相混合器、校色、色调曲线）与来源标识

#### Scenario: 未命中

- **WHEN** 请求一个本机不存在的配方名
- **THEN** 返回明确的"未找到"结果，并列出本机可用的配方名，供诊断

#### Scenario: 名称大小写与空白

- **WHEN** 请求名与存储名仅大小写或首尾空白不同
- **THEN** 系统仍应命中，因为配方名在文件内为定长 NUL 填充字段

### Requirement: 报告来源与格式代次

系统 SHALL 对每份解析结果标注来源与容器代次。

#### Scenario: 同名多来源

- **WHEN** 同一配方名在注册表与边车中各存在一份
- **THEN** 两份结果都被返回，并各自标注来源与代次，由调用方决定采用哪一份

### Requirement: 索引可枚举

系统 SHALL 支持枚举本机全部可用配方，供用户在界面中选择。

#### Scenario: 枚举

- **WHEN** 调用方请求列出全部配方
- **THEN** 返回去重后的配方名集合，每个名字附带其可用来源与代次
