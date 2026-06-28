# Code/Transfer_legacy

旧版数据格式转换脚本目录。

## 文件

| 文件 | 作用 |
| --- | --- |
| `HDR__14.m` | 早期 HDR 合成/转换脚本 |
| `png_to_mat.m` | PNG 到 MAT 的转换脚本 |

## 使用建议

该目录主要用于保留历史处理路径。新的实验数据建议优先使用：

```text
Code/Raw_data_processing/
```

如需继续使用 legacy 脚本，应在运行前检查路径、文件命名、图像尺寸、曝光参数和输出字段是否与当前 `himrecover` 接口一致。
