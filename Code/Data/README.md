# Code/Data

小型样例数据目录。

## 当前内容

该目录包含 FPM 示例 `.mat` 数据和少量示例图像。`Code/data_description.txt` 中记录了原始示例数据的采集条件：

- Nikon 2x objective，NA = 0.1；
- 样品面低分辨率像素尺寸 1.845 um；
- 32 x 32 LED 阵列，LED pitch = 4 mm；
- 15 x 15 LED 以 spiral-out 顺序点亮；
- LED 到样品距离 90.88 mm；
- 每个 `.mat` 文件包含 `imlow_HDR`、`theta`、`wlength`、`xint/yint`、`z` 等字段。

## 同步建议

适合随仓库保存的是小型、可公开、可复现实验的样例数据。新增数据前建议确认：

- 文件体积是否适合 Git；
- 数据来源和再分发权限是否清楚；
- README 或 manifest 是否记录了采集参数；
- 是否可以用外部数据链接替代直接提交二进制文件。
