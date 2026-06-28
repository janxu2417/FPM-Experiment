# Code/Led_intensity_correction_0529_legacy

0529 批次 LED 光强修正的旧版脚本目录。

## 文件

| 文件 | 作用 |
| --- | --- |
| `run_fit_led_intensity_0529.m` | 0529 光强拟合入口 |
| `led_roi_intensity_0529.m` | ROI 强度提取 |
| `fit_led_intensity_piecewise.m` | 分段拟合模型 |

## 使用建议

该目录保留用于回溯早期处理流程。新实验批次建议优先使用：

```text
Code/Diffuser_calibration/
Code/function/fpm_0619_shared_config.m
```

如果需要复用 0529 方法，应先把输入路径、LED 阵列大小、曝光时间、ROI 和输出命名改为显式 preset，避免与 0619 主流程混用。
