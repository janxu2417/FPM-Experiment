# Code/Diffuser_calibration

LED 光强标定目录。当前方法使用毛玻璃/扩散片图像估计不同 LED 的相对照明强度，并为 FPM 重建提供外部光强修正因子。

## 主要文件

| 文件 | 作用 |
| --- | --- |
| `run_diffuser_calibration.m` | 标定主入口 |
| `build_diffuser_rois.m` | 构建 LED 标定 ROI |
| `fit_led_intensity_smooth.m` | 对 LED 强度进行平滑/多项式拟合 |
| `compare_fit_models.m` | 比较不同拟合模型 |
| `summarize_led_qc.m` | 汇总标定质量控制指标 |
| `debug_diffuser_first25_images.m` | 早期 25 图调试脚本 |
| `analyze_centriod.m` | LED/ROI 质心分析 |

## 输出

标定结果通常写入：

```text
Code/Results/Calibration/<batch>/
```

其中 `.mat` 文件保存相对光强、拟合模型和诊断信息，`.png` 文件用于检查残差、ROI 和拟合曲线。该目录默认不进入 GitHub 主仓库。

## 与重建的关系

`Code/Main/FP_recover_code.m` 中的外部 LED 光强修正通过共享配置控制：

```matlab
cfg.calibration.use_external_led_intensity_correction
cfg.calibration.calib_file_name
```

物理上，标定得到的是相对强度 `I_rel`，而重建输入使用振幅，因此振幅修正应使用：

```matlab
A_corr = A_raw / sqrt(I_rel)
```
