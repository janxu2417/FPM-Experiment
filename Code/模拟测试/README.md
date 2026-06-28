# Code/模拟测试

FPM/DPC 仿真和 benchmark 目录。

## 主要文件

| 文件 | 作用 |
| --- | --- |
| `FPM_Benchmark.m` | FPM 分辨率和重建效果 benchmark |
| `FPM_resolution_benchmark.m` | 分辨率 benchmark 辅助脚本 |
| `fpm_measure_resolution_groups.m` | USAF 等分辨率目标的线组识别/测量 |
| `FP_image_sim.m` | FPM 图像形成仿真 |
| `FP_test_recover.m` | 仿真数据恢复测试 |
| `benchmark_nonideal_forward.py` | 非理想前向模型 benchmark |
| `build_simulation_report.py` | 仿真报告生成辅助脚本 |
| `DPC.m`、`Differential phase contrast.m` | DPC 相关测试脚本 |

## 输出管理

仿真生成的报告图、`.mat` 结果、`.ppm` 渲染文件和压缩包默认不提交到主仓库。需要公开时，建议只选择关键图像和指标，并在 README 或报告中说明：

- 仿真目标；
- 物镜 NA、照明 NA、波长和采样；
- 是否加入 LED 位置误差、强度非均匀、噪声、离焦或 pupil aberration；
- benchmark 使用的评价指标。

## 物理解释提醒

仿真 benchmark 的关键不是只看图像锐度，而是同时检查频谱覆盖、明暗场过渡、线剖面、相位稳定性和强度归一化。对于 FPM，照明角度错误和光强校正错误常会在高频拼接处表现为伪影或相位偏置。
