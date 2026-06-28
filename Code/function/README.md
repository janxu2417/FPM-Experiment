# Code/function

FPM 复现实验的共享函数目录。

## 主要文件

| 文件 | 作用 |
| --- | --- |
| `fpm_0619_shared_config.m` | 0619 批次预处理和重建共享配置入口 |
| `fpm_0619_build_manifest.m` | 生成预处理/重建 manifest，用于追踪输入和参数 |
| `fpm_sim_nonideal_config.m` | 非理想仿真的配置入口 |
| `LED_location.m` | 生成 LED 阵列坐标 |
| `k_vector.m` | 根据 LED 位置、阵列高度、玻璃参数等计算照明波矢 |
| `centerMatrix.m` | 矩阵中心裁剪/嵌入辅助函数 |
| `zernfun.m`、`gzn.m` | Zernike/像差相关函数 |

## 约定

- 新实验批次优先通过新增 preset 管理，而不是在多个脚本中重复硬编码参数。
- 预处理和重建都应尽量记录 manifest，便于复现实验时追踪数据来源。
- 几何参数的单位要保持一致：当前配置中 `H` 和 `LEDp` 使用 mm，波长和像素尺寸使用 m。
