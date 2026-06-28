# Code

本目录保存 FPM 复现实验的代码和小型样例数据。

## 子目录

| 目录 | 作用 |
| --- | --- |
| `Main/` | FPM 主重建入口、核心恢复算法、结果显示和历史入口脚本 |
| `Raw_data_processing/` | 原始采集图像读取、裁剪、曝光/HDR 处理和 `.mat` 输入生成 |
| `Diffuser_calibration/` | 基于毛玻璃/扩散片图像的 LED 光强标定 |
| `Led_intensity_correction_0529_legacy/` | 0529 批次旧版 LED 光强修正脚本 |
| `function/` | 几何、波矢、Zernike、共享配置和 manifest 辅助函数 |
| `Data/` | 小型样例数据和示例图像 |
| `Raw_input/` | 运行期原始采集图像和预处理输入，默认不提交到 GitHub |
| `Results/` | 运行期重建和标定输出，默认不提交到 GitHub |
| `Transfer_legacy/` | 旧版格式转换脚本 |
| `模拟测试/` | FPM/DPC 仿真、分辨率 benchmark 和报告生成脚本 |

## 推荐主流程

1. 配置实验 preset：

   ```matlab
   Code/function/fpm_0619_shared_config.m
   ```

2. 从原始采集图像生成重建输入：

   ```matlab
   run("Code/Raw_data_processing/read_transfer_0619.m")
   ```

3. 运行 FPM 重建：

   ```matlab
   run("Code/Main/FP_recover_code.m")
   ```

4. 查看输出：

   ```text
   Code/Results/<batch>/
   ```

## 重要约定

- `fpm_0619_shared_config.m` 是 0619 批次的共享参数入口，预处理和重建都应优先从这里读取批次、光学、几何、标定和重建参数。
- `Results/`、运行期可能生成的 `Code/Results/` 和 `Code/Raw_input/` 主要是本地实验资产，不建议直接提交到主仓库。
- 新增实验批次时，优先新增 preset，而不是复制多份硬编码脚本。
