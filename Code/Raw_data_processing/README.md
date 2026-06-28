# Code/Raw_data_processing

原始采集图像预处理目录。这里的脚本负责把相机采集的多 LED 图像序列整理成 FPM 重建入口可直接读取的 `.mat` 文件。

## 主要文件

| 文件 | 作用 |
| --- | --- |
| `read_transfer_0619.m` | 当前 0619 批次推荐预处理入口 |
| `read_transfer_0619_legacy_20260619.m` | 0619 批次旧版预处理脚本 |
| `read_transfer_0529.m` | 0529 批次预处理脚本 |
| `test_0507.m`、`test_0515.m` | 早期测试脚本 |

## 典型流程

1. 将原始采集图像放到：

   ```text
   Code/Raw_input/<batch>/<input_subdir>/
   ```

2. 在 `Code/function/fpm_0619_shared_config.m` 中确认：

   - `batch_folder`
   - `input_subdir`
   - `num_img`
   - `file_indices`
   - `exposure_ms`
   - ROI 配置
   - 波长、采样间隔、LED 几何参数

3. 运行：

   ```matlab
   run("Code/Raw_data_processing/read_transfer_0619.m")
   ```

4. 预处理输出写入：

   ```text
   Code/Raw_input/<batch>/*.mat
   ```

## 输出约定

重建脚本默认期待预处理输出中包含：

- `im_high_HDR`：预处理后的低分辨率振幅序列；
- `spsize`：样品面低分辨率像素尺寸；
- `wlength`：中心波长；
- `z`：离焦量；
- `save_manifest`：输入、参数、输出的记录信息。
