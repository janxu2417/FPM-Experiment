# Code/Results

运行期重建和标定输出目录。

`Code/function/fpm_0619_shared_config.m` 中当前设置：

```matlab
base.result_root = "Results";
```

因此 `Code/Main/FP_recover_code.m` 和标定脚本会默认把输出写入：

```text
Code/Results/<batch>/
Code/Results/Calibration/<batch>/
```

## 常见内容

- `Recover_*.png`：重建振幅/相位预览图；
- `Recover_*.mat`：高分辨率复振幅、pupil、迭代记录和 manifest；
- `diffuser_calib_*.mat`：LED 相对光强标定结果；
- `diffuser_summary_*.png`、`diffuser_residuals_*.png`：标定诊断图。

## 同步策略

该目录用于运行期输出，建议定期筛选关键 PNG 图像复制到根目录 `Results/`。大型 `.mat`、压缩包和渲染中间文件不进入主仓库；需要公开完整数据时，建议发布到 GitHub Releases、Zenodo、OSF 或网盘，并记录代码 commit、输入数据和 preset。
