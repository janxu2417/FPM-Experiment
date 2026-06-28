# Results

本目录用于保存已经归档到项目根目录的本地重建结果、标定结果、仿真结果和诊断图。

当前运行期脚本默认输出路径是 `Code/Results/`。本目录可作为历史结果或手工整理结果的归档区。两类结果目录都默认被 `.gitignore` 忽略，只保留 README。原因是 `.mat` 结果文件和中间图像可能很大，且通常可以由代码和输入数据重新生成。

## 常见结构

```text
Results/
├── 0619/
│   ├── Recover_*.png
│   └── Recover_*.mat
├── Calibration/
│   └── 0619/
│       ├── diffuser_calib_*.mat
│       ├── diffuser_summary_*.png
│       └── diffuser_residuals_*.png
└── sample/
```

## 发布建议

如果某个结果需要随论文、报告或复现实验公开，建议：

- 优先发布压缩后的关键图像和小型统计指标；
- 对大体积 `.mat` 使用 GitHub Releases、Zenodo、OSF 或网盘；
- 在 README 或 manifest 中记录生成该结果的代码 commit、输入数据、preset 和运行日期。
