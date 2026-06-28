# Results

本目录用于保存已经筛选后的重建结果图、标定诊断图、仿真结果图和实验进展记录。

当前运行期脚本默认输出路径是 `Code/Results/`。本目录作为手工整理后的轻量结果归档区，可进入 GitHub 主仓库，便于直接浏览复现实验效果。大型 `.mat`、压缩包和渲染中间文件仍由 `.gitignore` 排除。

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
