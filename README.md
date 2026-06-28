# FPM-Experiment

傅里叶叠层显微术（Fourier Ptychographic Microscopy, FPM）实验复现项目。当前仓库用于整理 FPM 的实验数据处理、LED 光强标定、重建算法、仿真验证、阅读笔记和阶段报告。

## 项目定位

本项目围绕 LED 阵列照明显微成像的 FPM 复现展开，重点包括：

- 低分辨率图像序列的读取、裁剪、HDR/曝光归一化和数据封装；
- 基于交替投影思想的 FPM 振幅/相位重建；
- LED 阵列位置、照明角度和光强非均匀性的建模与修正；
- 仿真样本、分辨率测试和非理想前向模型验证；
- FPM 和 DPC 相关论文的阅读笔记、翻译笔记和报告资料。

当前实验基本参数：

| 参数 | 当前值 |
| --- | --- |
| Objective NA | 0.15 |
| Wavelength | 628 nm |
| 原始放大倍数 | 5x |
| LED-sample distance | 约 4 cm |
| LED pitch | 4 mm |

`Code/Data` 中的公开样例数据来自原始 FPM 示例数据，参数见 [Code/data_description.txt](Code/data_description.txt)。本项目自己的 0619 实验批次参数集中在 [Code/function/fpm_0619_shared_config.m](Code/function/fpm_0619_shared_config.m)。

## 目录结构

```text
FPM-Experiment/
├── Code/                         # MATLAB/Python 代码和小型样例数据
│   ├── Main/                     # FPM 主重建入口与核心算法
│   ├── Raw_data_processing/      # 原始图像读取、裁剪、HDR/曝光预处理
│   ├── Diffuser_calibration/     # 毛玻璃/扩散片 LED 光强标定
│   ├── Led_intensity_correction_0529_legacy/
│   ├── function/                 # LED 几何、波矢、Zernike、共享配置等函数
│   ├── Data/                     # 小型样例数据
│   ├── Raw_input/                # 运行期原始采集图像和预处理输入，默认不纳入公开同步
│   ├── Results/                  # 运行期重建和标定输出，默认不纳入公开同步
│   └── 模拟测试/                 # 仿真、DPC 和分辨率 benchmark
├── Results/                      # 根目录本地结果归档，默认不纳入公开同步
├── Notes/                        # 论文笔记、实验笔记、翻译草稿
├── Original Paper/               # 原始论文 PDF 本地资料
├── Reference/                    # 综述、方法拓展和参考论文
├── Differential Phase Contrast (DPC)/
└── Report/                       # 阶段报告和文档源文件
```

更细的使用说明见各目录下的 `README.md`。

## 快速开始

推荐环境：

- MATLAB R2021a 或更新版本；
- Image Processing Toolbox；
- Python 3.10+，仅用于 `Code/模拟测试` 中的辅助报告脚本。

MATLAB 运行方式：

1. 打开 MATLAB，并将工作目录切换到仓库根目录或 `Code/Main`。
2. 运行主重建入口：

   ```matlab
   run("Code/Main/FP_recover_code.m")
   ```

3. 如需切换实验批次或参数 preset，优先修改：

   ```matlab
   Code/function/fpm_0619_shared_config.m
   ```

4. 如需从原始采集图像重新生成 `.mat` 输入，先运行：

   ```matlab
   run("Code/Raw_data_processing/read_transfer_0619.m")
   ```

注意：`FP_recover_code.m` 当前默认读取 0619 批次预处理输出，例如 `FPM_RawData_z z=_0.050_7_preproc_v1.mat`。该类原始输入和重建结果通常体积较大，默认作为本地实验资产管理。

## 数据与同步策略

为了让 GitHub 主仓库保持可读、可复现、可维护，本仓库采用以下策略：

- 源码、配置、说明文档和小型样例数据可以进入主仓库；
- `Results/` 和运行期可能生成的 `Code/Results/` 中的重建 `.mat`、标定 `.mat`、渲染图和中间结果默认不进入主仓库；
- `Code/Raw_input/` 中的相机原始采集图像默认不进入主仓库；
- 论文 PDF、Word 报告和翻译版 PDF 的版权状态需要单独确认，公开同步时建议只提交索引和阅读笔记；
- 大体积公开数据建议使用 GitHub Releases、Zenodo、OSF 或网盘，并在 README 中给出下载链接和校验信息。

相关边界说明见 [NOTICE.md](NOTICE.md)。

## GitHub 同步

远端仓库：

```bash
git remote add origin https://github.com/janxu2417/FPM-Experiment.git
```

如果远端已经存在：

```bash
git remote set-url origin https://github.com/janxu2417/FPM-Experiment.git
```

推荐同步流程：

```bash
git status --short
git add README.md NOTICE.md .gitignore .gitattributes .editorconfig
git add Code/**/README.md Results/README.md Notes/README.md "Original Paper/README.md" Reference/README.md "Differential Phase Contrast (DPC)/README.md" Report/README.md
git commit -m "docs: organize FPM experiment repository"
git push origin main
```

如果需要同步源码重构或新增代码，请先用 `git status --short` 检查将要提交的文件，避免把 `Results`、`Code/Results`、原始图像、论文 PDF 或 Word 报告误加入提交。

## 项目状态

这是一个实验复现和学习型研究仓库。代码中保留了阶段性脚本、legacy 版本和实验探索记录；正式复现实验时建议优先从 `Code/Main/FP_recover_code.m`、`Code/Raw_data_processing/read_transfer_0619.m` 和 `Code/function/fpm_0619_shared_config.m` 进入。
