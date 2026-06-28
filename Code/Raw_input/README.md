# Code/Raw_input

运行期原始采集图像和预处理输入目录。

`Code/function/fpm_0619_shared_config.m` 中当前设置：

```matlab
base.raw_root = "Raw_input";
```

因此 `Code/Raw_data_processing/read_transfer_0619.m` 会默认从以下路径读取原始图像：

```text
Code/Raw_input/<batch>/<input_subdir>/
```

并把预处理后的 `.mat` 输入写回：

```text
Code/Raw_input/<batch>/
```

## 推荐结构

```text
Code/Raw_input/
└── 0619/
    ├── z=_0.050_7/
    │   ├── 图像_1.tif
    │   ├── 图像_2.tif
    │   └── ...
    └── FPM_RawData_z z=_0.050_7_preproc_v1.mat
```

## 同步策略

该目录默认被 `.gitignore` 忽略，只保留本 README。原始图像和预处理 `.mat` 往往体积较大，应通过外部数据链接、GitHub Releases 或实验记录单独管理。
