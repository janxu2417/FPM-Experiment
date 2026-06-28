# Code/Main

FPM 主重建入口和核心恢复算法目录。

## 主要文件

| 文件 | 作用 |
| --- | --- |
| `FP_recover_code.m` | 当前推荐的一键重建入口，读取共享 preset、加载预处理 `.mat`、调用 `himrecover` 并保存结果 |
| `himrecover.m` | FPM 高频谱拼接与 pupil/object 迭代更新核心函数 |
| `show_recover.m` / `ShowResult.m` | 结果显示辅助脚本 |
| `generate_ideal_target.m` | 仿真或 benchmark 中的理想目标生成 |
| `FP_3_HDR.m` / `FP_read_HDR11.m` | 历史 HDR/读取入口 |
| `*_legacy_*.m`、`himrecover_wrong.m` | 阶段性保留版本，用于对比和回溯 |

## 运行入口

从仓库根目录运行：

```matlab
run("Code/Main/FP_recover_code.m")
```

脚本会自动把 `Code/` 加入 MATLAB path，并读取：

```matlab
Code/function/fpm_0619_shared_config.m
```

默认 preset 当前指向 0619 批次的 7x7 LED 重建。若输入 `.mat` 不存在，需要先运行 `Code/Raw_data_processing/read_transfer_0619.m`。

## 算法说明

`himrecover` 的输入是低分辨率振幅序列、LED 归一化波矢、物镜 NA、波长、低/高分辨率采样间隔、离焦量和迭代参数。核心步骤为：

1. 用低分辨率图像和上采样比例初始化高分辨率复振幅；
2. 根据每个 LED 的照明角度截取高分辨率傅里叶谱局部区域；
3. 乘以 pupil/defocus 传递函数，反变换得到当前估计的低分辨率图像；
4. 用实测振幅替换估计振幅，保留相位；
5. 回到频域更新目标谱和可选 pupil；
6. 循环迭代，输出高分辨率振幅、相位和 pupil。
