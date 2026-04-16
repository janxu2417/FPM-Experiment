%% FPM 模拟与实验数据对比诊断 (FPM Benchmarking)
clear; close all; clc;

% =========================================================================
% 1. 加载数据
% =========================================================================
fprintf('加载模拟与实验数据...\n');
load('FP_input_data_Sim.mat', 'imseqlow1'); % 载入模拟数据: imseqlow1
load('FPM_RawData_3HDR.mat', 'imlow_HDR');       % 载入实验数据: imlow_HDR

% 确保两者帧数一致
num_LED = min(size(imseqlow1, 3), size(imlow_HDR, 3));
imseqlow1 = imseqlow1(:, :, 1:num_LED);
imlow_HDR = imlow_HDR(:, :, 1:num_LED);

% 提取尺寸
[m, n, ~] = size(imseqlow1);
[me, ne, ~] = size(imlow_HDR);

% 寻找中心明场 LED（通常是最亮的一张）
% 通过计算每帧的平均能量来定位中心 LED
energy_sim = squeeze(mean(mean(imseqlow1, 1), 2));
energy_exp = squeeze(mean(mean(imlow_HDR, 1), 2));
%[~, center_idx] = max(energy_sim);

% =========================================================================
% 2. 数据全局归一化 (基于物理相对光强)
% 为了对比，将模拟和实验数据都除以各自中心明场图像的最大值
% =========================================================================
norm_factor_sim = max(max(imseqlow1(:, :, 1)));
norm_factor_exp = max(max(imlow_HDR(:, :, 1)));

Sim_Norm = imseqlow1 / norm_factor_sim;
Exp_Norm = imlow_HDR / norm_factor_exp;

% =========================================================================
% 3. 可视化 1：代表性角度的空间形态对比 (Spatial Comparison)
% =========================================================================
for center_idx = 1:2:10
    % 挑选三个代表性角度：中心明场、明暗场交界、大角度暗场
    idx_bright = center_idx;
    idx_mid    = center_idx + 24; % 偏移几个位置，具体可根据你的扫描顺序调整
    idx_dark   = center_idx + 39;              % 通常第一帧是最外圈的暗场
    figure('Name', '空间形态与暗场细节对比', 'Position', [100, 100, 1200, 600]);
    
    % 明场对比
    subplot(2, 3, 1); imagesc(Sim_Norm(:,:,idx_bright)); colormap gray; axis image off; 
    title(sprintf('模拟 - 明场 (LED %d)', idx_bright));
    subplot(2, 3, 4); imagesc(Exp_Norm(:,:,idx_bright)); colormap gray; axis image off; 
    title('实验 - 明场');
    
    % 过渡区对比
    subplot(2, 3, 2); imagesc(Sim_Norm(:,:,idx_mid)); colormap gray; axis image off; 
    title(sprintf('模拟 - 过渡区 (LED %d)', idx_mid));
    subplot(2, 3, 5); imagesc(Exp_Norm(:,:,idx_mid)); colormap gray; axis image off; 
    title('实验 - 过渡区');
    
    % 暗场对比 (使用对数显示，或者极小的显示范围，以看清底噪)
    subplot(2, 3, 3); imshow(Sim_Norm(:,:,idx_dark), []); colormap gray; axis image off; 
    title(sprintf('模拟 - 暗场 (LED %d)', idx_dark));
    subplot(2, 3, 6); imshow(Exp_Norm(:,:,idx_dark), []); colormap gray; axis image off; 
    title('实验 - 暗场');
end
% =========================================================================
% 4. 可视化 2：中心明场线剖面分析 (Line Profile Analysis)
% 诊断过曝、Gamma压缩以及照明不均匀性(Vignetting)
% =========================================================================
figure('Name', '线剖面与光度学诊断', 'Position', [150, 150, 1000, 400]);

% 取中心水平线
line_rows = round(m / 2);
line_rowe = round(me / 2);
profile_sim = Sim_Norm(line_rows, :, idx_bright);
profile_exp = Exp_Norm(line_rowe, :, idx_bright);

subplot(1, 2, 1);
plot(1:n, profile_sim, 'b-', 'LineWidth', 1.5); hold on;
plot(1:ne, profile_exp, 'r-', 'LineWidth', 1.5);
grid on; legend('理论模拟', '实际实验');
title('中心明场水平剖面 (Row=Center)');
xlabel('Pixel'); ylabel('归一化强度');
% 物理诊断提示：
% 1. 如果红线在顶部被"削平"，说明相机硬件过曝。
% 2. 如果红线比蓝线更"鼓"或反差更大，说明相机Gamma被开启了(非线性)。
% 3. 如果红线整体呈现抛物线衰减，说明 LED 并非理想平面波，存在高斯光束衰减。

% =========================================================================
% 5. 可视化 3：全局光强衰减曲线 (Energy Decay Curve)
% 诊断系统底噪、杂散光与动态范围
% =========================================================================
subplot(1, 2, 2);
% 将各帧总能量归一化后画在对数坐标下
energy_sim_norm = energy_sim / max(energy_sim);
energy_exp_norm = energy_exp / max(energy_exp);

semilogy(1:num_LED, energy_sim_norm, 'b.-', 'LineWidth', 1); hold on;
semilogy(1:num_LED, energy_exp_norm, 'r.-', 'LineWidth', 1);
grid on; legend('理论模拟', '实际实验');
title('多角度 LED 能量衰减曲线 (Log Scale)');
xlabel('LED 索引 (通常对应不同照明角度)'); ylabel('归一化全图平均能量');
% 物理诊断提示：
% 1. 在大角度（暗场）区域，如果红线（实验）远远高于蓝线（模拟），平滑地下不去了，
%    说明相机的暗电流噪声（底噪）或环境环境漏光（杂散光）掩盖了真实的微弱高频信号。
% 2. 如果红线下降的斜率和蓝线完全不同，再次证明相机的线性度被破坏（如对比度增强未关）。
