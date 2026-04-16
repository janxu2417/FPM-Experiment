%% 微分相衬 (DPC) 定量相位重建 MATLAB 示例
clear; clc; close all;

% 1. 系统参数设置
lambda = 0.513;       % 波长 (um)
NA_obj = 0.4;         % 物镜数值孔径
NA_illum = 0.4;       % 照明数值孔径 (sigma = 1)
pixel_size = 0.325;   % 像面实际像素尺寸 (um) (CCD像素/放大倍率)
N = 512;              % 图像大小
alpha = 1e-2;         % 蒂霍诺夫正则化参数

% 2. 坐标系建立
x = (-(N/2):(N/2)-1) * pixel_size;
[X, Y] = meshgrid(x, x);
df = 1 / (N * pixel_size);
f_axis = (-(N/2):(N/2)-1) * df;
[fx, fy] = meshgrid(f_axis, f_axis);
fr = sqrt(fx.^2 + fy.^2);

% 3. 生成光瞳函数 P(u) 和光源 S(u)
% 光瞳函数 (无像差理想系统)
P = double(fr <= NA_obj / lambda);

% 生成四个半圆光源 (Top, Bottom, Left, Right)
S_top = double((fr <= NA_illum / lambda) & (fy > 0));
S_bottom = double((fr <= NA_illum / lambda) & (fy < 0));
S_left = double((fr <= NA_illum / lambda) & (fx < 0));
S_right = double((fr <= NA_illum / lambda) & (fx > 0));

% 4. 计算相位传递函数 Hph (基于式 8)
% 定义计算单轴 WOTF 的匿名函数
% H(u) = i/B * [ (S*P) \otimes P - P \otimes (S*P) ] (\otimes 为相关运算)
calc_Hph = @(S, P) 1i * (ifftshift(ifft2(fft2(fftshift(S.*P)) .* conj(fft2(fftshift(P))))) - ...
                        ifftshift(ifft2(fft2(fftshift(P)) .* conj(fft2(fftshift(S.*P))))));

B = sum(S_top(:) .* P(:)); % 背景项 B (式 6)

% 计算两个轴的传递函数
H_y = calc_Hph(S_top - S_bottom, P) / B;
H_x = calc_Hph(S_right - S_left, P) / B;

% 5. 模拟生成待测样本 (相位物体)
% 读取或生成一个测试图像
phase_true = 0.5 * phantom('Modified Shepp-Logan', N); 
% 弱物体假设 o = 1 + i*phi
object = exp(1i * phase_true); 

% 6. 模拟 DPC 图像采集 (式 1, 2)
% 此处简化模拟：直接通过传递函数在频域合成，或进行空间域卷积
% 为准确起见，此处使用线性模型模拟：I_dpc = F^-1( Hph * F(phi) )
I_dpc_y = real(ifft2(fft2(phase_true) .* fftshift(H_y)));
I_dpc_x = real(ifft2(fft2(phase_true) .* fftshift(H_x)));

% 添加少量噪声
I_dpc_y = I_dpc_y + 0.01 * randn(N);
I_dpc_x = I_dpc_x + 0.01 * randn(N);

% 7. 定量相位重建 (式 13)
% 频域合成重建
numerator = conj(H_x) .* fftshift(fft2(I_dpc_x)) + conj(H_y) .* fftshift(fft2(I_dpc_y));
denominator = abs(H_x).^2 + abs(H_y).^2 + alpha;

phase_recovered = real(ifft2(ifftshift(numerator ./ denominator)));

% 8. 结果可视化
figure('Color', 'w', 'Position', [100, 100, 1200, 400]);
subplot(1,4,1); imshow(I_dpc_y, []); title('DPC 采集图像 (垂直轴)');
subplot(1,4,2); imshow(abs(H_y), []); title('相位传递函数 |H_y|');
subplot(1,4,3); imshow(phase_true, []); title('原始真实相位');
subplot(1,4,4); imshow(phase_recovered, []); title('重建相位 (Tikhonov)');

colormap jet;