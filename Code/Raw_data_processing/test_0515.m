close all;clear;clc;

%% ========== 1. 参数设置 ==========
% 低分辨率像素尺寸 (m) (对应你的 0.15NA 实际实验)
spsize  = 0.69e-6;
% % 上采样比例（假设 0.15NA 到 0.8NA 大约有 4~5 倍分辨率提升）
% pratio = 4;
% % 高分辨率像素尺寸 (m)
% psize = spsize / pratio;

% 低分辨率图像尺寸（像素）
m1 = 600;  % 高度
n1 = 600;  % 宽度

m = m1 * 2;
n = n1 * 2;

dark_noise_level = 0.005;

% 光学参数
wlength = 6.28e-07; % 波长 (m)
NA  = 0.15;     % 模拟的低分辨率物镜数值孔径 (修改为你的实际物镜NA)
z   = 2e-6;        % 离焦量 (m) 模拟理想情况先设为0

%% ========== 2. 导入真实高分辨率图像并合成复振幅 ==========
% 读取你拍摄的 0.8NA TIFF 图像
folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Raw_input', '0427');
filename = sprintf('0.15NA白光分辨率8192.tif');
data_name = fullfile(folder, filename);
I_hr = imread(data_name);
% 【注】：这里我用伪代码代替，请换成你实际的文件名
I_hr = double(I_hr(:,:,1));
fprintf('图 %s 的尺寸为 %d %d\n', data_name, size(I_hr));
fprintf('这幅图的物理最大光强值是: %d\n', max(I_hr(:)));
figure;
imshow(I_hr, []);
title('MATLAB 自动对比度显示的真实数据');