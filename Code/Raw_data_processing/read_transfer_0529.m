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

m = m1;
n = n1;

dark_noise_level = 0;

% 光学参数
wlength = 6.28e-07; % 波长 (m)
NA  = 0.15;     % 模拟的低分辨率物镜数值孔径 (修改为你的实际物镜NA)
z   = 2e-6;        % 离焦量 (m) 模拟理想情况先设为0

%% ========== 2. 导入真实高分辨率图像并合成复振幅 ==========

num_img = 25;
im_high_HDR = zeros(m, n, num_img, 'double');
im_up_HDR = zeros(m, n, num_img, 'double');
% im_med_HDR = zeros(m1, n1, num_img, 'double');

% 放在你的 MATLAB 脚本最前面
set(0, 'DefaultFigureVisible', 'on');

% exposure_ms = 42.307 * ones(1, num_img);
% exposure_ms([5,9]) = 50.325;
% exposure_ms(10:12) = 100.65;
% exposure_ms([13,17,21,25]) = 201.301;
% exposure_ms([14,15,16,22,24]) = 142.34;
% exposure_ms([18,19,20,23]) = 119.691;
exposure_ms = 100.65 * ones(1, num_img);
exposure_ms(1:9) = 42.307;
exposure_ms([13,17,21,25]) = 201.301;

str_z = "z=_0.050_high";
ptr = "";
batch_folder = "0619";
root = 'D:\1_徐前\物理\傅里叶叠层显微术\Code';
raw_root = "Raw_input";

for i = 1:num_img

    data_name = batch_folder + "\" + str_z + "\图像_" + int2str(i);
    data_dir = fullfile(raw_root, data_name + ".tif");
    % 读取你拍摄的 0.8NA TIFF 图像
    I_hr = imread(data_dir);
    I_hr = double(I_hr(:,:,1));
    fprintf('图 %s 的尺寸为 %d %d\n', data_name, size(I_hr));
    fprintf('这幅图的物理最大光强值是: %d\n', max(I_hr(:)));
    % figure;
    % imshow(I_hr, []);
    % title('MATLAB 自动对比度显示的真实数据');

    % 为了演示代码能跑通，如果找不到图，生成一个测试矩阵（你用时将其注释掉）
    %I_hr = double(imresize(rgb2gray(imread('cameraman.tif')), [4096, 3000]));

    % 0. 扣除相机固定本底噪声
    I_hr = I_hr - dark_noise_level;
    I_hr(I_hr < 0) = 0;

    % 1. 截取中心区域 (m x n) 对应 FPM 的视场
    center_y = round(size(I_hr, 1) / 2);
    center_x = round(size(I_hr, 2) / 2);
    row_start = center_y - m * 12 / 6 + 1;
    row_end = center_y + m * -6 / 6;
    col_start = center_x - n * 2 / 3 + 1;
    col_end = center_x + n * 1 / 3;
    I_crop = I_hr(row_start:row_end, col_start:col_end);

    % 2. 除以曝光时间
    I_crop = I_crop ./ exposure_ms(i);

    % 3. 从强度计算振幅 (Amplitude = sqrt(Intensity))
    amp = sqrt(I_crop);

    % 4. 物理模拟：合成"真实"的相位
    % 在生物细胞等透明样品中，相位延迟通常正比于样本厚度或吸收（振幅）
    % 这里我们假设一个随振幅变化的弱相位物体（最大相位差 pi）
    % phase = (1 - amp) * pi;
    phase = 0;
    % % 合成高分辨率基准复振幅 O
    O = amp .* exp(1j * phase);

    im_high_HDR(:,:,i) = amp;
    % im_med_HDR(:,:,i) = imresize(amp, 0.5, "bicubic", "Antialiasing", true);

    % if (i == 1 || i == 5 || i == 10)
    %     figure('Name', '基准真值 Ground Truth (Ideal)');
    %     imshow(abs(O), []); title('HR 振幅 (真实数据)');
    % end

    % subplot(122); imshow(im_med_HDR(:,:,i), []); title('HR 振幅 (压缩一半)');

    %     % subplot(122); imshow(angle(O), []); title('HR 相位 (合成物理模型)');
    % filename = sprintf('LED_%d_高分辨率原图.png', i);
    % saveas(gcf, fullfile(folder, filename));

    % himFT = fftshift(fft2(I_hr));
    % figure('Name', '高分辨率频谱分析');
    % % 使用 log(1 + abs(...)) 进行动态范围压缩
    % % subplot(121);
    % imshow(log(1 + abs(himFT(50:end-50, 50:end-50))), []);
    % title('Log Amplitude (对数振幅谱)');
    %
    % folder = fullfile(root, raw_root, batch_folder, '原图频谱', str_z + ptr);
    % if ~exist(folder, 'dir')
    %     mkdir(folder);
    % end
    %
    % filename = sprintf('LED_%d_高分辨率频谱.png', i);
    % saveas(gcf, fullfile(folder, filename));

    % subplot(122);
    % imshow(angle(himFT(50:end-50, 50:end-50)), []);
    % title('Phase (相位谱)');
end
out_name = "FPM_RawData_z " + str_z + ptr +".mat";
folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Raw_input', batch_folder);
save(fullfile(folder, out_name), 'z', 'NA', 'wlength', 'spsize', "im_high_HDR");