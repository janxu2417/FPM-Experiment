% =========================================================================
% 傅里叶叠层显微术(FPM) 双曝光视频处理脚本 (含冗余帧与无效帧跳过机制)
% 目标：基于索引提取长/短曝光视频 -> 中心300x300红通道 -> SURF全局对齐 -> 线性HDR融合
% =========================================================================

clear; clc; close all;

%% 1. 参数设置
videoFile_long = 'Raw_input\录像_02_long.avi';    % 长曝光视频路径
videoFile_short = 'Raw_input\录像_02_short.avi'; % 短曝光视频路径
exp_long = 1.0;                          % 长曝光时间(s)
exp_short = 0.1;                        % 短曝光时间(s)
exp_ratio = exp_long / exp_short;        % 曝光比例 (10倍)

crop_size = 300;                         % 裁剪大小 300x300
sat_thresh = 245;                        % 长曝光饱和阈值

% --- 新增：帧索引与步长控制 ---
start_idx_long = 7;   % 长曝光的第一个有效帧索引 (请根据实际视频修改)
start_idx_short = 48;  % 短曝光的第一个有效帧索引 (请根据实际视频修改)
frame_step_long = 2;       % 步长：2代表"取一帧，跳过一帧"(即每隔一帧为无效图)
frame_step_short = 2;

dark_noise_level = 0;    % 相机本底底噪平均值 (8位图通常在 1~3 左右，需根据你的相机暗电流测试测定)
valid_signal_thresh = 5; % 有效信号阈值。如果扣除底噪后，全图最大像素值仍低于此值，视为纯噪声帧

%% 2. 初始化视频读取器
v_long = VideoReader(videoFile_long);
v_short = VideoReader(videoFile_short);

% 计算裁剪区域的起始坐标 (居中裁剪)
height = v_long.Height;
width = v_long.Width;
r_start = floor((height - crop_size) / 2) + 1;
c_start = floor((width - crop_size) / 2) + 1;
rect_crop = [c_start, r_start, crop_size-1, crop_size-1];

% 预分配 3D 矩阵
imlow_HDR = zeros(crop_size, crop_size, 100, 'double');

% 直接利用指定的有效起始帧进行特征匹配
frame1_long = read(v_long, start_idx_long);
frame1_short = read(v_short, start_idx_short);

% 提取红通道并裁剪
img1_long_R = imcrop(frame1_long(:,:,1), rect_crop);
img1_short_R = imcrop(frame1_short(:,:,1), rect_crop);

figure;
set(gcf,'outerposition',get(0,'ScreenSize'))
subplot(121);imshow(double(img1_long_R),[]);title(['raw image ' num2str(1)]);
subplot(122);imshow(double(img1_short_R),[]);title(['raw image ' num2str(2)]);
%{
    3. 处理第一个有效帧：计算全局对齐变换矩阵 (tform)
% 方法 A：相位相关法（强烈推荐，专门针对纯平移位移，速度极快且受曝光差异影响小）
try
    tform = imregcorr(img1_short_R, img1_long_R, 'translation');
catch
    % 方法 B：基于强度的迭代配准（如果 imregcorr 效果不佳或需要考虑微小旋转，可作为备用）
    [optimizer, metric] = imregconfig('multimodal');
    % 'rigid' 表示只允许平移和旋转
    tform = imregtform(img1_short_R, img1_long_R, 'rigid', optimizer, metric);
end
% ================= 修改部分结束 =================

fprintf('全局对齐矩阵计算完成，共需处理 %d 个有效帧...\n', numValidFrames);
%}

%% 4. 循环读取有效帧、对齐与线性 HDR 融合
for i = 1:100
    % 计算当前有效帧在原视频中的真实索引
    curr_idx_long = start_idx_long + (i - 1) * frame_step_long;
    curr_idx_short = start_idx_short + (i - 1) * frame_step_short;

    % 1. 定点读取双边帧
    f_long = read(v_long, curr_idx_long);
    % 2. 提取红通道 (R通道) 并裁剪
    R_long = imcrop(f_long(:,:,1), rect_crop);
    % 转为 double
    img_L = double(R_long);

    if i <= 50
        f_short = read(v_short, curr_idx_short);
        R_short = imcrop(f_short(:,:,1), rect_crop);
        img_S = double(R_short);
    end

    if i > 50 || max(img_L(:)) <= sat_thresh
        hdr_img = img_L;
    else
        fprintf('  -> 提示: 第 %d 帧信号较强，采用短曝光。\n', i);
        hdr_img = img_S;
    end
    % 3. 空间对齐 (将短曝光配准到长曝光)
    % R_short_aligned = imwarp(R_short, tform, 'OutputView', imref2d(size(R_long)));

    %{
    % 4. FPM 专用线性 HDR 融合
    if i <= 50
        hdr_img = img_L;
        saturated_mask = (img_L > sat_thresh); % 寻找长曝光过曝区
        hdr_img(saturated_mask) = img_S(saturated_mask) * exp_ratio; % 短曝光线性替换
    else
        hdr_img = img_L;
    end
    %}
    % ---------------- 新增：低照度与噪声处理逻辑 ----------------
    % 1. 扣除相机固定本底噪声
    hdr_img = hdr_img - dark_noise_level;
    hdr_img(hdr_img < 0) = 0; % 物理光强不能为负，截断为0

    % 2. 极低亮度帧判定 (滤除纯噪声帧)
    % 如果整张图最亮的像素都低于有效信号阈值，说明没有捕捉到任何有效的散射光
    if max(hdr_img(:)) < valid_signal_thresh
        hdr_img(:) = 0; % 全图置零，告诉FPM算法该角度无高频信息更新
        fprintf('  -> 提示: 第 %d 帧信号极弱，已置零处理。\n', i);
    else
        % (可选) 针对极暗但仍有微弱信号的帧，可以使用中值滤波去除相机的"热像素(Hot Pixels)"
        % 因为热像素通常是孤立的极亮点，会严重干扰FPM相位恢复
        % --- 智能热像素剔除 (仅针对异常孤立点) ---
        % 计算整张图的中值滤波结果
        median_filtered = medfilt2(hdr_img, [3 3], "symmetric");

        % 如果某个像素比它周围环境的中值高出设定阈值（比如高出 50），才认为它是热像素
        noise_std = 80;
        if i > 25
            noise_std = 50;
        end
        hot_pixel_mask = (hdr_img - median_filtered) > noise_std;
        num_hot = nnz(hot_pixel_mask);
        if num_hot > 0
            hdr_img(hot_pixel_mask) = median_filtered(hot_pixel_mask);
            %fprintf('第 %d 帧共检测到 %d 个热像素，noise_std = %d\n', i, num_hot, noise_std);
        end
        % 只把判定为热像素的地方，替换为周围的中值，其他真实信号完全不碰
        % ----------------------------------------
    end
    % -----------------------------------------------------------
    %{
    if mod(i, 10) == 0 && i <= 50
        figure;
        set(gcf,'outerposition',get(0,'ScreenSize'))
        subplot(131);imshow(double(img_L),[]);title(['raw image ' num2str(1)]);
        subplot(132);imshow(double(img_S),[]);title(['raw image ' num2str(2)]);
        subplot(133);imshow(double(hdr_img),[]);title(['raw image ' num2str(3)]);
    end
    %}
    % 存入 3D 矩阵
    imlow_HDR(:,:,i) = sqrt(hdr_img);

    if mod(i, 20) == 0
        fprintf('已处理: %d / %d 帧\n', i, 100);
    end
end

%% 5. 保存为 .mat 文件
aberration = 0;
theta = 1.0;
xint = 0;
yint = 0;
z = 2.5000e-05;
wlength = 6.28e-07;
disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, num: ',num2str(size(imlow_HDR,3))]);
imshow(imlow_HDR(:,:,1),[]);
out_dir = 'Data\实验数据';
out_name = 'FPM_RawData.mat';
save('FPM_RawData.mat','aberration', 'theta', 'xint', 'yint', 'z', 'wlength', 'imlow_HDR');
disp('数据处理完成，已保存为 FPM_RawData.mat');

%}

%load("FPM_RawData.mat")
%% 5. 可视化展示 FPM_data 矩阵
k_show = 9;
rows = ceil(sqrt(k_show));
cols = ceil(k_show / rows);

for j = 0:25:75 % 挑选几个跨度展示即可
    figure('Color', 'w');
    for idx = 1:k_show
        if j+idx <= size(imlow_HDR,3)
            subplot(rows, cols, idx);
            img_show = imlow_HDR(:, :, j+idx);
            imagesc(img_show);
            colormap gray; axis image off;
            title(sprintf('Frame %d\nMax Val: %.1f', j+idx, max(img_show(:))), 'FontSize', 10);
        end
    end
    sgtitle(sprintf('FPM 数据预览 (第 %d 至 %d 帧)', j+1, j+k_show), 'FontSize', 14, 'FontWeight', 'bold');
end
%}