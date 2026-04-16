% =========================================================================
% 傅里叶叠层显微术(FPM) 三曝光视频处理脚本 
% 目标：长/中/短三曝光读取 -> 提取红通道 -> 三级HDR物理融合 -> 噪声剔除 -> 保存
% =========================================================================

clear; clc; close all;

%% 1. 参数设置
videoFile_long  = 'Raw_input\录像_03_long.avi';   % 长曝光视频
videoFile_mid   = 'Raw_input\录像_03_mid.avi';    % 中曝光视频 (新增)
videoFile_short = 'Raw_input\录像_03_short.avi';  % 短曝光视频

exp_long  = 1.0;    % 长曝光时间(s)
exp_mid   = 0.5;    % 中曝光时间(s)
exp_short = 0.1;   % 短曝光时间(s)

% 计算曝光倍率基准（全部向长曝光物理光强对齐）
ratio_ML = exp_long / exp_mid;    % 中曝光放大倍数 (如 10)
ratio_SL = exp_long / exp_short;  % 短曝光放大倍数 (如 100)

crop_size = 300;     % 裁剪大小 300x300
sat_thresh = 235;    % 饱和阈值 (接近255的值)

% --- 帧索引与步长控制 ---
start_idx_long  = 3;   
start_idx_mid   = 8;  % 中曝光的第一个有效帧索引 (需根据视频修改) 
start_idx_short = 15;  
frame_step_long  = 2;  
frame_step_mid   = 2;
frame_step_short = 5;

dark_noise_level = 0;    
valid_signal_thresh = 5; 

%% 2. 初始化视频读取器
v_long  = VideoReader(videoFile_long);
v_mid   = VideoReader(videoFile_mid);
v_short = VideoReader(videoFile_short);

% 计算裁剪区域的起始坐标 (居中裁剪)
height = v_long.Height;
width = v_long.Width;
r_start = floor((height - crop_size) / 2) + 1;
c_start = floor((width - crop_size) / 2) + 1;
rect_crop = [c_start, r_start, crop_size-1, crop_size-1]; 
%{
for i = 1:2:11
    curr_idx_mid   = start_idx_mid   + (i - 1) * frame_step_mid;
    f_mid = read(v_mid, curr_idx_mid);
    img_M = double(imcrop(f_mid(:,:,1), rect_crop));
    figure();imshow(img_M,[]);title(sprintf('Frame M %d\nMax Val: %.1f', i, max(img_M(:))), 'FontSize', 10);
end
%}

% 预分配 3D 矩阵
imlow_HDR = zeros(crop_size, crop_size, 100, 'double');

%% 3. 循环读取有效帧与三级 HDR 物理融合
for i = 1:100
    % 1. 计算真实索引
    curr_idx_long  = start_idx_long  + (i - 1) * frame_step_long;
    curr_idx_mid   = start_idx_mid   + (i - 1) * frame_step_mid;
    curr_idx_short = start_idx_short + (i - 1) * frame_step_short;
    
    % 2. 读取长曝光并处理
    f_long = read(v_long, curr_idx_long);
    img_L  = double(imcrop(f_long(:,:,1), rect_crop));

    % 默认 HDR 图像基底为长曝光
    hdr_img = img_L; 
    
    % 3. 三级 HDR 融合逻辑 (逐像素判定)
    % 步骤 A：寻找长曝光过曝的像素
    mask_L_sat = (img_L > sat_thresh);
    
    if any(mask_L_sat, 'all') % 仅在可能过曝的明场区域做融合
        % 读取中曝光
        f_mid = read(v_mid, curr_idx_mid);
        img_M = double(imcrop(f_mid(:,:,1), rect_crop));
        
        % 用中曝光(乘以倍率)替换长曝光的过曝区
        %hdr_img(mask_L_sat) = img_M(mask_L_sat) * ratio_ML;
        hdr_img = img_M * ratio_ML;
        
        % 步骤 B：检查中曝光是否也过曝了
        mask_M_sat = (img_M > sat_thresh);
        mask_need_short = mask_L_sat & mask_M_sat; % 必须是长、中都过曝的极亮区域
        
        if any(mask_need_short, 'all') && i <= 45
            % 读取短曝光
            f_short = read(v_short, curr_idx_short);
            img_S   = double(imcrop(f_short(:,:,1), rect_crop));
            
            % 用短曝光(乘以倍率)替换极端过曝区
            %hdr_img(mask_need_short) = img_S(mask_need_short) * ratio_SL;
            hdr_img = img_S * ratio_SL;
            fprintf('  -> 提示: 第 %d 帧信号较强，采用短曝光。\n', i);
        end
    end
    
    % ---------------- 低照度与噪声处理逻辑 ----------------
    % 扣除底噪
    hdr_img = hdr_img - dark_noise_level;
    hdr_img(hdr_img < 0) = 0; 

    % 极低亮度帧判定
    if max(hdr_img(:)) < valid_signal_thresh
        hdr_img(:) = 0; 
        fprintf('  -> 提示: 第 %d 帧信号极弱，已置零处理。\n', i);
    else
        % 智能热像素剔除
        median_filtered = medfilt2(hdr_img, [3 3], "symmetric");
        noise_std = 80;
        if i > 25
            noise_std = 50;
        end
        hot_pixel_mask = (hdr_img - median_filtered) > noise_std; 
        if any(hot_pixel_mask, 'all')
            hdr_img(hot_pixel_mask) = median_filtered(hot_pixel_mask);
        end
    end
    % -----------------------------------------------------------
    
    if i == 1 || i >=6 && i <=10
        figure();
        subplot(131);imshow(img_L,[]);title(sprintf('Frame L %d\nMax Val: %.1f', i, max(img_L(:))), 'FontSize', 10);
        subplot(132);imshow(img_M,[]);title(sprintf('Frame M %d\nMax Val: %.1f', i, max(img_M(:))), 'FontSize', 10);
        subplot(133);imshow(hdr_img,[]);title(sprintf('Frame HDR %d\nMax Val: %.1f', i, max(hdr_img(:))), 'FontSize', 10);
        sgtitle(sprintf('FPM 数据预览 (第 %d 帧)', i), 'FontSize', 14, 'FontWeight', 'bold');
    end
    

    % 取振幅存入 3D 矩阵 (注意：你的原始输入是强度，这里开根号转振幅非常正确)
    imlow_HDR(:,:,i) = sqrt(hdr_img);

    if mod(i, 20) == 0
        fprintf('已处理: %d / 100 帧\n', i);
    end
end

%% 4. 保存为 .mat 文件
aberration = 0;
theta = 1.0;
xint = 0;
yint = 0;
z = 2.5000e-05;
wlength = 6.28e-07;
disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, num: ',num2str(size(imlow_HDR,3))]);

out_name = 'FPM_RawData_3HDR.mat';
save(out_name,'aberration', 'theta', 'xint', 'yint', 'z', 'wlength', 'imlow_HDR');
disp('数据处理完成，已保存为 FPM_RawData_3HDR.mat');
%{
%% 5. 可视化展示 FPM_data 矩阵
k_show = 4; 
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