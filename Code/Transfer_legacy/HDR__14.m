%% HDR合成函数：将两张不同曝光的图像合成为14位HDR图像
function im_HDR_14bit = synthesize_HDR_14bit(im_short, im_long, t_short, t_long)
    % 输入参数：
    %   im_short, im_long: 两张不同曝光的图像矩阵 (double类型)
    %   t_short, t_long: 对应的曝光时间 (秒)
    % 输出：
    %   im_HDR_14bit: 合成后的14位HDR图像 (uint16类型，范围0-16383)

    % 确保图像是double类型
    im_short = double(im_short);
    im_long = double(im_long);

    % 获取图像尺寸
    [height, width] = size(im_short);

    % 初始化HDR图像（浮点类型）
    im_HDR_float = zeros(height, width);

    % 定义饱和阈值（假设8位图像，饱和阈值为250，留一些余量）
    saturation_threshold = 250;

    % 定义信噪比权重（可以根据实际调整）
    noise_floor = 10;  % 低于此值认为信噪比太低

    % 方法1：像素级优化选择（推荐用于FPM）
    for y = 1:height
        for x = 1:width
            % 获取两个曝光的像素值
            p_short = im_short(y, x);
            p_long = im_long(y, x);

            % 情况1：长曝光未饱和且信噪比足够
            if p_long < saturation_threshold && p_long > noise_floor
                im_HDR_float(y, x) = p_long / t_long;

            % 情况2：长曝光饱和或信噪比低，使用短曝光
            else
                im_HDR_float(y, x) = p_short / t_short;
            end
        end
    end

    % 对HDR图像进行归一化和14位量化
    min_val = prctile(im_HDR_float(:), 1);  % 1%分位数
    max_val = prctile(im_HDR_float(:), 99); % 99%分位数

    % 线性映射到14位范围 [0, 16383]
    im_HDR_14bit = uint16((im_HDR_float - min_val) / (max_val - min_val) * 16383);
end

%% 主程序：处理两组曝光图像并合成HDR
clear; close all; clc;

% 设置曝光时间
t1 = 0.005;  % 短曝光
t3 = 0.5;    % 长曝光

% 初始化两个图像矩阵
imlow_HDR_short = zeros(300, 300, 100);
imlow_HDR_long = zeros(300, 300, 100);

%% 读取短曝光图像序列 (0.005s)
disp('读取短曝光图像序列 (0.005s)...');
end_j = 0;
for i = 7:2:205
    base_name = 'Data\实验数据\result_5\frame_';
    num_str = num2str(i);
    filename = [base_name, num_str, '.png'];

    if ~exist(filename, 'file')
        warning(['文件不存在: ', filename]);
        continue;
    end

    rgbImage = imread(filename);
    redChannel = rgbImage(:, :, 1);  % R通道

    end_j = end_j + 1;
    imlow_HDR_short(:,:,end_j) = double(redChannel(601:900, 801:1100));
end

%% 读取长曝光图像序列 (0.5s)
disp('读取长曝光图像序列 (0.5s)...');
end_j = 0;
base_name_long = 'Data\实验数据\result_7\frame_';  %% 修改路径为长曝光图像所在位置
for i = 7:2:205
    num_str = num2str(i);
    filename = [base_name_long, num_str, '.png'];

    if ~exist(filename, 'file')
        warning(['文件不存在: ', filename]);
        continue;
    end

    rgbImage = imread(filename);
    redChannel = rgbImage(:, :, 1);  % R通道

    end_j = end_j + 1;
    imlow_HDR_long(:,:,end_j) = double(redChannel(601:900, 801:1100));
end

%% 对每一帧进行HDR合成
disp('开始HDR合成...');
num_frames = size(imlow_HDR_short, 3);
im_HDR_14bit_all = zeros(300, 300, num_frames, 'uint16');

for frame_idx = 1:num_frames
    im_short = imlow_HDR_short(:,:,frame_idx);
    im_long = imlow_HDR_long(:,:,frame_idx);

    im_HDR_14bit = synthesize_HDR_14bit(im_short, im_long, t1, t3);

    im_HDR_14bit_all(:,:,frame_idx) = im_HDR_14bit;

    if mod(frame_idx, 20) == 0
        fprintf('已处理 %d/%d 帧\n', frame_idx, num_frames);
    end
end

%% 可视化结果
disp('显示合成结果...');
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
imagesc(imlow_HDR_short(:,:,1));
title(sprintf('短曝光 (%.3fs)', t1));
colormap gray; colorbar;
clim([0 255]);

subplot(1, 3, 2);
imagesc(imlow_HDR_long(:,:,1));
title(sprintf('长曝光 (%.1fs)', t3));
colormap gray; colorbar;
clim([0 255]);

subplot(1, 3, 3);
imagesc(double(im_HDR_14bit_all(:,:,1)));
title('合成14位HDR图像');
colormap gray; colorbar;
clim([0 16383]);

sgtitle('双曝光HDR合成结果');

%% 保存结果
disp('保存HDR合成结果...');
save('HDR_14bit_results.mat', 'im_HDR_14bit_all', 't1', 't3');

disp('HDR合成完成！');