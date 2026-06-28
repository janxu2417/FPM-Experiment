%% FPM 绝对理想电子束靶标生成器 (严格匹配真实靶标的等比缩放特征)
clear; close all; clc;

%% 1. 物理坐标系设置
psize = 0.138e-6; % 物理像素尺寸 改为 69 nm (极高分辨率网格)
% 视场大小：宽 120 um, 高 200 um (根据等比排布动态调整了视野边界)
x = -50e-6 : psize : 70e-6;
y = -80e-6 : psize : 120e-6;
[X, Y] = meshgrid(x, y);

O_ideal = zeros(size(X)); % 初始化全黑背景

%% 2. 定义目标频率与空间布局
% 目标频率 (线对/mm)
f_list = [330, 300, 270, 240, 220];
labels = {'330', '300', '270', '240', '220'};

% 初始顶部中心位置 (最高频 330 线对的起始位置)
current_y_center = 100e-6;

%% 3. 循环生成每一组线对与数字 (等比缩放排布)
for i = 1:length(f_list)
    f_spatial = f_list(i) * 1000; % 转换为 线对/米
    T = 1 / f_spatial;            % 周期 (米)
    line_width = T / 2;           % 50% 占空比，线宽=半周期

    y_c = current_y_center;

    % --- (A) 左侧 5 条竖线 ---
    % 真实靶标中，竖线通常比较长。设高度为 8 个周期 (上下各 4T)
    x_start_v = -15e-6;
    x_end_v   = x_start_v + 4.5 * T;
    y_start_v = y_c - 3.5 * T;
    y_end_v   = y_c + 3.5 * T;

    mask_v = (X >= x_start_v & X <= x_end_v) & (Y >= y_start_v & Y <= y_end_v);
    O_ideal(mask_v) = double(mod(X(mask_v) - x_start_v, T) < line_width);

    % --- (B) 右侧 5 条横线 ---
    % 横线位于竖线右侧，底部与竖线底部平齐（或略高），高度为 4.5 个周期
    x_start_h = x_end_v + T; % 距离竖线右侧留出 1 个周期的空白
    x_end_h   = x_start_h + 7 * T; % 横线长度约为 7 个周期
    y_start_h = y_start_v; % 底部与竖线对齐
    y_end_h   = y_start_h + 4.5 * T;

    mask_h = (X >= x_start_h & X <= x_end_h) & (Y >= y_start_h & Y <= y_end_h);
    O_ideal(mask_h) = double(mod(Y(mask_h) - y_start_h, T) < line_width);

    % --- (C) 生成顶部数字 ---
    % 数字位于横线正上方，大小和间距严格与当前周期 T 成正比
    text_x = x_start_h;
    text_y = y_end_h + 0.5 * T; % 距离横线顶部留出 0.5T 空白
    O_ideal = add_digital_text(O_ideal, X, Y, labels{i}, text_x, text_y, 0.4 * T);

    % --- (D) 动态计算下一组的 y_center ---
    % 核心修改：间距随着周期变化。当前图案占用的总高度约 10T，加几倍 T 的留白
    if i < length(f_list)
        T_next = 1 / (f_list(i+1) * 1000);
        % 向下推移：包含当前图案的下半部分留白 + 下一个图案的上半部分空间
        current_y_center = current_y_center - 4.5 * T - 5.5 * T_next;
    end
end

%% 4. 后处理与显示
% 给予一个微弱的非零背景防止对数计算报错
O_ideal = max(O_ideal, 0.01);
O_ideal = flip(O_ideal, 1);
save('Raw_input\FP_ideal_data.mat', 'O_ideal');

figure('Name', '理想数学靶标', 'Position', [200, 100, 500, 800]);
imshow(O_ideal, [], 'XData', x*1e6, 'YData', y*1e6);
axis ij; axis on;
title('绝对理想 50% 占空比靶标 (严格物理等比排布)');
xlabel('X (\mu m)'); ylabel('Y (\mu m)');
colormap('gray');

%saveas(gcf, 'D:\1_徐前\物理\傅里叶叠层显微术\Code\模拟测试\ideal_target.png');

%% =================================================================
%% 局部辅助函数：用像素点阵绘制数字
function O = add_digital_text(O, X, Y, str, x0, y0, pixel_size)
    % 简易 5x3 像素点阵字库
    font.n0 = [1 1 1; 1 0 1; 1 0 1; 1 0 1; 1 1 1];
    font.n2 = [1 1 1; 0 0 1; 1 1 1; 1 0 0; 1 1 1];
    font.n3 = [1 1 1; 0 0 1; 1 1 1; 0 0 1; 1 1 1];
    font.n4 = [1 0 1; 1 0 1; 1 1 1; 0 0 1; 0 0 1];
    font.n7 = [1 1 1; 0 0 1; 0 1 0; 0 1 0; 0 1 0];

    curr_x = x0;
    for k = 1:length(str)
        char_k = str(k);
        field_name = ['n', char_k];

        if isfield(font, field_name)
            mat = font.(field_name);
            for r = 1:5
                for c = 1:3
                    if mat(r, c) == 1
                        % 坐标映射
                        px = curr_x + (c-1)*pixel_size;
                        py = y0 + (5-r)*pixel_size;
                        mask = (X >= px & X < px+pixel_size) & ...
                               (Y >= py & Y < py+pixel_size);
                        O(mask) = 1;
                    end
                end
            end
        end
        % 字符间距跟随像素大小等比缩放
        curr_x = curr_x + 4 * pixel_size;
    end
end