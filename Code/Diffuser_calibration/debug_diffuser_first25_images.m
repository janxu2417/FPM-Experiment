close all; clear; clc;
set(0, 'DefaultFigureVisible', 'off');

%% ===== 1. Basic settings =====
% 用于调试前 25 个 LED 的原始图像。
% 这个脚本不做 ROI 拟合，只做图像级检查，重点回答两个问题：
% 1. 2~5、21~25 号 LED 的异常是否直接体现在原图上？
% 2. 异常是按采集顺序出现，还是按 LED 空间位置成簇出现？

script_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(script_dir);
addpath(genpath(code_dir));

batch_folder = "0605";
image_subfolder = "毛玻璃01";
file_index_list = 1:25;
arraysize = 5;

num_img = numel(file_index_list);
if num_img ~= arraysize^2
    error('numel(file_index_list) must equal arraysize^2.');
end

%% ===== 2. Exposure settings =====
exposure_ms = 50.325 * ones(1, num_img);
exposure_ms([13, 17, 21, 22, 25]) = 59.845;
dark_noise_level = 0;

%% ===== 3. Display and debug settings =====
center_patch_size = 201;
save_figures = true;
mark_idx = [2, 3, 4, 5, 21, 22, 23, 24, 25];

%% ===== 4. Resolve paths =====
raw_root = fullfile(code_dir, 'Raw_input', batch_folder, image_subfolder);
if ~isfolder(raw_root)
    error('Input folder does not exist: %s', raw_root);
end

result_dir = fullfile(code_dir, 'Results', 'Calibration', batch_folder, 'Debug');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

image_files = strings(num_img, 1);
for i = 1:num_img
    image_files(i) = fullfile(raw_root, "图像_" + string(file_index_list(i)) + ".tif");
    if ~isfile(image_files(i))
        error('Missing image for LED index %d: %s', i, image_files(i));
    end
end

%% ===== 5. Pre-read image size and LED geometry =====
I_first = local_read_first_channel(image_files(1));
image_size = size(I_first);
center_rect = local_build_center_rect(image_size, center_patch_size);

[xlocation, ylocation] = LED_location(0, 0, arraysize);
xlocation = xlocation(:);
ylocation = ylocation(:);
x_unique = sort(unique(xlocation), 'ascend');
y_unique = sort(unique(ylocation), 'descend');

%% ===== 6. Read all images and compute simple image-level metrics =====
raw_mean = zeros(num_img, 1);
corr_mean = zeros(num_img, 1);
raw_max = zeros(num_img, 1);
corr_max = zeros(num_img, 1);
center_mean_corr = zeros(num_img, 1);
raw_p995 = zeros(num_img, 1);
corr_p995 = zeros(num_img, 1);

raw_cache = cell(num_img, 1);
corr_cache = cell(num_img, 1);

for i = 1:num_img
    I_raw = local_read_first_channel(image_files(i));
    I_raw = double(I_raw);

    I_corr = I_raw;
    if dark_noise_level ~= 0
        I_corr = I_corr - dark_noise_level;
    end
    I_corr(I_corr < 0) = 0;
    I_corr = I_corr ./ exposure_ms(i);

    raw_cache{i} = I_raw;
    corr_cache{i} = I_corr;

    raw_mean(i) = mean(I_raw(:));
    corr_mean(i) = mean(I_corr(:));
    raw_max(i) = max(I_raw(:));
    corr_max(i) = max(I_corr(:));
    raw_p995(i) = prctile(I_raw(:), 99.5);
    corr_p995(i) = prctile(I_corr(:), 99.5);

    patch = I_corr(center_rect(3):center_rect(4), center_rect(1):center_rect(2));
    center_mean_corr(i) = mean(patch(:));
end

raw_clim = [0, max(raw_p995)];
corr_clim = [0, max(corr_p995)];

%% ===== 7. Figure 1: raw images in acquisition order =====
fig1 = figure('Color', 'w', 'Position', [50, 50, 1450, 950]);
tiledlayout(arraysize, arraysize, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle('Raw images in acquisition order');

for i = 1:num_img
    nexttile;
    imagesc(raw_cache{i}, raw_clim);
    axis image off;
    colormap gray;
    hold on;
    rectangle('Position', local_rect_to_position(center_rect), 'EdgeColor', 'r', 'LineWidth', 0.8);
    title(sprintf('LED %d | file %d\nexp %.3f ms | max %.0f', ...
        i, file_index_list(i), exposure_ms(i), raw_max(i)), 'FontSize', 8, 'Interpreter', 'none');
    if any(i == mark_idx)
        text(20, 50, 'marked', 'Color', 'y', 'FontSize', 9, 'FontWeight', 'bold');
    end
    hold off;
end

%% ===== 8. Figure 2: corrected images in acquisition order =====
fig2 = figure('Color', 'w', 'Position', [80, 80, 1450, 950]);
tiledlayout(arraysize, arraysize, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle('Exposure-corrected images in acquisition order');

for i = 1:num_img
    nexttile;
    imagesc(corr_cache{i}, corr_clim);
    axis image off;
    colormap gray;
    hold on;
    rectangle('Position', local_rect_to_position(center_rect), 'EdgeColor', 'r', 'LineWidth', 0.8);
    title(sprintf('LED %d | file %d\nmean %.3f | center %.3f', ...
        i, file_index_list(i), corr_mean(i), center_mean_corr(i)), 'FontSize', 8, 'Interpreter', 'none');
    if any(i == mark_idx)
        text(20, 50, 'marked', 'Color', 'y', 'FontSize', 9, 'FontWeight', 'bold');
    end
    hold off;
end

%% ===== 9. Figure 3: raw images arranged by LED spatial position =====
fig3 = figure('Color', 'w', 'Position', [110, 110, 1450, 950]);
tiledlayout(arraysize, arraysize, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle('Raw images arranged by LED spatial position');

for i = 1:num_img
    row = find(y_unique == ylocation(i), 1, 'first');
    col = find(x_unique == xlocation(i), 1, 'first');
    tile_idx = (row - 1) * arraysize + col;

    nexttile(tile_idx);
    imagesc(raw_cache{i}, raw_clim);
    axis image off;
    colormap gray;
    hold on;
    rectangle('Position', local_rect_to_position(center_rect), 'EdgeColor', 'r', 'LineWidth', 0.8);
    title(sprintf('LED %d\n(x=%d, y=%d)', i, xlocation(i), ylocation(i)), ...
        'FontSize', 8, 'Interpreter', 'none');
    if any(i == mark_idx)
        text(20, 50, 'marked', 'Color', 'y', 'FontSize', 9, 'FontWeight', 'bold');
    end
    hold off;
end

%% ===== 10. Figure 4: corrected images arranged by LED spatial position =====
fig4 = figure('Color', 'w', 'Position', [140, 140, 1450, 950]);
tiledlayout(arraysize, arraysize, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle('Exposure-corrected images arranged by LED spatial position');

for i = 1:num_img
    row = find(y_unique == ylocation(i), 1, 'first');
    col = find(x_unique == xlocation(i), 1, 'first');
    tile_idx = (row - 1) * arraysize + col;

    nexttile(tile_idx);
    imagesc(corr_cache{i}, corr_clim);
    axis image off;
    colormap gray;
    hold on;
    rectangle('Position', local_rect_to_position(center_rect), 'EdgeColor', 'r', 'LineWidth', 0.8);
    title(sprintf('LED %d\n(x=%d, y=%d)', i, xlocation(i), ylocation(i)), ...
        'FontSize', 8, 'Interpreter', 'none');
    if any(i == mark_idx)
        text(20, 50, 'marked', 'Color', 'y', 'FontSize', 9, 'FontWeight', 'bold');
    end
    hold off;
end

%% ===== 11. Figure 5: image-level summary curves =====
fig5 = figure('Color', 'w', 'Position', [170, 170, 1400, 420]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on; grid on; box on;
plot(1:num_img, raw_mean, 'k-o', 'LineWidth', 1.2);
plot(1:num_img, corr_mean, 'r-s', 'LineWidth', 1.2);
plot(mark_idx, corr_mean(mark_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('LED index');
ylabel('Whole-image mean');
title('Whole-image mean intensity');
legend('Raw mean', 'Corrected mean', 'Marked LEDs', 'Location', 'best');

nexttile;
hold on; grid on; box on;
plot(1:num_img, center_mean_corr, 'b-o', 'LineWidth', 1.2);
plot(mark_idx, center_mean_corr(mark_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('LED index');
ylabel('Center-patch mean (corrected)');
title('Center-patch mean intensity');
legend('Center patch', 'Marked LEDs', 'Location', 'best');

nexttile;
hold on; grid on; box on;
plot(1:num_img, raw_max, 'k-o', 'LineWidth', 1.2);
plot(1:num_img, corr_max, 'm-s', 'LineWidth', 1.2);
plot(mark_idx, corr_max(mark_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('LED index');
ylabel('Max intensity');
title('Max DN / corrected max');
legend('Raw max', 'Corrected max', 'Marked LEDs', 'Location', 'best');

%% ===== 12. Console summary =====
debug_table = table((1:num_img).', file_index_list(:), exposure_ms(:), ...
    xlocation(:), ylocation(:), raw_mean, corr_mean, center_mean_corr, raw_max, corr_max, ...
    'VariableNames', {'led_idx', 'file_idx', 'exposure_ms', 'xloc', 'yloc', ...
    'raw_mean', 'corr_mean', 'center_mean_corr', 'raw_max', 'corr_max'});

disp(debug_table);

fprintf('\nMarked LEDs: %s\n', mat2str(mark_idx));
fprintf('Raw display clim       = [%.3f, %.3f]\n', raw_clim(1), raw_clim(2));
fprintf('Corrected display clim = [%.3f, %.3f]\n', corr_clim(1), corr_clim(2));

%% ===== 13. Save figures =====
if save_figures
    file_tag = sprintf('%s_%s_%dx%d', char(batch_folder), local_safe_name(image_subfolder), arraysize, arraysize);
    saveas(fig1, fullfile(result_dir, "debug_raw_order_" + file_tag + ".png"));
    saveas(fig2, fullfile(result_dir, "debug_corrected_order_" + file_tag + ".png"));
    saveas(fig3, fullfile(result_dir, "debug_raw_spatial_" + file_tag + ".png"));
    saveas(fig4, fullfile(result_dir, "debug_corrected_spatial_" + file_tag + ".png"));
    saveas(fig5, fullfile(result_dir, "debug_summary_curves_" + file_tag + ".png"));
end

%% ===== Local functions =====
function I_ch = local_read_first_channel(image_file)
I_raw = imread(image_file);
if ndims(I_raw) >= 3
    I_ch = I_raw(:, :, 1);
else
    I_ch = I_raw;
end
end

function rect = local_build_center_rect(image_size, patch_size)
if isscalar(patch_size)
    patch_h = patch_size;
    patch_w = patch_size;
else
    patch_h = patch_size(1);
    patch_w = patch_size(2);
end

center_x = round((image_size(2) + 1) / 2);
center_y = round((image_size(1) + 1) / 2);
half_w = floor((patch_w - 1) / 2);
half_h = floor((patch_h - 1) / 2);

x1 = center_x - half_w;
x2 = x1 + patch_w - 1;
y1 = center_y - half_h;
y2 = y1 + patch_h - 1;

if x1 < 1 || y1 < 1 || x2 > image_size(2) || y2 > image_size(1)
    error('Center patch is out of image bounds.');
end

rect = [x1, x2, y1, y2];
end

function pos = local_rect_to_position(rect)
pos = [rect(1), rect(3), rect(2) - rect(1) + 1, rect(4) - rect(3) + 1];
end

function text_out = local_safe_name(text_in)
text_out = regexprep(char(text_in), '[^\w\+\-]', '_');
end
