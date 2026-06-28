close all; clear; clc;

%% read_transfer_0619
% 0619 实验数据的一键预处理入口：
% 1. 选择一个明确的数据 preset
% 2. 统一处理分辨率、曝光表、ROI、输出命名与保存
% 3. 默认生成当前 preset 对应的 .mat，避免继续混用旧中间结果
% 4. 保留调试开关，便于人工核对 ROI、单帧内容和频谱
%
% cfg 是本脚本的单一配置入口。
% 输入目录、曝光表、ROI、输出命名等参数都通过 cfg 提供，
% 共享配置来源于：
%   Code/function/fpm_0619_shared_config.m
%
% MATLAB 结构体语法示例：
%   cfg = struct('field1', value1, 'field2', value2, ...)
% 访问字段时使用点号：
%   cfg.input_subdir
%   cfg.preproc.exposure_ms
%
% 本脚本中 cfg 的生成方式：
%   cfg = fpm_0619_shared_config("preproc", active_preset);
% 如果以后新增一套数据，优先在共享配置文件中新增 preset。

%% ===== 1. User entry =====
% 默认入口：当前实际使用的 0619/z=_0.050_7
active_preset = "0619_zneg0050_7_main";

% 共享配置函数位于 Code/function 下。
% 为了保证从不同工作目录运行时都能找到它，这里先把 Code 整体加入路径。
code_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(code_root));

% 如需切换到别的数据，只改这里即可：
% active_preset = "0619_zneg0050_up";
% active_preset = "0619_zneg0050_high2nd_multi";

% cfg 包含本次预处理所需的全部配置，例如：
%   cfg.input_subdir        : 原始图像所在子目录
%   cfg.output_tag          : 输出 .mat 的命名标签
%   cfg.preproc.exposure_ms : 每张图对应的曝光时间表
%   cfg.preproc.base_roi    : 基础 ROI 的几何定义
%   cfg.roi_configs         : 在基础 ROI 上附加的多 ROI 偏移配置
cfg = fpm_0619_shared_config("preproc", active_preset);

%% ===== 2. Optical / sampling parameters =====
spsize = cfg.optics.spsize;
m = cfg.optics.m;
n = cfg.optics.n;
dark_noise_level = cfg.optics.dark_noise_level;
wlength = cfg.optics.wlength;
NA = cfg.optics.NA_obj;
z = cfg.optics.z;

%% ===== 3. Debug switches =====
set(0, 'DefaultFigureVisible', 'on');

% Debug 开关说明
% 1) show_crop_overlay
%    作用：在第一张原始图上叠加显示 base ROI 和实际保存的 ROI。
%    检查对象：ROI 的几何位置是否正确。
%    适用场景：换新数据、修改 m/n、修改 base_roi、修改 row_shift/col_shift 时优先打开。
%
% 2) show_selected_amplitude
%    作用：显示选定 LED 编号对应的裁剪后振幅图 amp。
%    检查对象：ROI 裁剪后的图像内容是否合理，是否真的截到目标样品区域。
%    适用场景：怀疑 ROI 虽然位置大致正确，但框内内容亮度异常、边缘异常时。
%
% 3) show_selected_spectrum
%    作用：显示选定 LED 编号对应 ROI 的对数频谱。
%    检查对象：裁剪后数据的频谱分布是否异常。
%    适用场景：怀疑某个 ROI 输入质量差，或想比较不同 ROI 的频谱差异时。
%
% 4) selected_debug_led
%    作用：指定要重点检查的 LED 编号，只对 amplitude / spectrum 调试图生效。
%    常见做法：选中心附近、边缘附近和一张代表性图，例如 [1, 13, 25]。
%
% 5) pause_after_first_frame
%    作用：处理完第一张图后进入 keyboard 暂停，手动检查工作区变量。
%    适用场景：需要逐步确认 base_rect、roi_rects、I_crop、amp 等中间量时。
%
% 推荐调试顺序：
% 先开 show_crop_overlay 确认“框是否放对”；
% 再开 show_selected_amplitude 看“框里的内容是否对”；
% 最后如有需要再看 show_selected_spectrum。
show_crop_overlay = true;
show_selected_amplitude = true;
show_selected_spectrum = false;
selected_debug_led = [5, 10, 25];

% 如需逐段停下检查，把 false 改成 true。
pause_after_first_frame = false;

%% ===== 4. Resolve paths and validate configuration =====
% 下面这些语句演示了 cfg 的典型用法：
% 1. 从 cfg 里取字段
% 2. 用这些字段拼接输入/输出路径
% 3. 在正式读图前做一致性检查
num_img = cfg.preproc.num_img;
batch_dir = fullfile(code_root, cfg.raw_root, cfg.batch_folder);
input_dir = fullfile(batch_dir, cfg.input_subdir);
output_dir = fullfile(batch_dir, cfg.preproc.output_subdir);

if ~isfolder(input_dir)
    error('Input folder does not exist: %s', input_dir);
end

if ~isfolder(output_dir)
    mkdir(output_dir);
end

if numel(cfg.preproc.exposure_ms) ~= num_img
    error('exposure_ms length (%d) must match num_img (%d).', numel(cfg.preproc.exposure_ms), num_img);
end

num_roi = numel(cfg.roi_configs);
if num_roi < 1
    error('At least one ROI config is required.');
end

fprintf('\n=== 0619 preprocessing preset ===\n');
fprintf('preset       : %s\n', cfg.preset_name);
fprintf('input folder : %s\n', input_dir);
fprintf('output dir   : %s\n', output_dir);
fprintf('output tag   : %s\n', cfg.output_tag);
fprintf('ROI count    : %d\n\n', num_roi);

%% ===== 5. Prepare output storage =====
im_high_HDR_list = cell(1, num_roi);
for roi_idx = 1:num_roi
    im_high_HDR_list{roi_idx} = zeros(m, n, num_img, 'double');
end

%% ===== 6. Read images, apply exposure correction, crop ROIs =====
first_frame_info = struct();

for img_idx = 1:num_img
    data_path = fullfile(input_dir, "图像_" + string(cfg.preproc.file_indices(img_idx)) + ".tif");
    if ~isfile(data_path)
        error('Missing raw image: %s', data_path);
    end

    I_hr = imread(data_path);
    I_hr = double(I_hr(:, :, 1));

    fprintf('Image %2d/%2d : %s\n', img_idx, num_img, data_path);
    fprintf('    size = %d x %d, max = %.0f, exposure = %.3f ms\n', ...
        size(I_hr, 1), size(I_hr, 2), max(I_hr(:)), cfg.preproc.exposure_ms(img_idx));

    I_hr = I_hr - dark_noise_level;
    I_hr(I_hr < 0) = 0;

    [base_rect, roi_rects] = local_build_roi_rects(I_hr, m, n, cfg.preproc.base_roi, cfg.roi_configs);

    if img_idx == 1
        first_frame_info.base_rect = base_rect;
        first_frame_info.roi_rects = roi_rects;
        first_frame_info.image_size = size(I_hr);
    end

    for roi_idx = 1:num_roi
        rect = roi_rects(roi_idx);
        I_crop = I_hr(rect.row_start:rect.row_end, rect.col_start:rect.col_end);
        I_crop = I_crop ./ cfg.preproc.exposure_ms(img_idx);

        amp = sqrt(I_crop);
        im_high_HDR_list{roi_idx}(:, :, img_idx) = amp;

        if show_selected_amplitude && ismember(img_idx, selected_debug_led)
            figure('Name', sprintf('Amplitude LED %d / ROI %s', img_idx, rect.ptr), 'Color', 'w');
            imshow(amp, []);
            title(sprintf('LED %d, ROI %s amplitude', img_idx, rect.ptr), 'Interpreter', 'none');
        end

        if show_selected_spectrum && ismember(img_idx, selected_debug_led)
            himFT = fftshift(fft2(amp));
            crop_margin = 50;
            row_range = (1 + crop_margin):(size(himFT, 1) - crop_margin);
            col_range = (1 + crop_margin):(size(himFT, 2) - crop_margin);
            figure('Name', sprintf('Spectrum LED %d / ROI %s', img_idx, rect.ptr), 'Color', 'w');
            imshow(log(1 + abs(himFT(row_range, col_range))), []);
            title(sprintf('LED %d, ROI %s log spectrum', img_idx, rect.ptr), 'Interpreter', 'none');
        end
    end

    if pause_after_first_frame && img_idx == 1
        fprintf('Paused after first frame. Inspect variables, then type return.\n');
        keyboard;
    end
end

%% ===== 7. Optional ROI overlay debug =====
if show_crop_overlay
    data_path = fullfile(input_dir, "图像_" + string(cfg.preproc.file_indices(1)) + ".tif");
    I_preview = imread(data_path);
    I_preview = double(I_preview(:, :, 1));

    figure('Name', 'ROI overlay', 'Color', 'w');
    imshow(I_preview, []);
    hold on;

    rectangle('Position', local_rect_to_position(first_frame_info.base_rect), ...
        'EdgeColor', [0.2, 0.6, 1.0], 'LineWidth', 1.2, 'LineStyle', '--');

    color_map = lines(num_roi);
    for roi_idx = 1:num_roi
        rect = first_frame_info.roi_rects(roi_idx);
        rectangle('Position', local_rect_to_position(rect), ...
            'EdgeColor', color_map(roi_idx, :), 'LineWidth', 1.5);
        text(rect.col_start, rect.row_start - 20, rect.ptr, ...
            'Color', color_map(roi_idx, :), 'FontSize', 10, 'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    title(sprintf('ROI overlay: %s', cfg.preset_name), 'Interpreter', 'none');
    hold off;
end

%% ===== 8. Save outputs with explicit names =====
% manifest 统一使用共享 schema。
% 预处理和重建都使用同一套顶层字段设计，便于后续直接比较来源。
manifest_extra = struct();
manifest_extra.source.file_indices = cfg.preproc.file_indices(:);
manifest_extra.source.exposure_ms = cfg.preproc.exposure_ms(:);
manifest_extra.source.base_roi = cfg.preproc.base_roi;
manifest_extra.source.roi_configs = cfg.roi_configs;
manifest_extra.outputs.output_tag = cfg.output_tag;

save_manifest = fpm_0619_build_manifest(cfg, "preproc", manifest_extra);

for roi_idx = 1:num_roi
    ptr = cfg.roi_configs(roi_idx).ptr;
    roi_row_shift = cfg.roi_configs(roi_idx).row_shift;
    roi_col_shift = cfg.roi_configs(roi_idx).col_shift;
    im_high_HDR = im_high_HDR_list{roi_idx};

    out_name = cfg.output_names{roi_idx};
    out_path = fullfile(output_dir, out_name);

    save_manifest.outputs.mat_file = out_name;
    save_manifest.outputs.ptr = ptr;
    save_manifest.outputs.roi_row_shift = roi_row_shift;
    save_manifest.outputs.roi_col_shift = roi_col_shift;

    save(out_path, ...
        'z', 'NA', 'wlength', 'spsize', 'im_high_HDR', ...
        'ptr', 'roi_row_shift', 'roi_col_shift', 'save_manifest');

    fprintf('Saved: %s\n', out_path);
end

fprintf('\nPreprocessing finished for preset %s.\n', cfg.preset_name);

function [base_rect, roi_rects] = local_build_roi_rects(I_hr, m, n, base_roi, roi_configs)
% 输入：
%   I_hr        : 当前原始图像
%   m, n        : 目标裁剪尺寸
%   base_roi    : 基础 ROI 的相对中心定义
%   roi_configs : 多个 ROI 的平移配置
%
% 输出：
%   base_rect   : 基础 ROI 的实际像素坐标
%   roi_rects   : 每个实际 ROI 的坐标结构体数组
center_y = round(size(I_hr, 1) / 2);
center_x = round(size(I_hr, 2) / 2);

base_rect = struct();
base_rect.row_start = round(center_y + m * base_roi.row_start_factor + 1);
base_rect.row_end = round(center_y + m * base_roi.row_end_factor);
base_rect.col_start = round(center_x + n * base_roi.col_start_factor + 1);
base_rect.col_end = round(center_x + n * base_roi.col_end_factor);
base_rect.ptr = "base";

if base_rect.row_end - base_rect.row_start + 1 ~= m
    error('Base ROI height is not equal to m.');
end
if base_rect.col_end - base_rect.col_start + 1 ~= n
    error('Base ROI width is not equal to n.');
end

num_roi = numel(roi_configs);

% repmat(base_rect, num_roi, 1) 的意思是：
% 先复制出 num_roi 份同样字段的结构体，占好结构体数组的位置，
% 后面再逐个把每个 ROI 的坐标改掉。
roi_rects = repmat(base_rect, num_roi, 1);

for roi_idx = 1:num_roi
    rect = base_rect;
    rect.row_start = rect.row_start + roi_configs(roi_idx).row_shift;
    rect.row_end = rect.row_end + roi_configs(roi_idx).row_shift;
    rect.col_start = rect.col_start + roi_configs(roi_idx).col_shift;
    rect.col_end = rect.col_end + roi_configs(roi_idx).col_shift;
    rect.ptr = roi_configs(roi_idx).ptr;

    if rect.row_start < 1 || rect.row_end > size(I_hr, 1) || ...
            rect.col_start < 1 || rect.col_end > size(I_hr, 2)
        error('ROI %s exceeds image bounds: rows [%d, %d], cols [%d, %d].', ...
            rect.ptr, rect.row_start, rect.row_end, rect.col_start, rect.col_end);
    end

    roi_rects(roi_idx) = rect;
end
end

function pos = local_rect_to_position(rect)
pos = [rect.col_start, rect.row_start, ...
       rect.col_end - rect.col_start + 1, ...
       rect.row_end - rect.row_start + 1];
end
