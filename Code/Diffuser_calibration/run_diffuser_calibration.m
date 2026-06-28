close all; clear; clc;
set(0, 'DefaultFigureVisible', 'on');

%% ===== 1. Paths and data selection =====
% 主入口脚本：
% 用毛玻璃图像标定不同 LED 入射角下的相对照明强度。
% 核心步骤是：读图 -> 曝光校正 -> 生成 ROI -> 提取统计量 ->
% 选择 ROI 汇总方式 -> 按 radius / illumination NA 做平滑拟合。
%
% 如果结果异常，建议开启 debug_mode，在 MATLAB 编辑器中按 section
% 逐段运行，先检查 ROI 位置、图像亮区中心和各 ROI 的统计曲线。
% This script is the standard entry for diffuser-based LED intensity
% calibration. The core workflow is:
% 1) read images and correct exposure
% 2) estimate or set the ROI center
% 3) extract ROI statistics
% 4) compare several aggregation rules
% 5) fit smooth intensity curves versus radius / illumination NA
%
% When the result looks suspicious, enable debug_mode and run the script
% section by section in MATLAB Live Script / Editor. The script will then
% output ROI overlays, per-LED diagnostic plots and console tables.

script_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(script_dir);
addpath(genpath(code_dir));
% 将 Code 目录下的函数全部加入路径，便于调用 LED_location、k_vector 等函数

batch_folder = "0619";
image_subfolder = "光强修正";
file_index_list = 1:49;
arraysize = 7;
% file_index_list 是参与本次标定的图像编号列表
% arraysize 是 LED 子阵列边长，这里 5 表示 5x5

%% ===== 2. Exposure and geometry =====
num_img = numel(file_index_list);
% 1-9                       42.307
% 13，17，21，25  100.65
% 其余                      84.629
% exposure_ms = 50.325 * ones(1, num_img);
% exposure_ms([13, 17, 21, 22, 25]) = 59.845;
% if arraysize == 7
%     exposure_ms(26:35) = 100.65;
%     exposure_ms([31, 32, 36, 38:42, 45:48]) = 142.34;
%     exposure_ms([37, 43, 44, 49]) = 201.301;
% end

exposure_ms = zeros(1, num_img);

exposure_ms(1:9) = 42.307;
exposure_ms([10:12, 14:16, 18:20, 22:24]) = 84.629;
exposure_ms([13, 17, 21, 25]) = 100.65;
exposure_ms([26:30, 32:36, 38:42, 44:48]) = 148.645;
exposure_ms([31, 37, 43, 49]) = 210.218;

% 曝光时间数组，后续按 counts/ms 做归一化
% 若部分 LED 使用了不同曝光，这里逐项单独覆盖

if num_img ~= arraysize^2
    error('numel(file_index_list) must equal arraysize^2.');
end
if numel(exposure_ms) ~= num_img
    error('exposure_ms length must match numel(file_index_list).');
end

LEDp = 4;
H = 40.24;
nglass = 1.52;
t = 1;
theta = 1.0;
xstart = 0;
ystart = 0;
xint = 0;
yint = 0;
% 几何参数与 FPM 主恢复程序保持一致：
% LEDp 是 LED 间距，H 是 LED 到样品距离，
% nglass / t 是玻璃参数，theta 是 LED 阵列相对相机的旋转角

%% ===== 3. Preprocessing, ROI and debug =====
dark_noise_level = 0;
% 当前实验条件下暗噪声按 0 处理，但接口保留，便于以后切换为暗场扣除

% ROI sizes are intentionally configurable. If edge LEDs look abnormally
% strong, the first thing to check is whether the ROI is too large or not
% centered on the effective diffuser bright zone.
main_roi_size = [2001, 2001];
aux_roi_size = [151, 151];
aux_roi_offsets = [
    -320,   0;
     320,   0;
       0, -320;
       0,  320;
    -220, -220;
     220, -220;
    -220,  220;
     220,  220
];
exclude_idx = [];
% exclude_idx 用于指定“不参与拟合”的 LED 编号
% 常见用途是临时排除明显坏点或边界异常点

trim_fraction = 0;
max_dn_override = 900;
hist_tail_bins_pct = 1;
qc_thresholds = struct( ...
    'sat_ratio_warn', 1e-3, ...
    'hist_tail_warn', 2e-3, ...
    'roi_delta_warn', 0.08 ...
);
% trim_fraction: trimmed mean 截尾比例，降低局部亮斑的影响
% max_dn_override: 用于判断饱和的 DN 上限
% qc_thresholds: 汇总 QC 警告阈值

% Aggregation rule used for the final fitted curve.
% 'roi_robust_mean'  : robust mean of absolute ROI statistics
% 'roi_norm_mean'    : old-style comparison, normalize each ROI across LEDs first
% 'roi_norm_robust'  : robust mean after per-ROI normalization
% 'main_roi_only'    : use only the central large ROI
aggregation_mode = 'main_roi_only';
% aggregation_mode 决定最终采用哪条 LED 强度曲线参与拟合：
% main_roi_only   : 只用主 ROI
% roi_robust_mean : 所有 ROI 的绝对值直接做鲁棒平均
% roi_norm_mean   : 每个 ROI 先按自身均值归一化，再平均
% roi_norm_robust : 每个 ROI 先归一化，再做鲁棒平均

% Center selection:
% - leave manual_center_xy empty to estimate center from the corrected stack
% - fill [x, y] to force a known center during debugging
manual_center_xy = [];
auto_center_led_idx = 1:num_img;
center_threshold_quantile = 0.80;
% ROI 中心优先通过自动估计得到：
% 先对所有图像做曝光校正求平均图，再对较亮区域求加权重心
% 若自动中心不可信，可手动指定 manual_center_xy = [x, y]

% Debug controls. These are designed for step-by-step inspection in MATLAB.
debug_mode = true;
debug_led_indices = [];
debug_show_roi_overlay = true;
debug_save_roi_overlay = true;
debug_save_selected_led_figures = false;
debug_print_table = false;
debug_pause_after_roi_build = false;
debug_pause_after_roi_stats = false;
% debug_* 控制调试输出：
% - ROI overlay 图
% - 指定 LED 的单张调试图
% - 工作区调试表格
% - keyboard 暂停点

%% ===== 4. Fitting setup =====
% fit_method = 'smoothing_spline';
fit_method = 'poly4';
fit_spline_p = [];
% fit_method 是一维平滑拟合方法

%% ===== 5. Resolve files =====
raw_root = fullfile(code_dir, 'Raw_input', batch_folder, image_subfolder);
if ~isfolder(raw_root)
    error('Input folder does not exist: %s', raw_root);
end
% 原始数据目录：Code/Raw_input/<batch_folder>/<image_subfolder>

result_dir = fullfile(code_dir, 'Results', 'Calibration', batch_folder);
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end
debug_dir = fullfile(result_dir, 'Debug');
if debug_mode && ~exist(debug_dir, 'dir')
    mkdir(debug_dir);
end
% 结果目录：
% result_dir 存主结果，debug_dir 存 ROI 叠加图和单张 LED 调试图

image_files = strings(num_img, 1);
for i = 1:num_img
    image_files(i) = fullfile(raw_root, "图像_" + string(file_index_list(i)) + ".tif");
    if ~isfile(image_files(i))
        error('Missing image for LED index %d: %s', i, image_files(i));
    end
end

I_first_ch = local_read_first_channel(image_files(1));
image_size = size(I_first_ch);
% 先读取第一张图，只为获得图像尺寸，供 ROI 生成使用

%% ===== 6. Estimate effective bright-zone center =====
% 这一段的目的：
% 不直接假设“有效亮区中心 = 图像几何中心”，
% 而是先从平均校正图估计亮区中心，避免 ROI 放偏导致边缘 LED 偏高或偏低。
% The new workflow used geometric image center by default. If the diffuser
% bright zone is slightly shifted relative to the camera frame, edge LEDs
% can look artificially strong or weak. We therefore estimate the centroid
% of the mean corrected image and explicitly compare it with the geometric
% center.
[auto_center_xy, I_mean_corr] = local_estimate_center_from_stack( ...
    image_files, exposure_ms, dark_noise_level, auto_center_led_idx, center_threshold_quantile);

geom_center_xy = [round((image_size(2) + 1) / 2), round((image_size(1) + 1) / 2)];
if isempty(manual_center_xy)
    roi_center_xy = geom_center_xy;
    % center_source = 'auto_from_mean_image';
    center_source = 'from_geom_center';
else
    roi_center_xy = round(double(manual_center_xy(:)).');
    center_source = 'manual_override';
end

roi_info = build_diffuser_rois(image_size, main_roi_size, aux_roi_size, aux_roi_offsets, roi_center_xy);
% roi_info 中包含：
% 1 个主 ROI
% 若干围绕中心对称分布的辅助 ROI
% 每个 ROI 记录了矩形边界和偏移量

if debug_mode
    fprintf('\n=== ROI center debug ===\n');
    fprintf('Geometric image center     : [x=%d, y=%d]\n', geom_center_xy(1), geom_center_xy(2));
    fprintf('Auto estimated bright center: [x=%.2f, y=%.2f]\n', auto_center_xy(1), auto_center_xy(2));
    fprintf('ROI center used (%s)   : [x=%d, y=%d]\n', center_source, roi_center_xy(1), roi_center_xy(2));
end

if debug_show_roi_overlay || debug_save_roi_overlay
    fig_roi = figure('Color', 'w', 'Name', 'ROI overlay');
    local_plot_image_with_rois(I_mean_corr, roi_info, roi_center_xy, geom_center_xy, auto_center_xy);
    title(sprintf('Mean corrected image with ROIs (%s)', center_source), 'Interpreter', 'none');
end

%% ===== 7. Prepare ROI storage =====
num_aux = numel(roi_info.aux);
num_roi = 1 + num_aux;
% num_roi = 主 ROI 数 1 + 辅助 ROI 数

roi_trimmed_mean = zeros(num_roi, num_img);
roi_plain_mean = zeros(num_roi, num_img);
roi_std = zeros(num_roi, num_img);
roi_median = zeros(num_roi, num_img);
roi_p05 = zeros(num_roi, num_img);
roi_p95 = zeros(num_roi, num_img);

max_dn = zeros(num_img, 1);
sat_ratio = zeros(num_img, 1);
hist_tail_ratio = zeros(num_img, 1);
image_centroid_xy = zeros(num_img, 2);
% 预分配所有统计量数组，避免循环内动态扩容

%% ===== 8. Read images and extract ROI statistics =====
% 核心统计阶段：
% 每张图都执行以下步骤：
% 1. 读图
% 2. 检查最大 DN / 饱和比例 / 直方图尾部
% 3. 扣暗噪声（如果启用）
% 4. 按曝光时间归一化
% 5. 在每个 ROI 上提取多种统计量
% Core diagnosis:
% - If abnormal LEDs are caused by ROI position, then the centroid drift or
%   ROI overlay will usually reveal it immediately.
% - If abnormal LEDs are caused by local hot spots, then mean/trimmed
%   mean/median will disagree.
% - If abnormal LEDs are intrinsic to the data, all ROIs will show the same
%   trend, which is what we currently observe for LEDs 2-5 and 21-25.
for i = 1:num_img
    I_ch = local_read_first_channel(image_files(i));

    max_dn(i) = double(max(I_ch(:)));
    if max_dn_override > 0
        max_possible = max_dn_override;
    else
        max_possible = double(intmax(class(I_ch)));
    end
    sat_ratio(i) = nnz(double(I_ch) >= max_possible) / numel(I_ch);
    hist_tail_ratio(i) = local_hist_tail_ratio(I_ch, max_possible, hist_tail_bins_pct);

    I_corr = double(I_ch);
    if dark_noise_level ~= 0
        I_corr = I_corr - dark_noise_level;
    end
    I_corr(I_corr < 0) = 0;
    I_corr = I_corr ./ exposure_ms(i);
    % 曝光校正后，I_corr 的单位相当于 counts/ms
    % 不同曝光条件下的数据因此可以直接比较

    image_centroid_xy(i, :) = local_weighted_centroid(I_corr, center_threshold_quantile);
    % 对每张图单独估计亮区重心
    % 后面可用来检查不同 LED 下亮区是否明显漂移

    for r = 1:num_roi
        rect = roi_info.all(r).rect;
        patch = I_corr(rect(3):rect(4), rect(1):rect(2));
        patch_vec = patch(:);

        roi_trimmed_mean(r, i) = local_trimmed_mean(patch_vec, trim_fraction);
        roi_plain_mean(r, i) = mean(patch_vec);
        roi_std(r, i) = std(patch_vec);
        roi_median(r, i) = median(patch_vec);
        roi_p05(r, i) = prctile(patch_vec, 5);
        roi_p95(r, i) = prctile(patch_vec, 95);
    end
    % ROI 统计量说明：
    % trimmed_mean: 主统计量，抑制极端值
    % plain_mean  : 普通均值，用于对照
    % std         : ROI 内空间起伏
    % median/p05/p95: 分布形状辅助信息

    if debug_mode && debug_save_selected_led_figures && any(i == debug_led_indices)
        local_save_led_debug_figure(I_corr, roi_info, i, file_index_list(i), ...
            image_centroid_xy(i, :), debug_dir, image_subfolder);
    end
end

if debug_pause_after_roi_build
    fprintf('Paused after ROI build. Inspect roi_info / fig_roi, then type return.\n');
    keyboard;
end
% 若需要逐段调试，可在这里暂停，先检查 ROI 位置是否合理

%% ===== 9. Geometry =====
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
radius_mm = LEDp * hypot(xlocation(:) - xstart, ylocation(:) - ystart);
[kx, ky, NAt] = k_vector(xlocation(:).' - xstart, ylocation(:).' - ystart, ...
    H, LEDp, nglass, t, theta, xint, yint, arraysize^2);
illum_na = NAt(:);
% 将 LED 编号映射到两个几何变量：
% radius_mm: LED 在阵列平面内相对中心 LED 的距离
% illum_na : 该 LED 对应的照明 NA
% 后续会分别在这两个变量下拟合同一条测量曲线

%% ===== 10. ROI summary and aggregation comparison =====
main_trimmed_mean = roi_trimmed_mean(1, :);
aux_trimmed_mean = roi_trimmed_mean(2:end, :);
main_plain_mean = roi_plain_mean(1, :);
aux_plain_mean = roi_plain_mean(2:end, :);
main_roi_std = roi_std(1, :);
aux_roi_std = roi_std(2:end, :);
% 把主 ROI 和辅助 ROI 的统计量拆开，便于单独比较

% The old script first normalized each ROI across LEDs, then averaged over
% ROI. The new script originally averaged absolute ROI values directly.
% We now compute both explicitly, because disagreement between them helps
% distinguish "fixed ROI gain difference" from "true LED intensity change".
roi_trimmed_norm = roi_trimmed_mean ./ mean(roi_trimmed_mean, 2);
% 先对每个 ROI 沿 LED 方向归一化
% 这样不同 ROI 的“绝对亮度差”会被消除，只保留相对变化趋势

curve_main_roi_only = local_normalize_curve(main_trimmed_mean(:));
curve_roi_robust_mean = local_normalize_curve(local_column_trimmed_mean(roi_trimmed_mean, trim_fraction));
curve_roi_norm_mean = local_normalize_curve(mean(roi_trimmed_norm, 1).');
curve_roi_norm_robust = local_normalize_curve(local_column_trimmed_mean(roi_trimmed_norm, trim_fraction));
% 四条候选曲线：
% main_roi_only   : 只看中心主 ROI
% roi_robust_mean : 所有 ROI 的绝对值直接鲁棒平均
% roi_norm_mean   : 每个 ROI 先归一化，再做普通平均
% roi_norm_robust : 每个 ROI 先归一化，再做鲁棒平均

switch lower(string(aggregation_mode))
    case "main_roi_only"
        led_intensity_measured_norm = curve_main_roi_only;
    case "roi_robust_mean"
        led_intensity_measured_norm = curve_roi_robust_mean;
    case "roi_norm_mean"
        led_intensity_measured_norm = curve_roi_norm_mean;
    case "roi_norm_robust"
        led_intensity_measured_norm = curve_roi_norm_robust;
    otherwise
        error('Unknown aggregation_mode: %s', aggregation_mode);
end

led_intensity_measured = led_intensity_measured_norm * mean(local_column_trimmed_mean(roi_trimmed_mean, trim_fraction));
roi_cross_std = std(roi_trimmed_norm, 0, 1).';
% roi_cross_std 描述同一 LED 在不同 ROI 之间的一致性
% 若某些 LED 的 roi_cross_std 很大，说明 ROI 间差异明显，需要检查亮区均匀性或 ROI 位置

qc = summarize_led_qc(max_dn, sat_ratio, hist_tail_ratio, ...
    main_plain_mean, main_roi_std, aux_plain_mean, aux_roi_std, qc_thresholds);
qc.roi_cross_std = roi_cross_std;
qc.image_centroid_xy = image_centroid_xy;
qc.centroid_shift_px = sqrt(sum((image_centroid_xy - roi_center_xy).^2, 2));
% centroid_shift_px 是每张图的亮区重心相对 ROI 中心的偏移量
% 如果某些 LED 偏移很大，优先怀疑 ROI 放置不合适

aggregation_compare = struct();
aggregation_compare.main_roi_only = curve_main_roi_only;
aggregation_compare.roi_robust_mean = curve_roi_robust_mean;
aggregation_compare.roi_norm_mean = curve_roi_norm_mean;
aggregation_compare.roi_norm_robust = curve_roi_norm_robust;
aggregation_compare.selected_mode = aggregation_mode;
aggregation_compare.selected_curve = led_intensity_measured_norm;
% 保存所有汇总方式，后续可以离线比较而不必重新跑图像提取

if debug_print_table
    debug_table = table((1:num_img).', file_index_list(:), exposure_ms(:), ...
        image_centroid_xy(:, 1), image_centroid_xy(:, 2), ...
        main_trimmed_mean(:), curve_main_roi_only(:), curve_roi_robust_mean(:), ...
        curve_roi_norm_mean(:), qc.centroid_shift_px(:), ...
        'VariableNames', {'led_idx', 'file_idx', 'exposure_ms', 'centroid_x', 'centroid_y', ...
        'main_trimmed', 'main_norm', 'robust_norm', 'roi_norm_mean', 'centroid_shift_px'});
    disp(debug_table);
end
% debug_table 会直接出现在 MATLAB 工作区中，适合快速筛查异常 LED

if debug_pause_after_roi_stats
    fprintf('Paused after ROI statistics. Inspect roi_* variables, then type return.\n');
    keyboard;
end
% 若要详细比较不同统计量，可在这里暂停

%% ===== 11. Smooth fits =====
fit_radius = fit_led_intensity_smooth(led_intensity_measured_norm, radius_mm, ...
    'exclude_idx', exclude_idx, ...
    'normalize_mode', 'none', ...
    'fit_method', fit_method, ...
    'spline_p', fit_spline_p, ...
    'label', 'radius_mm');

fit_na = fit_led_intensity_smooth(led_intensity_measured_norm, illum_na, ...
    'exclude_idx', exclude_idx, ...
    'normalize_mode', 'none', ...
    'fit_method', fit_method, ...
    'spline_p', fit_spline_p, ...
    'label', 'illum_na');
% 用同一条实测曲线分别做：
% 1. 相对强度 vs radius 的拟合
% 2. 相对强度 vs illumination NA 的拟合

%% ===== 12. Package results =====
calib = struct();
% 将所有参数、几何量、ROI 统计量、QC、拟合结果打包到 calib 结构体
calib.metadata = struct( ...
    'batch_folder', char(batch_folder), ...
    'image_subfolder', char(image_subfolder), ...
    'file_index_list', file_index_list, ...
    'dark_noise_level', dark_noise_level, ...
    'trim_fraction', trim_fraction, ...
    'fit_method', fit_method, ...
    'fit_spline_p', fit_spline_p, ...
    'normalize_mode', 'none', ...
    'aggregation_mode', aggregation_mode, ...
    'center_source', center_source ...
);
calib.geometry = struct( ...
    'arraysize', arraysize, ...
    'LEDp', LEDp, ...
    'H', H, ...
    'nglass', nglass, ...
    't', t, ...
    'theta', theta, ...
    'xstart', xstart, ...
    'ystart', ystart, ...
    'xint', xint, ...
    'yint', yint, ...
    'xlocation', xlocation(:), ...
    'ylocation', ylocation(:), ...
    'kx', kx(:), ...
    'ky', ky(:), ...
    'illum_na', illum_na, ...
    'radius_mm', radius_mm, ...
    'geom_center_xy', geom_center_xy, ...
    'auto_center_xy', auto_center_xy, ...
    'roi_center_xy', roi_center_xy ...
);
calib.qc = qc;
calib.roi = struct( ...
    'roi_info', roi_info, ...
    'trimmed_mean', roi_trimmed_mean, ...
    'plain_mean', roi_plain_mean, ...
    'std', roi_std, ...
    'median', roi_median, ...
    'p05', roi_p05, ...
    'p95', roi_p95, ...
    'trimmed_norm', roi_trimmed_norm, ...
    'main_trimmed_mean', main_trimmed_mean(:), ...
    'aux_trimmed_mean', aux_trimmed_mean, ...
    'roi_cross_std', roi_cross_std ...
);
calib.measured = struct( ...
    'exposure_ms', exposure_ms(:), ...
    'led_intensity_measured', led_intensity_measured(:), ...
    'led_intensity_measured_norm', led_intensity_measured_norm(:), ...
    'exclude_idx', exclude_idx(:) ...
);
calib.aggregation_compare = aggregation_compare;
calib.fit_radius = fit_radius;
calib.fit_na = fit_na;
calib.led_intensity_measured_norm = led_intensity_measured_norm(:);
calib.led_intensity_fit_radius_norm = fit_radius.y_fitted(:);
calib.led_intensity_fit_na_norm = fit_na.y_fitted(:);
% 顶层保留常用曲线，便于其他脚本直接读取，不必层层索引

%% ===== 13. Save and plot =====
file_tag = sprintf('%s_%s_%dx%d_poly4', char(batch_folder), local_safe_name(image_subfolder), arraysize, arraysize);
save(fullfile(result_dir, "diffuser_calib_" + file_tag + ".mat"), 'calib');
% 保存主结果 mat 文件

if debug_mode && debug_save_roi_overlay
    saveas(fig_roi, fullfile(debug_dir, "roi_overlay_" + file_tag + ".png"));
end
% 保存 ROI 叠加图，用于直观核对 ROI 是否压在亮区上

fig1 = figure('Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% fig1: 主结果概览
% 包含：
% - 最终测量曲线与两种拟合
% - 不同 ROI 汇总方式对比
% - 各 ROI 归一化曲线
% - radius / NA 域下的散点与拟合

nexttile;
hold on; grid on; box on;
plot(1:num_img, led_intensity_measured_norm, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', [0.8 0.8 0.8]);
plot(1:num_img, fit_na.y_fitted, 'r-', 'LineWidth', 1.8);
plot(1:num_img, fit_radius.y_fitted, 'b--', 'LineWidth', 1.6);
xlabel('LED index');
ylabel('Relative intensity');
title(sprintf('Selected aggregation: %s', aggregation_mode), 'Interpreter', 'none');
legend('Measured', 'Fit by NA', 'Fit by radius', 'Location', 'best');

% nexttile;
% hold on; grid on; box on;
% plot(1:num_img, curve_main_roi_only, 'k-o', 'LineWidth', 1.0);
% plot(1:num_img, curve_roi_robust_mean, 'r-s', 'LineWidth', 1.0);
% plot(1:num_img, curve_roi_norm_mean, 'b-^', 'LineWidth', 1.0);
% plot(1:num_img, curve_roi_norm_robust, 'g-d', 'LineWidth', 1.0);
% xlabel('LED index');
% ylabel('Relative intensity');
% title('Aggregation comparison');
% legend('Main ROI only', 'ROI robust mean', 'ROI norm mean', 'ROI norm robust', 'Location', 'best');

nexttile;
hold on; grid on; box on;
plot(1:num_img, roi_trimmed_norm.', 'LineWidth', 0.8);
xlabel('LED index');
ylabel('ROI-normalized intensity');
title('Per-ROI normalized curves');

nexttile;
hold on; grid on; box on;
plot(radius_mm, fit_radius.y_measured, 'ko', 'MarkerFaceColor', [0.7 0.7 0.7]);
plot(radius_mm, fit_radius.y_fitted, 'b*');
xlabel('Radius (mm)');
ylabel('Relative intensity');
title('Radius-domain fit');

nexttile;
hold on; grid on; box on;
plot(illum_na, fit_na.y_measured, 'ko', 'MarkerFaceColor', [0.7 0.7 0.7]);
plot(illum_na, fit_na.y_fitted, 'r*');
xlabel('Illumination NA');
ylabel('Relative intensity');
title('NA-domain fit');

% nexttile;
% hold on; grid on; box on;
% plot(1:num_img, qc.sat_ratio, 'r-o', 'LineWidth', 1.2);
% plot(1:num_img, qc.hist_tail_ratio, 'b-s', 'LineWidth', 1.2);
% plot(1:num_img, qc.centroid_shift_px, 'k-^', 'LineWidth', 1.2);
% xlabel('LED index');
% ylabel('QC metric');
% title('Saturation / histogram tail / centroid shift');
% legend('sat ratio', 'tail ratio', 'centroid shift (px)', 'Location', 'best');
saveas(fig1, fullfile(result_dir, "diffuser_summary_" + file_tag + ".png"));

fig2 = figure('Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% fig2: 诊断图
% 包括拟合残差、同组离散度、图像亮区重心漂移

nexttile;
hold on; grid on; box on;
plot(1:num_img, fit_na.residual, 'o-', 'LineWidth', 1.2);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('LED index');
ylabel('Relative residual');
title('Residuals of NA fit');

nexttile;
hold on; grid on; box on;
plot(fit_radius.group_x, fit_radius.same_x_std, 'o-', 'LineWidth', 1.2);
% plot(fit_na.group_x, fit_na.same_x_std, 's-', 'LineWidth', 1.2);
xlabel('Grouped geometry variable');
ylabel('Std within same group');
title('Same-radius / same-NA spread');
legend('Radius grouped std', 'NA grouped std', 'Location', 'best');

% nexttile;
% hold on; grid on; box on;
% plot(1:num_img, image_centroid_xy(:, 1), 'r-o', 'LineWidth', 1.0);
% plot(1:num_img, image_centroid_xy(:, 2), 'b-s', 'LineWidth', 1.0);
% yline(roi_center_xy(1), 'r--');
% yline(roi_center_xy(2), 'b--');
% xlabel('LED index');
% ylabel('Centroid coordinate (px)');
% title('Image centroid drift');
% legend('centroid x', 'centroid y', 'ROI center x', 'ROI center y', 'Location', 'best');

saveas(fig2, fullfile(result_dir, "diffuser_residuals_" + file_tag + ".png"));

fprintf('\nCalibration saved to: %s\n', result_dir);
fprintf('Recommended first checks for abnormal LEDs:\n');
fprintf('1) roi_overlay_*.png\n');
fprintf('2) Debug/led_debug_*.png for LEDs %s\n', mat2str(debug_led_indices));
fprintf('3) debug_table in MATLAB workspace\n');
% 异常结果优先看这三类输出：ROI 叠加图、单张 LED 调试图、debug_table

%% ===== Local functions =====
function I_ch = local_read_first_channel(image_file)
% 读取 TIFF 第一通道，兼容灰度图和 RGB 图像
I_raw = imread(image_file);
if ndims(I_raw) >= 3
    I_ch = I_raw(:, :, 1);
else
    I_ch = I_raw;
end
end

function [center_xy, I_mean_corr] = local_estimate_center_from_stack(image_files, exposure_ms, dark_noise_level, use_idx, threshold_quantile)
% 将多张图做曝光校正后求平均图，再对平均图估计亮区重心
sample_count = numel(use_idx);
for k = 1:sample_count
    img_idx = use_idx(k);
    I_ch = local_read_first_channel(image_files(img_idx));
    I_corr = double(I_ch);
    if dark_noise_level ~= 0
        I_corr = I_corr - dark_noise_level;
    end
    I_corr(I_corr < 0) = 0;
    I_corr = I_corr ./ exposure_ms(img_idx);
    if k == 1
        I_sum = zeros(size(I_corr));
    end
    I_sum = I_sum + I_corr;
end
I_mean_corr = I_sum / sample_count;
center_xy = local_weighted_centroid(I_mean_corr, threshold_quantile);
end

function center_xy = local_weighted_centroid(I_in, threshold_quantile)
% 只对高于某个亮度分位数的像素求加权重心，减少背景噪声影响
I_use = double(I_in);
thr = quantile(I_use(:), threshold_quantile);
mask = I_use >= thr;
weights = I_use .* mask;
if ~any(weights(:))
    weights = I_use;
end
[yy, xx] = ndgrid(1:size(I_use, 1), 1:size(I_use, 2));
wsum = sum(weights(:));
center_x = sum(sum(weights .* xx)) / max(wsum, eps);
center_y = sum(sum(weights .* yy)) / max(wsum, eps);
center_xy = [center_x, center_y];
end

function local_plot_image_with_rois(I_img, roi_info, roi_center_xy, geom_center_xy, auto_center_xy)
% 在平均图上画出所有 ROI，并同时标出几何中心、自动中心和最终使用中心
imagesc(I_img);
axis image;
colormap gray;
colorbar;
hold on;
for r = 1:numel(roi_info.all)
    rect = roi_info.all(r).rect;
    rectangle('Position', [rect(1), rect(3), rect(2) - rect(1) + 1, rect(4) - rect(3) + 1], ...
        'EdgeColor', local_roi_color(r), 'LineWidth', 1.2);
    text(rect(1), rect(3) - 10, roi_info.all(r).name, 'Color', local_roi_color(r), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
end
plot(geom_center_xy(1), geom_center_xy(2), 'wo', 'MarkerSize', 9, 'LineWidth', 1.5);
plot(auto_center_xy(1), auto_center_xy(2), 'g+', 'MarkerSize', 12, 'LineWidth', 2.0);
plot(roi_center_xy(1), roi_center_xy(2), 'rx', 'MarkerSize', 12, 'LineWidth', 2.0);
legend('ROIs', 'Geom center', 'Auto center', 'ROI center', 'Location', 'best');
hold off;
end

function local_save_led_debug_figure(I_corr, roi_info, led_idx, file_idx, centroid_xy, debug_dir, image_subfolder)
% 保存指定 LED 的调试图：
% 左侧是校正后图像与 ROI 位置，右侧是像素直方图
if ~exist(debug_dir, 'dir')
    mkdir(debug_dir);
end
fig = figure('Color', 'w', 'Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(I_corr);
axis image;
colormap gray;
colorbar;
hold on;
for r = 1:numel(roi_info.all)
    rect = roi_info.all(r).rect;
    rectangle('Position', [rect(1), rect(3), rect(2) - rect(1) + 1, rect(4) - rect(3) + 1], ...
        'EdgeColor', local_roi_color(r), 'LineWidth', 1.0);
end
plot(centroid_xy(1), centroid_xy(2), 'r+', 'MarkerSize', 12, 'LineWidth', 1.8);
title(sprintf('LED %d, file %d, %s', led_idx, file_idx, image_subfolder), 'Interpreter', 'none');
hold off;

nexttile;
histogram(I_corr(:), 128, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
grid on;
xlabel('Counts / ms');
ylabel('Probability');
title('Pixel histogram');

saveas(fig, fullfile(debug_dir, sprintf('led_debug_%02d_file_%02d.png', led_idx, file_idx)));
close(fig);
end

function color = local_roi_color(idx)
% 为不同 ROI 分配可区分颜色
palette = lines(12);
color = palette(mod(idx - 1, size(palette, 1)) + 1, :);
end

function col_mean = local_column_trimmed_mean(data_matrix, trim_fraction)
% 对矩阵每一列分别求 trimmed mean
num_col = size(data_matrix, 2);
col_mean = zeros(num_col, 1);
for ii = 1:num_col
    col_mean(ii) = local_trimmed_mean(data_matrix(:, ii), trim_fraction);
end
end

function y_norm = local_normalize_curve(y)
% 归一化到均值为 1，便于不同曲线比较形状
y = double(y(:));
y_norm = y ./ mean(y);
end

function ratio = local_hist_tail_ratio(I_ch, max_possible, hist_tail_bins_pct)
% 用直方图高亮端的尾部占比，粗略判断是否接近饱和
hist_max = max_possible + 50;
num_bins = 256;
edges = linspace(0, hist_max, num_bins + 1);
data_use = double(I_ch(:));
data_use = data_use(data_use >= 10);
counts = histcounts(data_use, edges);
tail_bins = max(1, round(num_bins * hist_tail_bins_pct / 100));
if sum(counts) > 0
    ratio = sum(counts(end-tail_bins+1:end)) / sum(counts);
else
    ratio = 0;
end
end

function value = local_trimmed_mean(data, trim_fraction)
% 简单截尾均值：去掉两端一定比例极端值后再求均值
data = sort(double(data(:)));
n = numel(data);
trim_n = floor(n * trim_fraction / 2);
if 2 * trim_n >= n
    value = mean(data);
else
    value = mean(data(trim_n + 1:n - trim_n));
end
end

function text_out = local_safe_name(text_in)
% 将任意字符串转换为相对安全的文件名片段
text_out = regexprep(char(text_in), '[^\w\+\-]', '_');
end
