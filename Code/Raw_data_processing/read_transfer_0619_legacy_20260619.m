close all; clear; clc;

%% ========== 1. 参数设置 ==========
spsize = 0.69e-6;
% 设置高一档的分辨率
spsize = spsize / 2;

m1 = 600;
n1 = 600;

m = m1;
n = n1;

dark_noise_level = 0;

wlength = 6.28e-07;
NA = 0.15;
z = 2e-6;

%% ========== 2. 读入实验图像并裁剪多个 ROI ==========
num_img = 25;

set(0, 'DefaultFigureVisible', 'on');

% exposure setting for 0619 / z=_0.050_high
exposure_ms = 100.65 * ones(1, num_img);
exposure_ms(1:9) = 42.307;
exposure_ms([13, 17, 21, 25]) = 201.301;

str_z = "z=_0.050_high_2nd";
batch_folder = "0619";
raw_root = "Raw_input";

% Multi-slice ROI configuration.
% ptr determines the saved file suffix:
%   FPM_RawData_z <str_z><ptr>.mat
% row_shift / col_shift are offsets relative to the base ROI.
% Add more entries here if you want more slices.
roi_configs = struct( ...
    'ptr', {"_up", "", "_down"}, ...
    'row_shift', {-m * 6 / 6, 0, m}, ...
    'col_shift', {0, 0, 0});

num_roi = numel(roi_configs);
im_high_HDR_list = cell(1, num_roi);
for roi_idx = 1:num_roi
    im_high_HDR_list{roi_idx} = zeros(m, n, num_img, 'double');
end

code_root = fileparts(fileparts(mfilename('fullpath')));
batch_dir = fullfile(code_root, raw_root, batch_folder);

for i = 1:num_img
    data_path = fullfile(batch_dir, str_z, "图像_" + string(i) + ".tif");
    I_hr = imread(data_path);
    I_hr = double(I_hr(:,:,1));

    fprintf('图像 %s 的尺寸为 %d x %d\n', data_path, size(I_hr, 1), size(I_hr, 2));
    fprintf('图像 %d 的最大强度为 %.0f\n', i, max(I_hr(:)));

    I_hr = I_hr - dark_noise_level;
    I_hr(I_hr < 0) = 0;

    center_y = round(size(I_hr, 1) / 2);
    center_x = round(size(I_hr, 2) / 2);

    base_row_start = center_y - m * 9 / 6 + 1;
    base_row_end = center_y + m * -3 / 6;
    base_col_start = center_x - n * 3 / 3 + 1;
    base_col_end = center_x + n * 0 / 3;

    for roi_idx = 1:num_roi
        row_start = base_row_start + roi_configs(roi_idx).row_shift;
        row_end = base_row_end + roi_configs(roi_idx).row_shift;
        col_start = base_col_start + roi_configs(roi_idx).col_shift;
        col_end = base_col_end + roi_configs(roi_idx).col_shift;

        if row_start < 1 || row_end > size(I_hr, 1) || ...
                col_start < 1 || col_end > size(I_hr, 2)
            error('ROI %s 超出图像范围: rows [%d, %d], cols [%d, %d].', ...
                roi_configs(roi_idx).ptr, row_start, row_end, col_start, col_end);
        end

        I_crop = I_hr(row_start:row_end, col_start:col_end);
        I_crop = I_crop ./ exposure_ms(i);

        amp = sqrt(I_crop);
        im_high_HDR_list{roi_idx}(:,:,i) = amp;

        % if (i == 1 || i == 5 || i == 10)
        %     figure('Name', '基准真值 Ground Truth (Ideal)');
        %     imshow(abs(amp), []); title("HR 振幅 (真实数据)" + num2str(i));
        % end

    end

    % figure('Name', '基准真值 Ground Truth (Ideal)');
    % imshow(I_hr, []); title('HR 振幅 (真实数据)');
    % filename = sprintf('LED_%d_高分辨率原图.png', i);
    % saveas(gcf, fullfile(batch_dir, filename));
end

for roi_idx = 1:num_roi
    ptr = roi_configs(roi_idx).ptr;
    roi_row_shift = roi_configs(roi_idx).row_shift;
    roi_col_shift = roi_configs(roi_idx).col_shift;
    im_high_HDR = im_high_HDR_list{roi_idx};

    out_name = "FPM_RawData_z " + str_z + ptr + ".mat";
    save(fullfile(batch_dir, out_name), ...
        'z', 'NA', 'wlength', 'spsize', 'im_high_HDR', ...
        'ptr', 'roi_row_shift', 'roi_col_shift');

    fprintf('已保存 %s\n', out_name);
end
