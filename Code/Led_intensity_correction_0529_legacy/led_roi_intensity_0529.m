close all; clear; clc;
set(0, 'DefaultFigureVisible', 'on');

%% ===== 1. Params =====
num_img = 25;
batch_folder = "0529";
str_z = "光强修正";
raw_root = "Raw_input";
root = "D:\1_徐前\物理\傅里叶叠层显微术\Code";

% Camera offset (set to 0 if unknown)
dark_noise_level = 0;

% Histogram-based saturation check
max_dn_override = 900; % set 0 to use intmax(class)
hist_tail_bins_pct = 1; % top 1% bins used as tail
hist_tail_warn_pct = 0.2; % warn if tail ratio exceeds this percent

% ROI definitions on the raw (full-size) image: [x1 x2 y1 y2]
rois = [
    1010 1110 1000 1100;
    1210 1310 1000 1100;
    820  920  1000 1100;
    1010 1110 800  900;
    1010 1110 1200 1300
];

%% ===== 2. Exposure time (ms) =====
% LED 1-9: 50 ms
% LED 10-25: 200 ms
% LED 11, 23, 24: 100 ms
exposure_ms = 200 * ones(1, num_img);
exposure_ms(1:9) = 50;
exposure_ms([11, 23, 24]) = 100;

%% ===== 3. Read images and compute ROI mean intensity =====
roi_mean_raw = zeros(size(rois, 1), num_img);

for i = 1:num_img
    data_name = batch_folder + "\" + str_z + "\图像_" + int2str(i + 10);
    data_dir = fullfile(root, raw_root, data_name + ".tif");

    I_raw = imread(data_dir);
    I_ch = I_raw(:,:,1);
    if i == 1
        % figure('Color', 'w');
        % imshow(I_ch, []);
        % title('LED1 raw image (channel 1)');

        fprintf('这幅图的物理最大光强值是: %d\n', max(I_ch(:)));
        if max_dn_override > 0
            max_possible = max_dn_override;
        else
            max_possible = double(intmax(class(I_ch)));
        end
        hist_max = max_possible + 50;
        sat_ratio = nnz(I_ch >= max_possible) / numel(I_ch);
        fprintf('LED1 saturation ratio (at max DN): %.4f%%\n', sat_ratio * 100);
        if sat_ratio > 0.1
            warning('LED1 may be saturated (%.4f%% pixels at max DN).', sat_ratio * 100);
        end

        num_bins = 256;
        edges = linspace(0, hist_max, num_bins + 1);
        data_use = double(I_ch(:));
        data_use = data_use(data_use >= 10);
        counts = histcounts(data_use, edges);
        tail_bins = max(1, round(num_bins * hist_tail_bins_pct / 100));
        if sum(counts) > 0
            tail_ratio = sum(counts(end-tail_bins+1:end)) / sum(counts);
        else
            tail_ratio = 0;
        end
        fprintf('LED1 histogram tail ratio (top %.2f%% bins): %.4f%%\n', ...
            hist_tail_bins_pct, tail_ratio * 100);
        if tail_ratio > hist_tail_warn_pct / 100
            warning('LED1 may be saturated (%.4f%% in top %.2f%% bins).', ...
                tail_ratio * 100, hist_tail_bins_pct);
        end

        figure('Color', 'w', 'Visible', 'on');
        h = histogram(data_use, edges, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
        xlim([edges(1), edges(end)]);
        grid on;
        xlabel('DN');
        ylabel('Count');
        title('LED1 histogram (raw)');
        if any(h.Values)
            ylim([0, max(h.Values) * 1.1]);
        end
    end

    I = double(I_ch);
    I = I - dark_noise_level;
    I(I < 0) = 0;

    for r = 1:size(rois, 1)
        x1 = rois(r, 1); x2 = rois(r, 2);
        y1 = rois(r, 3); y2 = rois(r, 4);
        patch = I(y1:y2, x1:x2);
        roi_mean_raw(r, i) = mean(patch(:));
    end
end

%% ===== 4. Exposure correction =====
roi_mean_corr = roi_mean_raw ./ reshape(exposure_ms, 1, []);

%% ===== 5. Normalize and average across ROIs =====
roi_norm = roi_mean_corr ./ mean(roi_mean_corr, 2);
led_intensity_avg = mean(roi_norm, 1);
led_intensity_avg_norm = led_intensity_avg ./ mean(led_intensity_avg);

%% ===== 6. Save results =====
out_name = "LED_intensity_roi_0529.mat";
out_dir = fullfile(root, raw_root, batch_folder);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save(fullfile(out_dir, out_name), ...
    'roi_mean_raw', 'roi_mean_corr', 'roi_norm', ...
    'led_intensity_avg', 'led_intensity_avg_norm', ...
    'rois', 'exposure_ms');

%% ===== 7. Plot =====
led_idx = 1:num_img;
figure('Color', 'w');
hold on; grid on;
for r = 1:size(rois, 1)
    plot(led_idx, roi_norm(r, :), '-o', 'LineWidth', 1.5);
end
xlabel('LED index');
ylabel('Intensity (counts/ms)');
title('ROI intensity vs LED index (exposure corrected)');
legend('ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI5', 'Location', 'best');

figure('Color', 'w');
plot(led_idx, led_intensity_avg_norm, '-o', 'LineWidth', 1.5);
grid on;
xlabel('LED index');
ylabel('Normalized mean intensity');
title('Mean LED intensity (ROI-normalized, averaged)');
