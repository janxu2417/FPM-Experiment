%% 模拟 FP 成像过程：加入部分相干与 pupil/MTF 高频衰减
clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'function'));

%% 1. 基本参数
spsize = 0.69e-6;
pratio = 5;
psize = spsize / pratio;

m1 = 600;
n1 = 600;
m = m1 * pratio;
n = n1 * pratio;

dark_noise_level = 0;
wlength = 628e-9;
NA_obj = 0.15;
z = 0;

sim_cfg = fpm_sim_nonideal_config();

%% 2. 导入高分辨率靶标并构造复振幅
load(fullfile('Raw_input', 'FP_ideal_data.mat'), 'O_ideal');
I_hr = centerMatrix(O_ideal, m, n);
I_hr = max(I_hr, 0.001);

fprintf('HR 靶标最大强度: %.6f\n', max(I_hr(:)));
figure('Name', 'HR target');
imshow(I_hr, []);
title('Ground truth intensity');

I_hr = I_hr - dark_noise_level;
I_hr(I_hr < 0) = 0;
I_hr = I_hr ./ max(I_hr(:));

center_y = round(size(I_hr, 1) / 2);
center_x = round(size(I_hr, 2) / 2);
I_crop = I_hr(center_y - m/2 + 1:center_y + m/2, center_x - n/2 + 1:center_x + n/2);

amp = sqrt(I_crop);
phase = (1 - amp) * pi;
O = amp .* exp(1j * phase);

figure('Name', 'Ground truth complex field');
subplot(1, 2, 1); imshow(abs(O), []); title('Amplitude');
subplot(1, 2, 2); imshow(angle(O), []); title('Phase');

%% 3. LED 阵列参数
xstart = 18;
ystart = 20;
arraysize = 7;
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);

H = 40.24;
LEDp = 4;
nglass = 1.52;
t = 1;
xint = 0.0;
yint = 0.0;
theta = 1.0;

[kx_norm, ky_norm, NAt] = k_vector( ...
    xlocation - xstart, ylocation - ystart, ...
    H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

LED_num = arraysize^2;
k0 = 2 * pi / wlength;
dkx = 2 * pi / (psize * n);
dky = 2 * pi / (psize * m);

%% 4. 基础 pupil 与 MTF 衰减
[KX, KY] = meshgrid( ...
    (-n1/2:n1/2-1) * (2 * pi / (n1 * spsize)), ...
    (-m1/2:m1/2-1) * (2 * pi / (m1 * spsize)));
K = sqrt(KX.^2 + KY.^2);
kc = NA_obj * k0;
rho = K ./ max(kc, eps);

pupil_support = double(K <= kc);
if sim_cfg.pupil_mtf.enable
    pupil_amp = exp(-((rho ./ sim_cfg.pupil_mtf.rolloff_sigma) .^ sim_cfg.pupil_mtf.rolloff_power));
else
    pupil_amp = ones(size(rho));
end
pupil_amp = pupil_amp .* pupil_support;
P = pupil_amp .* exp(1j * K.^2 * z / (2 * k0));

%% 5. 正向模型
O_ft = fftshift(fft2(O));
imseqlow1 = zeros(m1, n1, LED_num);

fprintf('开始生成非理想 FPM 模拟序列...\n');
for i = 1:LED_num
    led_center_mm = [ ...
        (xlocation(i) - xstart) * LEDp + xint, ...
        (ylocation(i) - ystart) * LEDp + yint];
    sample_list = local_build_led_subsources(led_center_mm, sim_cfg.partial_coherence);

    intensity_acc = zeros(m1, n1);
    for jj = 1:size(sample_list, 1)
        sample_x_mm = sample_list(jj, 1);
        sample_y_mm = sample_list(jj, 2);
        sample_w = sample_list(jj, 3);

        [kx_loc, ky_loc, ~] = k_vector( ...
            sample_x_mm / LEDp, sample_y_mm / LEDp, ...
            H, LEDp, nglass, t, theta, 0.0, 0.0, 1);

        kxc = round((n + 1) / 2 - (kx_loc * k0) / dkx);
        kyc = round((m + 1) / 2 - (ky_loc * k0) / dky);
        kyl = round(kyc - (m1 - 1) / 2);
        kyh = round(kyc + (m1 - 1) / 2);
        kxl = round(kxc - (n1 - 1) / 2);
        kxh = round(kxc + (n1 - 1) / 2);

        if kyl < 1 || kxl < 1 || kyh > m || kxh > n
            continue;
        end

        O_j = O_ft(kyl:kyh, kxl:kxh);
        V = ifft2(ifftshift(O_j .* P));
        intensity_acc = intensity_acc + sample_w * abs(V).^2;
    end

    im_amp = sqrt(max(intensity_acc, 0));
    if sim_cfg.partial_coherence.enable_image_blur
        im_amp = imgaussfilt(im_amp, sim_cfg.partial_coherence.image_sigma_px, 'FilterSize', sim_cfg.partial_coherence.filter_size_px);
    end
    imseqlow1(:, :, i) = im_amp;
end
fprintf('序列生成完成。\n');

save('FP_input_data_Sim.mat', ...
    'imseqlow1', 'theta', 'spsize', 'wlength', 'NA_obj', 'z', 'xint', 'yint', ...
    'sim_cfg', 'pupil_amp');

%% 6. 简单对比输出
figure('Name', 'System resolution comparison');
subplot(1, 2, 1); imshow(abs(O), []); title('Ground truth');
subplot(1, 2, 2); imshow(imseqlow1(:, :, 1), []); title('Center brightfield');

if exist(fullfile('Raw_input', '0619', 'FPM_RawData_z z=_0.050_7_preproc_v1.mat'), 'file')
    load(fullfile('Raw_input', '0619', 'FPM_RawData_z z=_0.050_7_preproc_v1.mat'), 'im_high_HDR');
    figure('Name', 'Center brightfield vs experiment');
    subplot(1, 2, 1); imshow(imseqlow1(:, :, 1), []); title('Simulation');
    subplot(1, 2, 2); imshow(im_high_HDR(:, :, 1), []); title('Experiment');
end

%% 7. 分辨率 benchmark
result = FPM_resolution_benchmark( ...
    'FP_input_data_Sim.mat', ...
    fullfile('Raw_input', '0619', 'FPM_RawData_z z=_0.050_7_preproc_v1.mat'));
disp(result.summary_table);

function sample_list = local_build_led_subsources(led_center_mm, cfg)
if ~cfg.enable
    sample_list = [led_center_mm(1), led_center_mm(2), 1];
    return;
end

grid_n = cfg.grid_n;
half_size = cfg.source_size_mm / 2;
axis_mm = linspace(-half_size, half_size, grid_n);
[dx, dy] = meshgrid(axis_mm, axis_mm);

if cfg.gaussian_weight_sigma_mm > 0
    w = exp(-0.5 * ((dx / cfg.gaussian_weight_sigma_mm).^2 + (dy / cfg.gaussian_weight_sigma_mm).^2));
else
    w = ones(size(dx));
end
w = w / sum(w(:));

sample_list = [ ...
    led_center_mm(1) + dx(:), ...
    led_center_mm(2) + dy(:), ...
    w(:)];
end
