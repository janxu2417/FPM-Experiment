close all; clear; clc;

%% FP_recover_code
% FPM 重建的一键入口：
% 1. 用 cfg 统一管理输入数据、光强校正、重建参数和输出命名
% 2. 默认直接读取新的预处理输出：
%    FPM_RawData_z z=_0.050_preproc_v1.mat
% 3. 保留必要的显示/调试开关，避免主流程继续混用旧文件名
% 4. 与预处理脚本共享同一个配置源：
%    Code/function/fpm_0619_shared_config.m
%
% 这里的 cfg 与 read_transfer_0619.m 保持同风格：
%   cfg = struct('字段名', 字段值, ...)
% 使用时通过点号访问：
%   cfg.data_name
%   cfg.recon.upsmp_ratio
%   cfg.calibration.calib_file_name
%
% 如果后续新增另一套重建入口，优先在 fpm_0619_shared_config.m 里新增 case。

%% ===== 1. User entry =====
% 默认入口：读取 read_transfer_0619.m 当前默认生成的新 .mat 文件
active_preset = "0619_zneg0050_7_preproc_v3";
% active_preset = "0619_zneg0050_7_led_corr_v2"
% 共享配置函数位于 Code/function 下。
% 为了保证从不同工作目录运行时都能找到它，这里先把 Code 整体加到路径里。
code_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(code_root));

% 可选示例：
% active_preset = "0619_zneg0050_preproc_v1_no_calib";
% active_preset = "0619_zneg0050_preproc_v1_up";

cfg = fpm_0619_shared_config("recover", active_preset);

%% ===== 2. Path setup =====
set(0, 'DefaultFigureVisible', 'on');
data_dir = fullfile(code_root, cfg.raw_root, cfg.batch_folder, cfg.data_name);
result_dir = fullfile(code_root, cfg.result_root, cfg.batch_folder);

if ~isfile(data_dir)
    error('Input .mat file does not exist: %s', data_dir);
end

if ~isfolder(result_dir)
    mkdir(result_dir);
end

fprintf('\n=== FPM recovery preset ===\n');
fprintf('preset       : %s\n', cfg.preset_name);
fprintf('input .mat   : %s\n', data_dir);
fprintf('result dir   : %s\n', result_dir);
fprintf('result tag   : %s\n\n', cfg.result_tag);

%% ===== 3. Load preprocessed data =====
% 约定：预处理脚本保存的主数据字段名是 im_high_HDR。
% 在重建入口中统一映射到 imlow_HDR，保持与 himrecover 接口兼容。
data_struct = load(data_dir);
if ~isfield(data_struct, 'im_high_HDR')
    error('Input .mat does not contain im_high_HDR: %s', data_dir);
end

imlow_HDR = data_struct.im_high_HDR;

% 如果预处理输出里带有 save_manifest，则打印出来用于核对来源。
if isfield(data_struct, 'save_manifest')
    if isfield(data_struct.save_manifest, 'outputs') && isfield(data_struct.save_manifest.outputs, 'output_tag')
        fprintf('preproc tag   : %s\n', data_struct.save_manifest.outputs.output_tag);
    elseif isfield(data_struct.save_manifest, 'output_tag')
        fprintf('preproc tag   : %s\n', data_struct.save_manifest.output_tag);
    end

    if isfield(data_struct.save_manifest, 'source') && isfield(data_struct.save_manifest.source, 'input_subdir')
        fprintf('preproc input : %s\n\n', data_struct.save_manifest.source.input_subdir);
    elseif isfield(data_struct.save_manifest, 'input_subdir')
        fprintf('preproc input : %s\n\n', data_struct.save_manifest.input_subdir);
    end
end

% 优先使用数据文件中的采样/光学参数，避免重建端与预处理端手动写两遍。
required_fields = {'spsize', 'wlength', 'z'};
for k = 1:numel(required_fields)
    field_name = required_fields{k};
    if ~isfield(data_struct, field_name)
        error('Input .mat missing required field %s: %s', field_name, data_dir);
    end
end

spsize = data_struct.spsize;
wlength = data_struct.wlength;
z = data_struct.z;

fprintf('%s size = %d x %d x %d\n', cfg.data_name, size(imlow_HDR));

%% ===== 4. Optional external LED intensity correction =====
% 物理含义：
% imlow_HDR 里存的是振幅 A，而光强标定给的是相对强度 I_rel。
% 因此校正时要使用：
%   A_corr = A_raw / sqrt(I_rel)
% if cfg.calibration.use_external_led_intensity_correction
%     calib_path = fullfile(code_root, cfg.result_root, 'Calibration', ...
%         cfg.calibration.calib_batch_folder, cfg.calibration.calib_file_name);
%
%     calib_data = load(calib_path, 'calib');
%     if ~isfield(calib_data, 'calib')
%         error('Calibration file does not contain calib struct: %s', calib_path);
%     end
%     if ~isfield(calib_data.calib, 'led_intensity_fit_na_norm')
%         error('Calibration file missing led_intensity_fit_na_norm: %s', calib_path);
%     end
%
%     led_intensity_rel = double(calib_data.calib.led_intensity_fit_na_norm(:));
%     if numel(led_intensity_rel) ~= size(imlow_HDR, 3)
%         error('LED calibration length (%d) does not match number of images (%d).', ...
%             numel(led_intensity_rel), size(imlow_HDR, 3));
%     end
%     if any(led_intensity_rel <= 0) || any(~isfinite(led_intensity_rel))
%         error('LED calibration contains non-positive or non-finite values.');
%     end
%
%     amp_scale = 1 ./ sqrt(led_intensity_rel);
%     for img_idx = 1:size(imlow_HDR, 3)
%         imlow_HDR(:, :, img_idx) = imlow_HDR(:, :, img_idx) .* amp_scale(img_idx);
%     end
% end

%% ===== 5. Debug display =====
show_raw_mode = 'center';  % 'center' 或 'all'

if strcmp(show_raw_mode, 'center')
    figure('Color', 'w');
    set(gcf, 'outerposition', get(0, 'ScreenSize'));
    imshow(imlow_HDR(:, :, 1), []);
    title('raw image 1');
elseif strcmp(show_raw_mode, 'all')
    for k = 1:size(imlow_HDR, 3)
        figure(1);
        set(gcf, 'outerposition', get(0, 'ScreenSize'));
        imshow(imlow_HDR(:, :, k), []);
        title(['raw image ' num2str(k)]);
        pause(0.1);
    end
end

%% ===== 6. Experiment geometry =====
% 这里用 cfg.geometry / cfg.recon 管理几何参数与算法参数，
% 便于以后切换不同物镜、LED 阵列或迭代策略。
aberration = 0;
theta = 1.0;
xint = 0;
yint = 0;

xstart = 0;
ystart = 0;
arraysize = cfg.geometry.arraysize;
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
[kx, ky, NAt] = k_vector(xlocation - xstart, ylocation - ystart, ...
    cfg.geometry.H, cfg.geometry.LEDp, cfg.geometry.nglass, cfg.geometry.t, ...
    theta, xint, yint, arraysize^2);

%% ===== 7. FPM reconstruction =====
NA = cfg.recon.NA;
upsmp_ratio = cfg.recon.upsmp_ratio;
psize = spsize / upsmp_ratio;

opts.loopnum = cfg.recon.loopnum;
opts.alpha = cfg.recon.alpha;
opts.beta = cfg.recon.beta;
opts.gamma_obj = cfg.recon.gamma_obj;
opts.gamma_p = cfg.recon.gamma_p;
opts.eta_obj = cfg.recon.eta_obj;
opts.eta_p = cfg.recon.eta_p;
opts.T = cfg.recon.T;
opts.aberration = aberration;
opts.use_internal_intensity_correction = cfg.calibration.use_external_led_intensity_correction;

% used_idx 是真正送入重建的 LED 编号集合。
% 例如 1:25 表示 25 张全用；1:2:25 表示隔一张取一张。
used_idx = cfg.recon.used_idx;
imlow_used = imlow_HDR(:, :, used_idx);
kx_used = kx(used_idx);
ky_used = ky(used_idx);

[him, tt, fprobe, imlow_HDR1] = himrecover( ...
    imlow_used, kx_used, ky_used, NA, wlength, spsize, psize, z, opts);

%% ===== 8. Result display and save =====
figure('Color', 'w');
set(gcf, 'outerposition', get(0, 'ScreenSize'));
subplot(121);
imshow(abs(him(50:end-50, 50:end-50)), []);
title('Amplitude High res');
subplot(122);
imshow(angle(him(50:end-50, 50:end-50)), []);
title('Phase High res');

preview_name = "Recover_" + cfg.result_tag + ".png";
saveas(gcf, fullfile(result_dir, preview_name));
result_mat_name = "Recover_" + cfg.result_tag + ".mat";

disp(['Wavelength: ', num2str(wlength * 1e+9), ' nm, Loop: ', num2str(opts.loopnum)]);
disp(['Maximum illumination NA = ', num2str(max(NAt(used_idx)))]);

% 与预处理端一样，重建结果也用同一套 manifest schema。
manifest_extra = struct();
manifest_extra.source.upstream_manifest = [];
if isfield(data_struct, 'save_manifest')
    manifest_extra.source.upstream_manifest = data_struct.save_manifest;
end
manifest_extra.recon.used_idx = used_idx(:);
manifest_extra.outputs.result_tag = cfg.result_tag;
manifest_extra.outputs.preview_png = preview_name;
manifest_extra.outputs.result_mat = result_mat_name;

recover_manifest = fpm_0619_build_manifest(cfg, "recover", manifest_extra);

save(fullfile(result_dir, result_mat_name), ...
    'him', 'fprobe', 'tt', 'imlow_HDR1', 'recover_manifest');

fprintf('Saved preview: %s\n', fullfile(result_dir, preview_name));
fprintf('Saved result : %s\n', fullfile(result_dir, result_mat_name));
