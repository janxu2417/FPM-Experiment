function cfg = fpm_0619_shared_config(mode, preset_name)
%FPM_0619_SHARED_CONFIG 0619 FPM 实验的共享配置入口
% 用法：
%   cfg = fpm_0619_shared_config("preproc", "0619_zneg0050_main");
%   cfg = fpm_0619_shared_config("recover", "0619_zneg0050_preproc_v1");
%
% 设计目标：
% 1. 把 batch_folder、光学参数、几何参数、标定文件名等共用信息集中管理
% 2. 让预处理和重建脚本只保留主流程，减少重复函数和重复硬编码
% 3. 如果后续新增一套 0619 数据，只在这里补 preset

base = local_base_config();

switch char(mode)
    case 'preproc'
        cfg = local_preproc_config(base, preset_name);
    case 'recover'
        cfg = local_recover_config(base, preset_name);
    otherwise
        error('Unknown config mode: %s', mode);
end
end

function base = local_base_config()
base = struct();

% ===== project / path =====
base.batch_folder = "0619";
base.raw_root = "Raw_input";
base.result_root = "Results";

% ===== optical / sampling =====
base.optics = struct( ...
    'spsize', 0.69e-6, ...
    'm', 600, ...
    'n', 600, ...
    'dark_noise_level', 0, ...
    'wlength', 6.28e-07, ...
    'NA_obj', 0.15, ...
    'z', 0);

% ===== preproc shared defaults =====
base.preproc = struct( ...
    'num_img', 49, ...
    'file_indices', 1:49, ...
    'output_subdir', "", ...
    'base_roi', struct( ...
        'row_start_factor', -6 / 6, ...
        'row_end_factor', 0 / 6, ...
        'col_start_factor', -1 / 2, ...
        'col_end_factor', 1 / 2), ...
    'exposure_ms', local_exposure_49_zneg0050_7());
    % 'exposure_ms', local_exposure_first25_high());

% ===== calibration =====
base.calibration = struct( ...
    'use_external_led_intensity_correction', true, ...
    'calib_batch_folder', "0619", ...
    'calib_file_name', "diffuser_calib_0619_光强修正_7x7_poly4.mat");
    % 'calib_file_name', "diffuser_calib_0605_毛玻璃01_7x7.mat");

% ===== experiment geometry =====
base.geometry = struct( ...
    'arraysize', 7, ...
    'H', 40.24, ...
    'LEDp', 4, ...
    'nglass', 1.52, ...
    't', 1);

% ===== reconstruction defaults =====
base.recon = struct( ...
    'NA', 0.15, ...
    'upsmp_ratio', 4, ...
    'loopnum', 10, ...
    'alpha', 1, ...
    'beta', 1, ...
    'gamma_obj', 1, ...
    'gamma_p', 1, ...
    'eta_obj', 0.2, ...
    'eta_p', 0.2, ...
    'T', 1, ...
    'used_idx', 1:(base.geometry.arraysize ^ 2));
end

function cfg = local_preproc_config(base, preset_name)
cfg = base;
cfg.mode = "preproc";

switch char(preset_name)
    case 'template_new_batch_main'
        % 新一批数据模板：
        % 1. 先在 local_base_config 里按本批数据修改：
        %    base.batch_folder, base.optics, base.preproc, base.geometry
        % 2. 再复制这个 case，改 preset_name / input_subdir / output_tag / roi_configs
        % 3. 最后在 read_transfer_06xx.m 里切 active_preset
        %
        % 必改项至少包括：
        %   base.batch_folder        -> Raw_input\06xx
        %   cfg.input_subdir         -> 子文件夹 z=xxx_xx
        %   base.optics.m, n, spsize -> 图像分辨率/采样参数
        %   base.geometry.arraysize  -> LED 阵列边长
        cfg.preset_name = "template_new_batch_main";
        cfg.input_subdir = "z=example";
        cfg.output_tag = "z=example_preproc_v1";
        cfg.roi_configs = struct('ptr', {""}, 'row_shift', {0}, 'col_shift', {0});

    case '0619_zneg0050_main'
        cfg.preset_name = "0619_zneg0050_main";
        cfg.input_subdir = "z=_0.050";
        cfg.output_tag = "z=_0.050_preproc_v1";
        cfg.roi_configs = struct('ptr', {""}, 'row_shift', {0}, 'col_shift', {0});

    case '0619_zneg0050_up'
        cfg.preset_name = "0619_zneg0050_up";
        cfg.input_subdir = "z=_0.050";
        cfg.output_tag = "z=_0.050_preproc_v1";
        cfg.roi_configs = struct('ptr', {"_up"}, 'row_shift', {-base.optics.m}, 'col_shift', {0});

    case '0619_zneg0050_high2nd_multi'
        cfg.preset_name = "0619_zneg0050_high2nd_multi";
        cfg.input_subdir = "z=_0.050_high_2nd";
        cfg.output_tag = "z=_0.050_high_2nd_preproc_v1";
        cfg.roi_configs = struct( ...
            'ptr', {"_up", "", "_down"}, ...
            'row_shift', {-base.optics.m, 0, base.optics.m}, ...
            'col_shift', {0, 0, 0});

    case '0619_zneg0050_7_main'
        % 具体例子：对应目录 Code\Raw_input\0619\z=_0.050_7
        % 当前观察到文件名格式为 图像_1.tif, 图像_2.tif, ...
        % 如果这一批最终不是 25 张，而是 7x7 = 49 张，
        % 需要同时修改 local_base_config 里的：
        %   base.preproc.num_img
        %   base.preproc.file_indices
        %   base.geometry.arraysize
        cfg.preset_name = "0619_zneg0050_7_main";
        cfg.input_subdir = "z=_0.050_7";
        cfg.output_tag = "z=_0.050_7_preproc_v1";
        cfg.roi_configs = struct('ptr', {""}, 'row_shift', {0}, 'col_shift', {0});

    otherwise
        error('Unknown preproc preset: %s', preset_name);
end

cfg.output_names = local_build_preproc_output_names(cfg.output_tag, cfg.roi_configs);
end

function cfg = local_recover_config(base, preset_name)
cfg = base;
cfg.mode = "recover";

switch char(preset_name)
    case 'template_new_batch_recover'
        % 与 template_new_batch_main 对应的重建模板。
        % 通常需要改：
        %   cfg.preproc_output_tag
        %   cfg.preproc_ptr
        %   cfg.result_tag
        % 如果 LED 阵列边长改了，也要在 local_base_config 里改 base.geometry.arraysize。
        cfg.preset_name = "template_new_batch_recover";
        cfg.preproc_output_tag = "z=example_preproc_v1";
        cfg.preproc_ptr = "";
        cfg.calibration.use_external_led_intensity_correction = true;
        cfg.result_tag = "z=example_preproc_v1";

    case '0619_zneg0050_preproc_v1'
        cfg.preset_name = "0619_zneg0050_preproc_v1";
        cfg.preproc_output_tag = "z=_0.050_preproc_v1";
        cfg.preproc_ptr = "";
        cfg.calibration.use_external_led_intensity_correction = true;
        cfg.result_tag = "z=_0.050_preproc_v1";

    case '0619_zneg0050_preproc_v1_no_calib'
        cfg.preset_name = "0619_zneg0050_preproc_v1_no_calib";
        cfg.preproc_output_tag = "z=_0.050_preproc_v1";
        cfg.preproc_ptr = "";
        cfg.calibration.use_external_led_intensity_correction = false;
        cfg.result_tag = "z=_0.050_preproc_v1_no_calib";

    case '0619_zneg0050_preproc_v1_up'
        cfg.preset_name = "0619_zneg0050_preproc_v1_up";
        cfg.preproc_output_tag = "z=_0.050_preproc_v1";
        cfg.preproc_ptr = "_up";
        cfg.calibration.use_external_led_intensity_correction = true;
        cfg.result_tag = "z=_0.050_preproc_v1_up";

    case '0619_zneg0050_7_preproc_v3'
        % 使用poly4拟合led_intensity
        % 对应预处理输出：
        %   FPM_RawData_z z=_0.050_7_preproc_v1.mat
        cfg.preset_name = "0619_zneg0050_7_preproc_v1";
        cfg.preproc_output_tag = "z=_0.050_7_preproc_v1";
        cfg.preproc_ptr = "";
        cfg.calibration.use_external_led_intensity_correction = true;
        cfg.result_tag = "z=_0.050_7_preproc_v3";

    case '0619_zneg0050_7_led_corr_v2'
        % 对应预处理输出：
        %   FPM_RawData_z z=_0.050_7_preproc_v1.mat
        cfg.preset_name = "0619_zneg0050_7_preproc_v1";
        cfg.preproc_output_tag = "z=_0.050_7_preproc_v1";
        cfg.preproc_ptr = "";
        cfg.calibration.use_external_led_intensity_correction = false;
        cfg.result_tag = "z=_0.050_7_led_corr_v2";

    otherwise
        error('Unknown recover preset: %s', preset_name);
end

cfg.data_name = local_build_preproc_output_name(cfg.preproc_output_tag, cfg.preproc_ptr);
end

function exposure_ms = local_exposure_first25_high()
exposure_ms = 100.65 * ones(1, 25);
exposure_ms(1:9) = 42.307;
exposure_ms([13, 17, 21, 25]) = 201.301;
end

function exposure_ms = local_exposure_49_zneg0050_7()
exposure_ms = zeros(1, 49);

exposure_ms(1:9) = 44.189;
exposure_ms([10:12, 14:16, 18:20, 22:24]) = 88.378;
exposure_ms([13, 17, 21, 25, 27:29, 33:35, 39:41, 45:47]) = 210.218;
exposure_ms([26, 30:32, 36:38, 42:44, 48:49]) = 499.995;

end

function output_names = local_build_preproc_output_names(output_tag, roi_configs)
num_roi = numel(roi_configs);
output_names = cell(1, num_roi);
for k = 1:num_roi
    output_names{k} = local_build_preproc_output_name(output_tag, roi_configs(k).ptr);
end
end

function out_name = local_build_preproc_output_name(output_tag, ptr)
out_name = "FPM_RawData_z " + output_tag + ptr + ".mat";
end
