close all; clear; clc;

% Compare several 1D fit models using an existing diffuser calibration .mat.
% This script does not re-read raw images. It only re-fits the saved
% measured LED intensity curve against the saved geometry variable.

script_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(script_dir);
addpath(genpath(code_dir));

%% ===== User config =====
mat_path = fullfile(code_dir, 'Results', 'Calibration', '0619', ...
    'diffuser_calib_0619_光强修正_7x7_v2.mat');

% Compare against 'illum_na' or 'radius_mm'.
fit_domain = 'radius_mm';

% Use saved exclude_idx by default.
use_saved_exclude_idx = true;
manual_exclude_idx = [];

model_specs = { ...
    struct('name', 'spline p=0.98',  'fit_method', 'smoothing_spline', 'spline_p', 0.98), ...
    struct('name', 'spline p=0.995', 'fit_method', 'smoothing_spline', 'spline_p', 0.995), ...
    struct('name', 'spline p=0.999', 'fit_method', 'smoothing_spline', 'spline_p', 0.999), ...
    struct('name', 'poly3',          'fit_method', 'poly3',            'spline_p', []), ...
    struct('name', 'poly4',          'fit_method', 'poly4',            'spline_p', []) ...
};

%% ===== Load saved calibration =====
S = load(mat_path);
if ~isfield(S, 'calib')
    error('File does not contain variable "calib": %s', mat_path);
end
calib = S.calib;

y_measured = calib.measured.led_intensity_measured_norm(:);
switch lower(string(fit_domain))
    case "illum_na"
        fit_x = calib.geometry.illum_na(:);
        x_label = 'Illumination NA';
        output_tag = 'na';
    case "radius_mm"
        fit_x = calib.geometry.radius_mm(:);
        x_label = 'Radius (mm)';
        output_tag = 'radius';
    otherwise
        error('Unknown fit_domain: %s', fit_domain);
end

if use_saved_exclude_idx && isfield(calib.measured, 'exclude_idx')
    exclude_idx = calib.measured.exclude_idx(:).';
else
    exclude_idx = manual_exclude_idx;
end

num_model = numel(model_specs);
num_led = numel(y_measured);
fit_results = cell(num_model, 1);
fit_rmse_list = zeros(num_model, 1);
residual_rmse_list = zeros(num_model, 1);
max_abs_residual_list = zeros(num_model, 1);
mean_abs_residual_list = zeros(num_model, 1);

%% ===== Re-fit models =====
for k = 1:num_model
    spec = model_specs{k};
    fit_results{k} = fit_led_intensity_smooth(y_measured, fit_x, ...
        'exclude_idx', exclude_idx, ...
        'normalize_mode', 'none', ...
        'fit_method', spec.fit_method, ...
        'spline_p', spec.spline_p, ...
        'label', fit_domain);

    residual_k = fit_results{k}.residual(:);
    valid_mask = isfinite(residual_k);
    abs_residual_k = abs(residual_k(valid_mask));

    fit_rmse_list(k) = sqrt(mean((y_measured(valid_mask) - fit_results{k}.y_fitted(valid_mask)).^2));
    residual_rmse_list(k) = sqrt(mean(residual_k(valid_mask).^2));
    max_abs_residual_list(k) = max(abs_residual_k);
    mean_abs_residual_list(k) = mean(abs_residual_k);
end

%% ===== Print metrics =====
model_name = strings(num_model, 1);
for k = 1:num_model
    model_name(k) = string(model_specs{k}.name);
end

metrics_table = table(model_name, fit_rmse_list, residual_rmse_list, ...
    mean_abs_residual_list, max_abs_residual_list, ...
    'VariableNames', {'model', 'fit_rmse', 'residual_rmse', ...
    'mean_abs_residual', 'max_abs_residual'});
disp(metrics_table);

%% ===== Plot comparison =====
fig = figure('Color', 'w', 'Position', [120, 120, 1000, 650]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on; grid on; box on;
plot(1:num_led, y_measured, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.8 0.8 0.8]);
for k = 1:num_model
    plot(1:num_led, fit_results{k}.y_fitted, '-', 'LineWidth', 1.4);
end
xlabel('LED index');
ylabel('Relative intensity');
title(sprintf('Measured vs fitted (%s)', fit_domain), 'Interpreter', 'none');
legend_entries = [{'measured'}, cellfun(@(s) s.name, model_specs, 'UniformOutput', false)];
legend(legend_entries, 'Location', 'best');

nexttile;
hold on; grid on; box on;
for k = 1:num_model
    plot(1:num_led, fit_results{k}.residual, 'o-', 'LineWidth', 1.1);
end
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('LED index');
ylabel('Relative residual');
title(sprintf('Residual comparison (%s)', fit_domain), 'Interpreter', 'none');
legend(cellfun(@(s) s.name, model_specs, 'UniformOutput', false), 'Location', 'best');

nexttile;
bar(categorical(model_name), residual_rmse_list);
grid on; box on;
ylabel('Residual RMSE');
title('Residual RMSE by model');

nexttile;
yyaxis left
bar(categorical(model_name), max_abs_residual_list);
ylabel('Max |residual|');
yyaxis right
plot(categorical(model_name), mean_abs_residual_list, 'ko-', 'LineWidth', 1.3, 'MarkerFaceColor', 'k');
ylabel('Mean |residual|');
grid on; box on;
title('Residual summary by model');

sgtitle(sprintf('Diffuser fit model comparison: %s', fit_domain), 'Interpreter', 'none');

%% ===== Save figures =====
result_dir = fileparts(mat_path);
saveas(fig, fullfile(result_dir, sprintf('compare_fit_models_%s.png', output_tag)));

fprintf('\nSaved figure:\n%s\n', fullfile(result_dir, sprintf('compare_fit_models_%s.png', output_tag)));
