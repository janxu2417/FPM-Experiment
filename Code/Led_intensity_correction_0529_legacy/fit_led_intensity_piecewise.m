function result = fit_led_intensity_piecewise(led_intensity, varargin)
%FIT_LED_INTENSITY_PIECEWISE
% Piecewise fit of relative LED intensity for an odd-size NxN FPM LED set.
%
% This function assumes:
% 1) brightfield LEDs are a central odd-size core region,
% 2) darkfield LEDs are the remaining outer region,
% 3) the LED relative intensity is mainly a continuous function of the
%    radial distance from the array center.
%
% Input:
%   led_intensity : 25x1 or 1x25 vector, relative LED intensity after
%                   exposure correction.
%
% Name-value pairs:
%   'arraysize'      : LED array side length; [] means infer from data
%   'xstart'         : center LED x index, default 0
%   'ystart'         : center LED y index, default 0
%   'LEDp'           : LED pitch in mm, default 4
%   'core_size'      : side length of the central brightfield core, default 3
%   'exclude_idx'    : indices removed from fitting, default []
%   'bf_idx'         : manually specified brightfield indices, default []
%   'df_idx'         : manually specified darkfield indices, default []
%   'poly_order_bf'  : polynomial order for brightfield fit, default 2
%   'poly_order_df'  : polynomial order for darkfield fit, default 2
%   'normalize_mode' : 'mean1', 'center1', or 'none', default 'mean1'
%   'plot_result'    : true/false, default true
%
% Output:
%   result : structure containing geometry, fit coefficients, the spliced
%            fitted curve, and fitted intensity of all LEDs.

parser = inputParser;
parser.addRequired('led_intensity', @(x) isnumeric(x) && isvector(x));
parser.addParameter('arraysize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 3));
parser.addParameter('xstart', 0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ystart', 0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('LEDp', 4, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('core_size', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('exclude_idx', [], @(x) isnumeric(x) && isvector(x));
parser.addParameter('bf_idx', [], @(x) isnumeric(x) && isvector(x));
parser.addParameter('df_idx', [], @(x) isnumeric(x) && isvector(x));
parser.addParameter('poly_order_bf', 2, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('poly_order_df', 2, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('normalize_mode', 'mean1', @(x) ischar(x) || isstring(x));
parser.addParameter('plot_result', true, @(x) islogical(x) || isnumeric(x));
parser.parse(led_intensity, varargin{:});
params = parser.Results;

led_intensity = double(led_intensity(:));
num_led = numel(led_intensity);

if isempty(params.arraysize)
    arraysize = round(sqrt(num_led));
else
    arraysize = params.arraysize;
end

if arraysize^2 ~= num_led
    error('led_intensity length must equal arraysize^2.');
end
if mod(arraysize, 2) ~= 1
    error('arraysize must be odd for a centered LED array.');
end

core_size = params.core_size;
if core_size > arraysize || mod(core_size, 2) ~= 1
    error('core_size must be an odd integer and cannot exceed arraysize.');
end

addpath(genpath(fileparts(mfilename('fullpath'))));
[xlocation, ylocation] = LED_location(params.xstart, params.ystart, arraysize);
xlocation = xlocation(:);
ylocation = ylocation(:);

radius_led = hypot(xlocation - params.xstart, ylocation - params.ystart);
radius_mm = params.LEDp * radius_led;

all_idx = (1:num_led).';
half_core = (core_size - 1) / 2;
auto_bf_idx = find(abs(xlocation - params.xstart) <= half_core & ...
                   abs(ylocation - params.ystart) <= half_core);

manual_bf_idx = unique(params.bf_idx(:));
manual_df_idx = unique(params.df_idx(:));
manual_bf_idx = manual_bf_idx(manual_bf_idx >= 1 & manual_bf_idx <= num_led);
manual_df_idx = manual_df_idx(manual_df_idx >= 1 & manual_df_idx <= num_led);

if ~isempty(manual_bf_idx) || ~isempty(manual_df_idx)
    if isempty(manual_bf_idx) || isempty(manual_df_idx)
        error('bf_idx and df_idx should either both be empty or both be provided.');
    end
    inner_idx = manual_bf_idx;
    outer_idx = manual_df_idx;
else
    inner_idx = auto_bf_idx;
    outer_idx = setdiff(all_idx, inner_idx, 'stable');
end

exclude_idx = unique(params.exclude_idx(:));
exclude_idx = exclude_idx(exclude_idx >= 1 & exclude_idx <= num_led);

bf_idx = setdiff(inner_idx, exclude_idx, 'stable');
df_idx = setdiff(outer_idx, exclude_idx, 'stable');

if isempty(bf_idx) || isempty(df_idx)
    error('After exclusion, brightfield or darkfield data becomes empty.');
end

[rbf_u, ibf_mean, ibf_std, ibf_count] = local_group_by_radius(radius_mm(bf_idx), led_intensity(bf_idx));
[rdf_u, idf_mean, idf_std, idf_count] = local_group_by_radius(radius_mm(df_idx), led_intensity(df_idx));

order_bf = min(params.poly_order_bf, max(numel(rbf_u) - 1, 0));
order_df = min(params.poly_order_df, max(numel(rdf_u) - 1, 0));

pbf = polyfit(rbf_u, ibf_mean, order_bf);
pdf = polyfit(rdf_u, idf_mean, order_df);

bf_fun = @(r) polyval(pbf, r);
df_fun = @(r) polyval(pdf, r);

r_left = max(rbf_u);
r_right = min(rdf_u);

if r_right < r_left
    error('Brightfield and darkfield fitting domains overlap unexpectedly.');
end

r_grid = linspace(min(radius_mm), max(radius_mm), 400);
fit_grid = local_splice_curve(r_grid, bf_fun, df_fun, r_left, r_right);
fit_led = local_splice_curve(radius_mm, bf_fun, df_fun, r_left, r_right);

switch lower(string(params.normalize_mode))
    case "mean1"
        scale = mean(fit_led);
    case "center1"
        center_idx = find(radius_mm == 0, 1, 'first');
        scale = fit_led(center_idx);
    case "none"
        scale = 1;
    otherwise
        error('Unknown normalize_mode.');
end

if scale == 0
    error('Normalization scale is zero.');
end

fit_grid = fit_grid ./ scale;
fit_led = fit_led ./ scale;
led_intensity_norm = led_intensity ./ scale;
ibf_mean = ibf_mean ./ scale;
idf_mean = idf_mean ./ scale;
ibf_std = ibf_std ./ scale;
idf_std = idf_std ./ scale;

result = struct();
result.led_index = all_idx;
result.led_intensity_input = led_intensity;
result.led_intensity_input_norm = led_intensity_norm;
result.xlocation = xlocation;
result.ylocation = ylocation;
result.radius_led = radius_led;
result.radius_mm = radius_mm;
result.arraysize = arraysize;
result.core_size = core_size;
result.inner_idx = inner_idx;
result.outer_idx = outer_idx;
result.exclude_idx = exclude_idx;
result.bf_idx = bf_idx;
result.df_idx = df_idx;
result.poly_order_bf = order_bf;
result.poly_order_df = order_df;
result.poly_bf = pbf;
result.poly_df = pdf;
result.transition_radius_mm = [r_left, r_right];
result.rbf_unique_mm = rbf_u;
result.rdf_unique_mm = rdf_u;
result.ibf_mean = ibf_mean;
result.idf_mean = idf_mean;
result.ibf_std = ibf_std;
result.idf_std = idf_std;
result.ibf_count = ibf_count;
result.idf_count = idf_count;
result.fit_grid_radius_mm = r_grid;
result.fit_grid_intensity = fit_grid;
result.fit_led_intensity = fit_led;

if params.plot_result
    local_plot_result(result);
end

end

function [ru, imean, istd, icount] = local_group_by_radius(r, val)
[ru, ~, group_id] = uniquetol(r, 1e-10, 'DataScale', 1);
num_group = numel(ru);
imean = zeros(num_group, 1);
istd = zeros(num_group, 1);
icount = zeros(num_group, 1);

for k = 1:num_group
    data_k = val(group_id == k);
    imean(k) = mean(data_k);
    if numel(data_k) > 1
        istd(k) = std(data_k);
    else
        istd(k) = 0;
    end
    icount(k) = numel(data_k);
end

[ru, sort_idx] = sort(ru);
imean = imean(sort_idx);
istd = istd(sort_idx);
icount = icount(sort_idx);
end

function y = local_splice_curve(r, bf_fun, df_fun, r_left, r_right)
y = zeros(size(r));
mask_bf = r <= r_left;
mask_df = r >= r_right;
mask_mid = ~mask_bf & ~mask_df;

y(mask_bf) = bf_fun(r(mask_bf));
y(mask_df) = df_fun(r(mask_df));

if any(mask_mid)
    if abs(r_right - r_left) < eps
        y(mask_mid) = 0.5 * bf_fun(r(mask_mid)) + 0.5 * df_fun(r(mask_mid));
    else
        w = (r_right - r(mask_mid)) ./ (r_right - r_left);
        y(mask_mid) = w .* bf_fun(r(mask_mid)) + (1 - w) .* df_fun(r(mask_mid));
    end
end
end

function local_plot_result(result)
figure('Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on; grid on; box on;
plot(result.radius_mm(result.bf_idx), result.led_intensity_input_norm(result.bf_idx), ...
    'bo', 'MarkerFaceColor', [0.2 0.5 1], 'DisplayName', 'Brightfield data');
plot(result.radius_mm(result.df_idx), result.led_intensity_input_norm(result.df_idx), ...
    'rs', 'MarkerFaceColor', [1 0.4 0.4], 'DisplayName', 'Darkfield data');
if ~isempty(result.exclude_idx)
    plot(result.radius_mm(result.exclude_idx), result.led_intensity_input_norm(result.exclude_idx), ...
        'kx', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Excluded data');
end
errorbar(result.rbf_unique_mm, result.ibf_mean, result.ibf_std, ...
    'b-', 'LineWidth', 1.2, 'DisplayName', 'BF radial mean');
errorbar(result.rdf_unique_mm, result.idf_mean, result.idf_std, ...
    'r-', 'LineWidth', 1.2, 'DisplayName', 'DF radial mean');
plot(result.fit_grid_radius_mm, result.fit_grid_intensity, ...
    'k-', 'LineWidth', 1.8, 'DisplayName', 'Spliced fit');
xline(result.transition_radius_mm(1), '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xline(result.transition_radius_mm(2), '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xlabel('LED radial distance (mm)');
ylabel('Relative intensity');
title('LED intensity piecewise fit');
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
plot(result.led_index, result.led_intensity_input_norm, 'ko', ...
    'MarkerFaceColor', [0.7 0.7 0.7], 'DisplayName', 'Measured');
plot(result.led_index, result.fit_led_intensity, 'm-', ...
    'LineWidth', 1.8, 'DisplayName', 'Fitted');
xlabel('LED index');
ylabel('Relative intensity');
title('Measured vs fitted LED intensity');
legend('Location', 'best');
end
