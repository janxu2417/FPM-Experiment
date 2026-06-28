function fit_result = fit_led_intensity_smooth(led_intensity, fit_x, varargin)
%FIT_LED_INTENSITY_SMOOTH Smoothly fit LED relative intensity against a scalar geometry variable.
% Default fit_method is 'smoothing_spline'. If spline toolbox support is
% unavailable, the function falls back to smoothed pchip interpolation.

parser = inputParser;
parser.addRequired('led_intensity', @(x) isnumeric(x) && isvector(x));
parser.addRequired('fit_x', @(x) isnumeric(x) && isvector(x));
parser.addParameter('exclude_idx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
parser.addParameter('normalize_mode', 'mean1', @(x) ischar(x) || isstring(x));
parser.addParameter('fit_method', 'smoothing_spline', @(x) ischar(x) || isstring(x));
parser.addParameter('spline_p', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1));
parser.addParameter('label', 'radius_mm', @(x) ischar(x) || isstring(x));
parser.parse(led_intensity, fit_x, varargin{:});
params = parser.Results;

y = double(led_intensity(:));
x = double(fit_x(:));
num_led = numel(y);
if numel(x) ~= num_led
    error('fit_x length must match led_intensity length.');
end

exclude_idx = double(params.exclude_idx(:));
exclude_idx = exclude_idx(isfinite(exclude_idx));
exclude_idx = round(exclude_idx);
exclude_idx = unique(exclude_idx(exclude_idx >= 1 & exclude_idx <= num_led));
include_mask = true(num_led, 1);
include_mask(exclude_idx) = false;

x_use = x(include_mask);
y_use = y(include_mask);
[x_group, y_group, y_group_std, y_group_count] = local_group_stats(x_use, y_use);

fit_method = lower(string(params.fit_method));
[y_fit_group, y_fit_all, model_info] = local_fit_grouped_curve(x_group, y_group, x, fit_method, params.spline_p);

scale = local_get_scale(y_fit_all, x, params.normalize_mode);
if scale == 0
    error('Normalization scale is zero.');
end

y_norm = y ./ scale;
y_fit_all_norm = y_fit_all ./ scale;
y_group_norm = y_group ./ scale;
y_group_std_norm = y_group_std ./ scale;

residual = nan(num_led, 1);
residual(include_mask) = (y_norm(include_mask) - y_fit_all_norm(include_mask)) ./ y_fit_all_norm(include_mask);

fit_result = struct();
fit_result.label = char(params.label);
fit_result.fit_method = char(model_info.fit_method_used);
fit_result.model_info = model_info;
fit_result.exclude_idx = exclude_idx;
fit_result.include_mask = include_mask;
fit_result.x = x;
fit_result.y_measured = y_norm;
fit_result.y_fitted = y_fit_all_norm;
fit_result.residual = residual;
fit_result.group_x = x_group;
fit_result.group_y = y_group_norm;
fit_result.group_std = y_group_std_norm;
fit_result.group_count = y_group_count;
fit_result.same_x_std = y_group_std_norm;
fit_result.same_x_count = y_group_count;
fit_result.same_radius_std = y_group_std_norm;
fit_result.same_na_std = y_group_std_norm;
fit_result.normalization_mode = char(params.normalize_mode);
fit_result.normalization_scale = scale;

end

function [x_group, y_group, y_group_std, y_group_count] = local_group_stats(x, y)
[x_group, ~, group_id] = uniquetol(x, 1e-10, 'DataScale', 1);
num_group = numel(x_group);
y_group = zeros(num_group, 1);
y_group_std = zeros(num_group, 1);
y_group_count = zeros(num_group, 1);
for k = 1:num_group
    yk = y(group_id == k);
    y_group(k) = mean(yk);
    if numel(yk) > 1
        y_group_std(k) = std(yk);
    else
        y_group_std(k) = 0;
    end
    y_group_count(k) = numel(yk);
end
[x_group, sort_idx] = sort(x_group);
y_group = y_group(sort_idx);
y_group_std = y_group_std(sort_idx);
y_group_count = y_group_count(sort_idx);
end

function [y_fit_group, y_fit_all, model_info] = local_fit_grouped_curve(x_group, y_group, x_all, fit_method, spline_p)
model_info = struct();
model_info.fit_method_requested = char(fit_method);
model_info.fit_method_used = char(fit_method);
model_info.fallback_used = false;

if numel(x_group) == 1
    y_fit_group = y_group;
    y_fit_all = y_group(1) * ones(size(x_all));
    model_info.fit_method_used = 'constant';
    model_info.window = 1;
    return;
end

switch fit_method
    case "smoothing_spline"
        if exist('csaps', 'file') == 2
            if isempty(spline_p)
                p = max(0.7, min(0.98, 1 - 1 / max(numel(x_group), 3)));
            else
                p = double(spline_p);
            end
            pp = csaps(x_group.', y_group.', p);
            y_fit_group = fnval(pp, x_group.').';
            y_fit_all = fnval(pp, x_all.').';
            model_info.csaps_p = p;
            model_info.pp = pp;
        else
            y_fit_group = local_smoothed_pchip(x_group, y_group);
            y_fit_all = local_predict_from_anchor(x_group, y_fit_group, x_all, 'pchip');
            model_info.fit_method_used = 'smoothed_pchip_fallback';
            model_info.fallback_used = true;
        end
    case "pchip"
        y_fit_group = local_predict_from_anchor(x_group, y_group, x_group, 'pchip');
        y_fit_all = local_predict_from_anchor(x_group, y_group, x_all, 'pchip');
    case "poly2"
        order = min(2, numel(x_group) - 1);
        coef = polyfit(x_group, y_group, order);
        y_fit_group = polyval(coef, x_group);
        y_fit_all = polyval(coef, x_all);
        model_info.poly_order = order;
        model_info.poly_coef = coef;
    case "poly3"
        order = min(3, numel(x_group) - 1);
        coef = polyfit(x_group, y_group, order);
        y_fit_group = polyval(coef, x_group);
        y_fit_all = polyval(coef, x_all);
        model_info.poly_order = order;
        model_info.poly_coef = coef;
    case "poly4"
        order = min(4, numel(x_group) - 1);
        coef = polyfit(x_group, y_group, order);
        y_fit_group = polyval(coef, x_group);
        y_fit_all = polyval(coef, x_all);
        model_info.poly_order = order;
        model_info.poly_coef = coef;
    otherwise
        error('Unsupported fit_method: %s', fit_method);
end
end

function y_fit = local_smoothed_pchip(x_group, y_group)
if numel(y_group) <= 2
    y_fit = y_group;
    return;
end

window = min(5, numel(y_group));
if mod(window, 2) == 0
    window = window - 1;
end
kernel = ones(window, 1) / window;
ypad = [repmat(y_group(1), floor(window / 2), 1); y_group; repmat(y_group(end), floor(window / 2), 1)];
y_smooth = conv(ypad, kernel, 'valid');
y_fit = local_predict_from_anchor(x_group, y_smooth, x_group, 'pchip');
end

function yq = local_predict_from_anchor(x_anchor, y_anchor, xq, method)
if numel(x_anchor) == 1
    yq = y_anchor(1) * ones(size(xq));
else
    yq = interp1(x_anchor, y_anchor, xq, method, 'extrap');
end
end

function scale = local_get_scale(y_fit, x, normalize_mode)
switch lower(string(normalize_mode))
    case "mean1"
        scale = mean(y_fit);
    case "center1"
        center_idx = find(abs(x) == min(abs(x)), 1, 'first');
        scale = y_fit(center_idx);
    case "none"
        scale = 1;
    otherwise
        error('Unknown normalize_mode.');
end
end
