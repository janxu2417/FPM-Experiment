close all; clear; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));

% Replace this vector with your measured LED relative intensity.
% For example, you can paste led_intensity_avg_norm from ROI calibration.
load("Raw_input\0529\LED_intensity_roi_0529.mat")

% led_intensity_avg_norm = [
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
%     1.0000
% ];

% If arraysize is left empty, it is inferred from numel(led_intensity_avg_norm).
arraysize = 5;

% Brightfield is assumed to be the central odd-size square core.
% For 5x5 use 3; for 7x7 you may still use 3, or increase to 5 if your
% brightfield region is experimentally reliable over the larger core.
core_size = 3;

% If some LEDs near the brightfield-darkfield boundary should be excluded
% from fitting, fill their indices here.
exclude_idx = [11,15,19,23,24];

result = fit_led_intensity_piecewise(led_intensity_avg_norm, ...
    'arraysize', arraysize, ...
    'xstart', 0, ...
    'ystart', 0, ...
    'LEDp', 4, ...
    'core_size', core_size, ...
    'exclude_idx', exclude_idx, ...
    'poly_order_bf', 2, ...
    'poly_order_df', 2, ...
    'normalize_mode', 'mean1', ...
    'plot_result', true);

disp('Fitted LED relative intensity:');
disp(result.fit_led_intensity(:).');
