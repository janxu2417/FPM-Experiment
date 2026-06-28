function [RMSE] = similarity(lim,him)
%RMSE quantifies how close the reconstruction quality of the image is to the true value
%   此处显示详细说明

% 计算幅度比误差和相位比误差
amp_ratio = abs(lim) ./ abs(him);
phase_ratio = angle(lim) ./ angle(him);

% 避免除零和 NaN
amp_ratio(isinf(amp_ratio) | isnan(amp_ratio)) = 0;
phase_ratio(isinf(phase_ratio) | isnan(phase_ratio)) = 0;

% 计算 RMSE
RMSE_amp = sqrt(mean((amp_ratio - 1).^2, 'all'));   % 幅度 RMSE
RMSE_phase = sqrt(mean((phase_ratio - 1).^2, 'all')); % 相位 RMSE

RMSE = [RMSE_amp, RMSE_phase];

end