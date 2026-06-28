%% ===== Analyze centroid_x / centroid_y vs kx / ky / LED index =====
clear; clc; close all;

% 1. Load calibration result
calib_file = fullfile( ...
    'D:\1_徐前\物理\傅里叶叠层显微术\Code\Results\Calibration\0605', ...
    'diffuser_calib_0605_毛玻璃01_5x5.mat');
S = load(calib_file, 'calib');
calib = S.calib;

% 2. Read variables
centroid_x = calib.qc.image_centroid_xy(:, 1);
centroid_y = calib.qc.image_centroid_xy(:, 2);
led_idx = (1:numel(centroid_x)).';

kx = calib.geometry.kx(:);
ky = calib.geometry.ky(:);

% Remove DC offset, easier to compare phase
cx0 = centroid_x - mean(centroid_x);
cy0 = centroid_y - mean(centroid_y);
kx0 = kx - mean(kx);
ky0 = ky - mean(ky);

% Normalize to unit std for shape comparison
cxn = cx0 / std(cx0);
cyn = cy0 / std(cy0);
kxn = kx0 / std(kx0);
kyn = ky0 / std(ky0);

% 3. Correlation check
corr_cx_kx = local_corrcoef(cxn, kxn);
corr_cx_ky = local_corrcoef(cxn, kyn);
corr_cy_kx = local_corrcoef(cyn, kxn);
corr_cy_ky = local_corrcoef(cyn, kyn);

fprintf('corr(centroid_x, kx) = %.4f\n', corr_cx_kx);
fprintf('corr(centroid_x, ky) = %.4f\n', corr_cx_ky);
fprintf('corr(centroid_y, kx) = %.4f\n', corr_cy_kx);
fprintf('corr(centroid_y, ky) = %.4f\n', corr_cy_ky);

% 4. Plot by LED index
figure('Color', 'w', 'Position', [100, 100, 1400, 900]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on; grid on; box on;
plot(led_idx, cxn, 'r-o', 'LineWidth', 1.2);
plot(led_idx, kxn, 'k--', 'LineWidth', 1.2);
xlabel('LED index');
ylabel('Normalized value');
title('centroid_x vs kx');
legend('centroid_x', 'kx', 'Location', 'best');

nexttile;
hold on; grid on; box on;
plot(led_idx, cyn, 'b-o', 'LineWidth', 1.2);
plot(led_idx, kyn, 'k--', 'LineWidth', 1.2);
xlabel('LED index');
ylabel('Normalized value');
title('centroid_y vs ky');
legend('centroid_y', 'ky', 'Location', 'best');

nexttile;
hold on; grid on; box on;
plot(led_idx, cxn, 'r-o', 'LineWidth', 1.2);
plot(led_idx, cyn, 'b-s', 'LineWidth', 1.2);
xlabel('LED index');
ylabel('Normalized value');
title('centroid_x and centroid_y');
legend('centroid_x', 'centroid_y', 'Location', 'best');

% 5. Scatter against kx / ky
nexttile;
hold on; grid on; box on;
scatter(kx, centroid_x, 50, led_idx, 'filled');
xlabel('kx');
ylabel('centroid_x');
title(sprintf('centroid_x vs kx, corr = %.3f', corr_cx_kx));
colorbar;

nexttile;
hold on; grid on; box on;
scatter(ky, centroid_y, 50, led_idx, 'filled');
xlabel('ky');
ylabel('centroid_y');
title(sprintf('centroid_y vs ky, corr = %.3f', corr_cy_ky));
colorbar;

% 6. 2D centroid trajectory
nexttile;
hold on; grid on; box on; axis equal;
plot(cx0, cy0, 'k-o', 'LineWidth', 1.2, 'MarkerFaceColor', [0.8 0.8 0.8]);
text(cx0, cy0, string(led_idx), 'FontSize', 8, 'Color', 'r');
xlabel('centroid_x - mean');
ylabel('centroid_y - mean');
title('2D centroid trajectory');

sgtitle('Centroid / k-vector comparison');

% 7. Estimate phase difference by FFT first harmonic
X = fft(cxn);
Y = fft(cyn);

% Ignore DC term, search strongest nonzero harmonic
[~, idx_fx] = max(abs(X(2:floor(end/2))));
[~, idx_fy] = max(abs(Y(2:floor(end/2))));
idx_fx = idx_fx + 1;
idx_fy = idx_fy + 1;

phase_x = angle(X(idx_fx));
phase_y = angle(Y(idx_fy));
phase_diff = angle(exp(1j * (phase_y - phase_x)));

fprintf('dominant harmonic index for centroid_x = %d\n', idx_fx);
fprintf('dominant harmonic index for centroid_y = %d\n', idx_fy);
fprintf('phase_x = %.4f rad\n', phase_x);
fprintf('phase_y = %.4f rad\n', phase_y);
fprintf('phase difference (y - x) = %.4f rad\n', phase_diff);
fprintf('phase difference / pi = %.4f pi\n', phase_diff / pi);

function r = local_corrcoef(a, b)
a = double(a(:));
b = double(b(:));

valid = isfinite(a) & isfinite(b);
a = a(valid);
b = b(valid);

a = a - mean(a);
b = b - mean(b);

den = sqrt(sum(a.^2) * sum(b.^2));
if den == 0
    r = NaN;
else
    r = sum(a .* b) / den;
end
end