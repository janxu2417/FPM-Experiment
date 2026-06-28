%% FPM 模拟与实验对比 benchmark
clear; close all; clc;

sim_mat_path = 'FP_input_data_Sim.mat';
exp_mat_path = fullfile('Raw_input', '0619', 'FPM_RawData_z z=_0.050_7_preproc_v1.mat');

s = load(sim_mat_path, 'imseqlow1');
e = load(exp_mat_path);

imseqlow1 = s.imseqlow1;
if isfield(e, 'im_high_HDR')
    exp_stack = e.im_high_HDR;
elseif isfield(e, 'imlow_HDR')
    exp_stack = e.imlow_HDR;
else
    error('Experimental .mat does not contain im_high_HDR or imlow_HDR: %s', exp_mat_path);
end

num_LED = min(size(imseqlow1, 3), size(exp_stack, 3));
sim = imseqlow1(:, :, 1:num_LED);
exp = exp_stack(:, :, 1:num_LED);

sim = sim ./ max(sim(:, :, 1), [], 'all');
exp = exp ./ max(exp(:, :, 1), [], 'all');

energy_sim = squeeze(mean(mean(sim, 1), 2));
energy_exp = squeeze(mean(mean(exp, 1), 2));

idx_bright = 1;
idx_mid = 12;
idx_dark = 23;

figure('Name', 'Spatial comparison');
subplot(2, 3, 1); imagesc(sim(:, :, idx_bright)); axis image off; colormap gray; title('Sim brightfield');
subplot(2, 3, 2); imagesc(sim(:, :, idx_mid)); axis image off; colormap gray; title('Sim transition');
subplot(2, 3, 3); imagesc(sim(:, :, idx_dark)); axis image off; colormap gray; title('Sim darkfield');
subplot(2, 3, 4); imagesc(exp(:, :, idx_bright)); axis image off; colormap gray; title('Exp brightfield');
subplot(2, 3, 5); imagesc(exp(:, :, idx_mid)); axis image off; colormap gray; title('Exp transition');
subplot(2, 3, 6); imagesc(exp(:, :, idx_dark)); axis image off; colormap gray; title('Exp darkfield');

figure('Name', 'Line profile and energy');
subplot(1, 2, 1);
plot(sim(round(end / 2), :, idx_bright), 'b-', 'LineWidth', 1.2); hold on;
plot(exp(round(end / 2), :, idx_bright), 'r-', 'LineWidth', 1.2);
grid on; legend('Simulation', 'Experiment');
title('Center brightfield line profile');
xlabel('Pixel'); ylabel('Normalized amplitude');

subplot(1, 2, 2);
semilogy(energy_sim ./ max(energy_sim), 'b.-', 'LineWidth', 1.0); hold on;
semilogy(energy_exp ./ max(energy_exp), 'r.-', 'LineWidth', 1.0);
grid on; legend('Simulation', 'Experiment');
title('LED energy decay');
xlabel('LED index'); ylabel('Normalized energy');

result = FPM_resolution_benchmark(sim_mat_path, exp_mat_path);
disp(result.summary_table);
