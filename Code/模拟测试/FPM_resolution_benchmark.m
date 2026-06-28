function result = FPM_resolution_benchmark(sim_mat_path, exp_mat_path)
% Quantify center brightfield resolution on the synthetic line-pair target.

if nargin < 1 || isempty(sim_mat_path)
    sim_mat_path = 'FP_input_data_Sim.mat';
end
if nargin < 2 || isempty(exp_mat_path)
    exp_mat_path = 'FPM_RawData_3HDR.mat';
end

s = load(sim_mat_path, 'imseqlow1', 'spsize', 'sim_cfg');
sim_img = s.imseqlow1(:, :, 1);
sim_img = sim_img ./ max(sim_img(:));
sim_metrics = fpm_measure_resolution_groups(sim_img, s.spsize);

has_exp = exist(exp_mat_path, 'file') == 2;
if has_exp
    e = load(exp_mat_path);
    if isfield(e, 'im_high_HDR')
        exp_stack = e.im_high_HDR;
    elseif isfield(e, 'imlow_HDR')
        exp_stack = e.imlow_HDR;
    else
        error('Experimental .mat does not contain im_high_HDR or imlow_HDR: %s', exp_mat_path);
    end
    exp_img = exp_stack(:, :, 1);
    exp_img = exp_img ./ max(exp_img(:));
    exp_metrics = fpm_measure_resolution_groups(exp_img, s.spsize);
else
    exp_metrics = [];
end

freqs = [sim_metrics.freq_lpmm].';
sim_avg = [sim_metrics.contrast_avg].';
sim_max = [sim_metrics.contrast_max].';

if has_exp
    exp_avg = [exp_metrics.contrast_avg].';
    exp_max = [exp_metrics.contrast_max].';
else
    exp_avg = nan(size(sim_avg));
    exp_max = nan(size(sim_max));
end

threshold = 0.03;
sim_resolved = freqs(sim_max >= threshold);
exp_resolved = freqs(exp_max >= threshold);

result = struct();
result.threshold = threshold;
result.sim_metrics = sim_metrics;
result.exp_metrics = exp_metrics;
result.sim_limit_lpmm = local_pick_limit(sim_resolved);
result.exp_limit_lpmm = local_pick_limit(exp_resolved);
result.summary_table = table(freqs, sim_avg, sim_max, exp_avg, exp_max, ...
    'VariableNames', {'freq_lpmm', 'sim_avg', 'sim_max', 'exp_avg', 'exp_max'});

fprintf('Resolution benchmark threshold = %.3f\n', threshold);
fprintf('Simulation resolved up to %.0f lp/mm\n', result.sim_limit_lpmm);
if has_exp
    fprintf('Experiment resolved up to %.0f lp/mm\n', result.exp_limit_lpmm);
end
end

function out = local_pick_limit(resolved_freqs)
if isempty(resolved_freqs)
    out = NaN;
    return;
end
out = min(resolved_freqs);
end
