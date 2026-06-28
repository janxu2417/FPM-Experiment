function sim_cfg = fpm_sim_nonideal_config()
% Shared non-ideal settings for FPM simulation under current 0.15 NA setup.

% sim_cfg.partial_coherence：LED 扩展光源导致的部分相干退化
sim_cfg.partial_coherence.enable = true;
sim_cfg.partial_coherence.source_size_mm = 0.4;
sim_cfg.partial_coherence.grid_n = 3;
sim_cfg.partial_coherence.gaussian_weight_sigma_mm = 0.18;
% 像面上重复模糊
sim_cfg.partial_coherence.enable_image_blur = false;
sim_cfg.partial_coherence.image_sigma_px = 0.3;
sim_cfg.partial_coherence.filter_size_px = 5;

% sim_cfg.pupil_mtf：物镜 pupil/MTF 的高频滚降
sim_cfg.pupil_mtf.enable = true;
sim_cfg.pupil_mtf.rolloff_sigma = 1.20;
sim_cfg.pupil_mtf.rolloff_power = 1.4;
end
