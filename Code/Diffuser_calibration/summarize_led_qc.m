function qc = summarize_led_qc(max_dn, sat_ratio, hist_tail_ratio, main_roi_mean, main_roi_std, aux_roi_mean, aux_roi_std, thresholds)
%SUMMARIZE_LED_QC Build QC flags and summary statistics for all LEDs.

num_led = numel(max_dn);
aux_mean_spread = max(aux_roi_mean, [], 1) - min(aux_roi_mean, [], 1);
aux_std_mean = mean(aux_roi_std, 1);
main_to_aux_delta = abs(main_roi_mean - mean(aux_roi_mean, 1));

is_saturated = sat_ratio > thresholds.sat_ratio_warn;
is_tail_heavy = hist_tail_ratio > thresholds.hist_tail_warn;
is_roi_inconsistent = main_to_aux_delta > thresholds.roi_delta_warn;

flag_text = strings(num_led, 1);
for k = 1:num_led
    flags = strings(0, 1);
    if is_saturated(k)
        flags(end + 1) = "saturated";
    end
    if is_tail_heavy(k)
        flags(end + 1) = "tail_heavy";
    end
    if is_roi_inconsistent(k)
        flags(end + 1) = "roi_inconsistent";
    end
    if isempty(flags)
        flag_text(k) = "ok";
    else
        flag_text(k) = strjoin(flags, ",");
    end
end

qc = struct();
qc.max_dn = max_dn(:);
qc.sat_ratio = sat_ratio(:);
qc.hist_tail_ratio = hist_tail_ratio(:);
qc.main_roi_mean = main_roi_mean(:);
qc.main_roi_std = main_roi_std(:);
qc.aux_roi_mean = aux_roi_mean;
qc.aux_roi_std = aux_roi_std;
qc.aux_mean_spread = aux_mean_spread(:);
qc.aux_std_mean = aux_std_mean(:);
qc.main_to_aux_delta = main_to_aux_delta(:);
qc.thresholds = thresholds;
qc.is_saturated = is_saturated(:);
qc.is_tail_heavy = is_tail_heavy(:);
qc.is_roi_inconsistent = is_roi_inconsistent(:);
qc.flag_text = flag_text;
end
