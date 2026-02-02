function [feature_out, wave_m] = get_waveform_features(t_wave, wave_all)
% waveform features (amp_trough, t_trough, amp_peak, t_peak, peak_trough_ratio, peak_trough_delay, fwhh)
[~, ~, tf_outlier] = rmoutliers(wave_all, 1);
wave_all(tf_outlier) = nan;
wave_m = mean(wave_all, 1, 'omitnan');

v_baseline = mean(wave_m(t_wave<0.25));
wave_all = wave_all - v_baseline;
wave_m   = wave_m - v_baseline;
if -min(wave_m) < max(wave_m) % regular spike, negative phase is larger than the  positive phase
    wave_all = - wave_all;
    wave_m   = - wave_m;
end

[amp_trough, ind_trough] = min(wave_all, [], 2);
amp_trough = mean(amp_trough, 'omitnan');
t_trough   = t_wave(ind_trough);
t_trough   = mean(t_trough, 'omitnan');

[amp_peak, ind_peak] = max(wave_all(:, t_wave>t_trough), [], 2);
amp_peak = mean(amp_peak, 'omitnan');
t_peak   = t_wave(ind_peak);
t_peak   = mean(t_peak, 'omitnan') + t_wave(find(t_wave<=t_trough, 1, 'last'));

peak_trough_ratio = abs(amp_peak / amp_trough);
peak_trough_delay = t_peak - t_trough;
fwhh = get_fwhh(t_wave, wave_m);

feature_out.amp_trough = amp_trough;
feature_out.t_trough   = t_trough;
feature_out.amp_peak = amp_peak;
feature_out.t_peak   = t_peak;
feature_out.peak_trough_ratio = peak_trough_ratio;
feature_out.peak_trough_delay = peak_trough_delay;
feature_out.fwhh = fwhh;

end
