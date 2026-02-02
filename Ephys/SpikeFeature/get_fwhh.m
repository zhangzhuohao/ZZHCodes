function [fwhh, trough_amplitude, peak_amplitude] = get_fwhh(t, v)
    % Input:
    % t - time vector (in ms or s)
    % v - voltage vector (in µV or mV)
    %
    % Output:
    % fwhh - full width at half height (in same units as t)
    % trough_amplitude - amplitude of negative trough
    % peak_amplitude - amplitude of positive peak

    % Find baseline (assuming it's near zero, taking mean of first 10% of data)
    baseline = mean(v(1:floor(length(v)*0.1)));
    
    % Find trough (minimum value)
    [trough_amplitude, trough_idx] = min(v);
    trough_amplitude = trough_amplitude - baseline;
    
    % Find peak after trough (maximum value after trough)
    [peak_amplitude, peak_idx] = max(v(trough_idx:end));
    peak_idx = peak_idx + trough_idx - 1;
    peak_amplitude = peak_amplitude - baseline;
    
    % Calculate half height level
    half_height = baseline + (trough_amplitude / 2);
    
    % Find points where waveform crosses half height
    % Before trough
    pre_crossings = find(diff(sign(v(1:trough_idx) - half_height)) ~= 0);
    if isempty(pre_crossings)
        disp('No crossing found before trough');
        fwhh = nan;
        return
    end
    t1_idx = pre_crossings(end);
    
    % After trough
    post_crossings = find(diff(sign(v(trough_idx:end) - half_height)) ~= 0);
    if isempty(post_crossings)
        disp('No crossing found after trough');
        fwhh = nan;
        return
    end
    t2_idx = trough_idx + post_crossings(1);
    
    % Interpolate to get more precise crossing times
    t1 = interp1(v(t1_idx:t1_idx+1), t(t1_idx:t1_idx+1), half_height);
    t2 = interp1(v(t2_idx-1:t2_idx), t(t2_idx-1:t2_idx), half_height);
    
    % Calculate full width at half height
    fwhh = t2 - t1;
end