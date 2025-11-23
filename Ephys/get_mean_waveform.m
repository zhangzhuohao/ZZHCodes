function waveform_mean = get_mean_waveform(waveforms)

[n_channel, n_spike, n_sample] = size(waveforms);
waveform_mean = nan(n_channel, n_sample);

for i = 1:n_spike
    waveform_i = waveforms(i,:,:);
    waveform_i = reshape(waveform_i, n_spike, n_sample);
    waveform_i = rmoutliers(waveform_i, 'mean', 2);

    waveform_mean(i,:) = mean(waveform_i, 1, 'omitnan');
end
end