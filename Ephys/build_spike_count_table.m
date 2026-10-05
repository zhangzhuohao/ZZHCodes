function spk_table = build_spike_count_table( ...
        spike_times_ms, trials, event_names, event_columns, ...
        count_range_ms, bin_size_ms)
%BUILD_SPIKE_COUNT_TABLE Build compact SRT trial-event-bin spike counts.
%
% Inputs
%   spike_times_ms : Complete sorted spike times for one unit, in ms.
%   trials         : NWB trials table. Absolute event times are in seconds.
%   event_names    : Output event labels, for example "cent_in".
%   event_columns  : Corresponding trials columns, for example "t_cent_in".
%   count_range_ms : Struct containing a [start end] range for each event.
%   bin_size_ms    : Spike-count bin width in ms.
%
% Output
%   spk_table      : One row per trial, event, and relative-time bin.

assert(istable(trials), 'trials must be a table.');
spike_times_ms = double(spike_times_ms(:));
assert(all(isfinite(spike_times_ms)) ...
    && all(diff(spike_times_ms) >= 0), ...
    'spike_times_ms must be finite and sorted.');
assert(numel(event_names) == numel(event_columns), ...
    'event_names and event_columns must have equal lengths.');
assert(all(ismember(event_columns, ...
    string(trials.Properties.VariableNames))), ...
    'trials is missing one or more requested event columns.');
assert(isscalar(bin_size_ms) && isfinite(bin_size_ms) ...
    && bin_size_ms > 0, 'bin_size_ms must be positive.');

n_trials = height(trials);
n_events = numel(event_names);
assert(n_trials <= intmax('uint16'), ...
    'Trial count exceeds uint16 storage range.');
n_rows = 0;
event_edges = cell(n_events, 1);
for i_event = 1:n_events
    event = char(event_names(i_event));
    assert(isfield(count_range_ms, event), ...
        'count_range_ms is missing event %s.', event);
    event_range_ms = count_range_ms.(event);
    assert(numel(event_range_ms) == 2 ...
        && all(isfinite(event_range_ms)) ...
        && event_range_ms(1) < event_range_ms(2), ...
        'Event range must be a finite [start end] vector: %s.', event);
    edges = event_range_ms(1):bin_size_ms:event_range_ms(2);
    if edges(end) < event_range_ms(2)
        edges(end + 1) = event_range_ms(2);
    end
    event_edges{i_event} = edges;
    valid_events = isfinite(trials.(char(event_columns(i_event))));
    n_rows = n_rows + nnz(valid_events) * (numel(edges) - 1);
end

trial_row = zeros(n_rows, 1, 'uint16');
event_label = strings(n_rows, 1);
event_time_s = nan(n_rows, 1);
bin_start_ms = nan(n_rows, 1);
bin_end_ms = nan(n_rows, 1);
spike_count = zeros(n_rows, 1, 'uint16');

next_row = 1;
for i_trial = 1:n_trials
    for i_event = 1:n_events
        this_event_s = double(trials.(char(event_columns(i_event)))(i_trial));
        if ~isfinite(this_event_s)
            continue
        end

        edges_rel_ms = event_edges{i_event};
        n_bins = numel(edges_rel_ms) - 1;
        rows = next_row:(next_row + n_bins - 1);
        edges_abs_ms = 1000 * this_event_s + edges_rel_ms;
        in_window = spike_times_ms >= edges_abs_ms(1) ...
            & spike_times_ms <= edges_abs_ms(end);
        counts = histcounts(spike_times_ms(in_window), edges_abs_ms);
        assert(all(counts <= intmax('uint16')), ...
            'A spike-count bin exceeds uint16 storage range.');

        trial_row(rows) = uint16(i_trial);
        event_label(rows) = event_names(i_event);
        event_time_s(rows) = this_event_s;
        bin_start_ms(rows) = edges_rel_ms(1:end-1);
        bin_end_ms(rows) = edges_rel_ms(2:end);
        spike_count(rows) = uint16(counts);
        next_row = next_row + n_bins;
    end
end
assert(next_row == n_rows + 1, ...
    'Spike-count table preallocation did not match populated rows.');

spk_table = table( ...
    trials.id(double(trial_row)), trial_row, categorical(event_label), ...
    event_time_s, bin_start_ms, bin_end_ms, ...
    (bin_start_ms + bin_end_ms) / 2, ...
    event_time_s + bin_start_ms / 1000, ...
    event_time_s + bin_end_ms / 1000, spike_count, ...
    'VariableNames', { ...
    'trial_id', 'trial_row', 'event', 'event_time_s', ...
    'bin_start_ms', 'bin_end_ms', 'bin_center_ms', ...
    'bin_start_time_s', 'bin_end_time_s', 'spike_count'});

covariates = ["outcome", "foreperiod", "stage", "cued", ...
    "port_correct", "port_lateral", "port_chosen"];
assert(all(ismember(covariates, ...
    string(trials.Properties.VariableNames))), ...
    'trials is missing one or more required spike-count covariates.');
for covariate = covariates
    name = char(covariate);
    value = trials.(name)(double(trial_row));
    if covariate == "outcome"
        value = categorical(string(value));
    end
    spk_table.(name) = value;
end
end
