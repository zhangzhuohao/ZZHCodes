function SDFOut = ComputeSDFWarped(r, id)
%
% figure out which unit it is
spk_note = r.Units.SpikeNotes;
ind_unit = find(spk_note(:, 1)==id(1) & spk_note(:, 2)==id(2));
spkwave  = r.Units.SpikeTimes(ind_unit).wave;
spktimes = r.Units.SpikeTimes(ind_unit).timings; % in ms
spktimes_ = zeros(1, ceil(max(spktimes)));
spktimes_(round(spktimes)) = 1;
[ar, lags]  = xcorr(spktimes_, 25);
ar(lags==0) = 0;

spks = spkwave;
tspk = (1:size(spks, 2))/30; % in ms
spk_mean = mean(spks, 1);
spk_std  = std(spks, 0, 1);

%% Task information
FPs    = r.BehaviorClass.TargetFP; % you have to use BuildR2023 or BuildR4Tetrodes2023 to have this included in r.
FPs    = FPs(FPs>0);
NumFPs = length(FPs);

Ports    = r.BehaviorClass.LeftRight;
NumPorts = length(Ports);

%% Gather all
rb = r.Behavior;

id_initin  = find(strcmp(rb.Labels, 'PokeInitIn'));
id_initout = find(strcmp(rb.Labels, 'PokeInitOut'));
id_centin  = find(strcmp(rb.Labels, 'PokeCentIn'));
id_trigger = find(strcmp(rb.Labels, 'Trigger'));
id_centout = find(strcmp(rb.Labels, 'PokeCentOut'));
id_choice  = find(strcmp(rb.Labels, 'PokeChoiceIn'));

% event times in original order
t_initin  = rb.EventTimings(rb.EventMarkers==id_initin);
t_initout = rb.EventTimings(rb.EventMarkers==id_initout);
t_centin  = rb.EventTimings(rb.EventMarkers==id_centin);
t_trigger = rb.EventTimings(rb.EventMarkers==id_trigger);
t_centout = rb.EventTimings(rb.EventMarkers==id_centout);
t_choice  = rb.EventTimings(rb.EventMarkers==id_choice);

% trial information
trial_id = rb.TrialID;
outcome  = rb.Outcome;
FP       = rb.Foreperiods;
port     = rb.PortCorrect;
FP(outcome=="Probe") = Inf;

% align t_trigger and t_choice to correspond t_cent_in
n_trials = length(t_centin);
t_choice_aligned  = nan(n_trials, 1);
t_trigger_aligned = nan(n_trials, 1);
for i = 1:n_trials
    % choice
    if i < n_trials
        % choice poke should occur between cent-out and next init-in
        i_choice = find(t_choice>t_centout(i) & t_choice<t_initin(i+1), 1, 'first');
        if ~isempty(i_choice)
            t_choice_aligned(i) = t_choice(i_choice);
        end
    else
        % in the last trial, choice poke should occur after cent-out
        i_choice = find(t_choice>t_centout(i), 1, 'first');
        if ~isempty(i_choice)
            t_choice_aligned(i) = t_choice(i_choice);
        end
    end

    % trigger should occur between cent-in and cent-out
    i_trigger = find(t_trigger>t_centin(i) & t_trigger<t_centout(i), 1, 'first');
    if ~isempty(i_trigger)
        t_trigger_aligned(i) = t_trigger(i_trigger);
    end
end
t_choice  = t_choice_aligned;
t_trigger = t_trigger_aligned;

% get behavior times (in ms)
ST = t_centin - t_initout;
RT = t_centout - t_trigger;
HD = t_centout - t_centin;
MT = t_choice - t_centout;

%% get sdf of the whole trial
% from .5 sec before init-out, to 2.5 sec after choice (2.5 sec after cent-out in case there was no choice)
% the whole duration of time warp
pre_  = .5; % 0.5 sec before init-out
pre_max = 10; % if the rat took very long time for shuttling, trim it to 10 sec
pre_min = 2.5; % if the rat was too fast, set it to 2.5 sec
post_ = 3; % take 3 sec after choice/cent-out at first

t_spk_times = cell(n_trials, 1);
t_sdf_trial = cell(n_trials, 1);
sdf_trial   = cell(n_trials, 1);

for i = 1:n_trials
    i_pre = min([ST(i)/1000 + pre_, pre_max]);
    i_pre = max([i_pre, pre_min]);
    i_pre = round(i_pre*1000);

    if ~isnan(t_choice(i))
        total_dur = i_pre + round(t_choice(i)-t_centin(i)) + post_*1000; % the unit is ms, [pre, from cent-in to choice, post]
        this_spk_train = spktimes(spktimes>=t_centout(i)-i_pre & spktimes<=t_choice(i)+post_*1000);
    else
        total_dur = i_pre + round(t_centout(i)-t_centin(i)) + post_*1000; % the unit is ms, [pre, from cent-in to cent-out, post]
        this_spk_train = spktimes(spktimes>=t_centout(i)-i_pre & spktimes<=t_centout(i)+post_*1000);
    end

    % map spike times to spike matrix, then calculate spike density function
    tsdf = (0:total_dur-1) - i_pre; % in ms
    spkmat = zeros(1, total_dur); % each entry represents for 1 ms
    if ~isempty(this_spk_train)
        % align to cent-in
        this_spk_train = this_spk_train - t_centin(i);
        [~, ind_spikes] = intersect(round(tsdf), round(this_spk_train));
        spkmat(ind_spikes) = 1;
    end
    % convert the spike train to sdf
    sdfout = sdf(tsdf/1000, spkmat, 20); % spkout = sdf(tspk, spkin, kernel_width)

    t_spk_times{i} = this_spk_train;
    t_sdf_trial{i} = tsdf;
    sdf_trial{i}   = sdfout';
end

%% save whole trial sdf to strcut SDFOut
SDFOut.trial_id = trial_id;

SDFOut.t_initin  = t_initin;
SDFOut.t_initout = t_initout;
SDFOut.t_centin  = t_centin;
SDFOut.t_trigger = t_trigger;
SDFOut.t_centout = t_centout;
SDFOut.t_choice  = t_choice;

SDFOut.FP = FP;
SDFOut.Port = port;
SDFOut.Outcome = outcome;
SDFOut.ST = ST;
SDFOut.RT = RT;
SDFOut.HD = HD;
SDFOut.MT = MT;

SDFOut.t_spk_times = t_spk_times;
SDFOut.t_sdf = t_sdf_trial;
SDFOut.sdf = sdf_trial;

SDFTable = struct2table(SDFOut);

% basal information
this_unit = r.BehaviorClass.Subject + '|Session' + r.BehaviorClass.Session + '|Ch' + num2str(id(1)) + '|Unit' + num2str(id(2));
fprintf('\n*** %s ***\n', this_unit);

SDFOut.meta.unit_name = this_unit;
SDFOut.meta.subject   = r.BehaviorClass.Subject;
SDFOut.meta.session   = r.BehaviorClass.Session;
SDFOut.meta.Channel   = num2str(id(1));
SDFOut.meta.Unit      = num2str(id(2));

SDFOut.spike_waveform.tspk     = tspk;
SDFOut.spike_waveform.spk_mean = spk_mean;
SDFOut.spike_waveform.spk_std  = spk_std;

%% Warp sdf for each condition
% warp correct trials only
SDFCorrect = SDFTable(SDFOut.Outcome=="Correct", :);

% the whole duration of time warp
pre_  = 2.5; % 2.5 sec before cent-in, about the lower boundary of shuttle time
post_ = 3; % take 3 sec after choice at first
post_keep = 2.5; % should be shorter than post_, to avoid length error

latency = 2.5; % max cent-out to choice latency (max movement time)

% trim whole trial sdf to time warp duration
n_correct = height(SDFCorrect);
ind_valid = false(n_correct, 1);
for i = 1:n_correct
    i_centin  = SDFCorrect.t_centin(i);
    i_centout = SDFCorrect.t_centout(i);
    i_choice  = SDFCorrect.t_choice(i);
    i_mt = i_choice - i_centout;
    if ~isnan(i_choice) && i_mt<latency*1000
        % if there is a choice poke after cent-out && there is a choice poke within a certain latency after cent-out
        ind_valid(i) = true;

        total_dur = pre_*1000 + round(i_choice-i_centin) + post_*1000; % the unit is ms, [pre, from cent-in to choice, post]
        this_spk_train = spktimes(spktimes>=i_centin-pre_*1000 & spktimes<=i_choice+post_*1000);

        % map spike times to spike matrix, then calculate spike density function
        tsdf = (0:total_dur-1) - pre_*1000; % in ms
        spkmat = zeros(1, total_dur); % each entry represents for 1 ms
        if ~isempty(this_spk_train)
            % align to cent-in
            this_spk_train = this_spk_train - i_centin;
            [~, ind_spikes] = intersect(round(tsdf), round(this_spk_train));
            spkmat(ind_spikes) = 1;
        end
        % convert the spike train to sdf
        sdfout = sdf(tsdf/1000, spkmat, 20); % spkout = sdf(tspk, spkin, kernel_width)

        SDFCorrect.t_spk_times(i) = {this_spk_train};
        SDFCorrect.t_sdf(i) = {tsdf};
        SDFCorrect.sdf(i) = {sdfout'};
    end
end
SDFCorrect = SDFCorrect(ind_valid, :);

% Warp for each condition
fprintf('\n** Performing warp for each condition **\n');
warped_time_points = cell(NumFPs, NumPorts);
ioc_sdf_warped = cell(NumFPs, NumPorts);
t_warped       = cell(NumFPs, NumPorts);
for i = 1:NumFPs
    i_fp = FPs(i);
    for j = 1:NumPorts
        j_port = Ports(j);

        id_ij  = SDFCorrect.FP==i_fp & SDFCorrect.Port==j_port;
        SDF_ij = SDFCorrect(id_ij, :);

        median_hold_duration = median(SDF_ij.t_centout - SDF_ij.t_centin); % median, in ms
        median_movement_time = median(SDF_ij.t_choice - SDF_ij.t_centout); % median, in ms

        % jt_template has the following structures: end of FP, median hold
        % duration, median movement time
        jt_template = [0 i_fp*1000 median_hold_duration median_movement_time+median_hold_duration]; % [cent-in, trigger, cent-out, choice] in ms
        warped_time_points{i, j} = jt_template;
        dt = 1; % 1 ms
        jt_target_time = (0:dt:median_movement_time+median_hold_duration);
        sprintf('Critical time points are %2.2f\n', jt_template(1), jt_template(2), jt_template(3), jt_template(4));

        t_warptarget_first  = jt_target_time(jt_target_time>=i_fp*1000 & jt_target_time<jt_template(3)); % from trigger to cent-out
        t_warptarget_second = jt_target_time(jt_target_time>=jt_template(3) & jt_target_time<jt_template(4)); % from cent-out to choice

        ioc_sdf_warped{i, j} = [];
        for k = 1:sum(id_ij) % for each trial
            % jt is [cent-in, cent-out, choice]
            k_centin  = SDF_ij.t_centin(k);
            k_centout = SDF_ij.t_centout(k);
            k_choice  = SDF_ij.t_choice(k);

            jt = SDF_ij(k, ["t_centin", "t_centout", "t_choice"]); % normalize so that the first point is time 0
            jt = table2array(jt);
            jt = jt - jt(1);

            tsdf = SDF_ij.t_sdf{k}; % this is the time of sdf, defined previously as: tspk = (0:total_dur-1)-pre_*1000; % in ms
            jsdf = SDF_ij.sdf{k}; % this is the sdf

            not_warped = jsdf(tsdf<i_fp*1000);

            towarp_first   = jsdf(tsdf>=i_fp*1000 & tsdf<jt(2)); % from trigger to cent-out
            t_towarp_first = tsdf(tsdf>=i_fp*1000 & tsdf<jt(2)); % from trigger to cent-out
            towarp_second   = jsdf(tsdf>=jt(2) & tsdf<jt(3)); % from cent-out to choice
            t_towarp_second = tsdf(tsdf>=jt(2) & tsdf<jt(3)); % from cent-out to choice

            not_warped2 = jsdf(tsdf>=jt(3));
            not_warped2 = not_warped2(1:post_keep*1000); % max 2.5 sec after choice poke

            dt = 1; % 1 ms
            if ~isempty(t_towarp_first) && ~isempty(t_towarp_second)
                sdf_warped_first  = SpikesGPS.SRT.warp_sdf(t_towarp_first, towarp_first, t_warptarget_first); % input is V, X, and duration to warp
                sdf_warped_second = SpikesGPS.SRT.warp_sdf(t_towarp_second, towarp_second, t_warptarget_second);
                new_sdf = [not_warped, sdf_warped_first, sdf_warped_second, not_warped2];
                ioc_sdf_warped{i, j} = [ioc_sdf_warped{i, j}; new_sdf];
            end
        end
        t_warped{i, j} = (-pre_*1000:dt:length(new_sdf)-pre_*1000-dt);
    end
end

% Calculate mean and CI of warped sdf
sdf_warped_ci = cell(NumFPs, NumPorts);
sdf_warped_mean = cell(NumFPs, NumPorts);
for i = 1:NumFPs
    for j = 1:NumPorts
        if size(ioc_sdf_warped{i,j}, 1)>10
            sdf_warped_ci{i,j} = bootci(1000, @mean, ioc_sdf_warped{i,j});
        end
        sdf_warped_mean{i,j} = mean(ioc_sdf_warped{i,j}, 1);
    end
end

SDFWarp.t_points = warped_time_points;
SDFWarp.t_points_mean = ["CentIn", "Trigger", "CentOut", "Choice"];
SDFWarp.pre_post_durations = [pre_, post_, post_keep];

SDFWarp.t_warp   = t_warped;
SDFWarp.sdf_warp = ioc_sdf_warped;
SDFWarp.sdf_mean = sdf_warped_mean;
SDFWarp.sdf_ci   = sdf_warped_ci;

SDFOut.SDFWarp = SDFWarp;

%% Another warp method
% Pool all FPs for each Port
% Approach to Cent-In: no change
% Cent-In to Cent-In + 500 ms: no change
% Trigger - 250 ms to Trigger: no change
% Trigger to Cent-Out: warp
% Cent-Out to Choice: warp

% the whole duration of time warp
% around cent-in
pre_centin  = 2.5; % 2.5 sec before cent-in, about the lower boundary of shuttle time
post_centin = 0.5; % 0.5 sec after cent-in

% around trigger
pre_trigger = 0.25; % take 0.25 sec after choice at first
post_choice = 3; % take 3 sec after choice at first
post_keep = 2.5; % should be shorter than post_, to avoid length error

latency = 2.5; % max cent-out to choice latency (max movement time)

% Gather spike times and sdfs
% initialize variables
SDFCorrect.t_spk_centin = SDFCorrect.t_spk_times;
SDFCorrect.t_sdf_centin = SDFCorrect.t_sdf;
SDFCorrect.sdf_centin   = SDFCorrect.sdf;
SDFCorrect.t_spk_trigger = SDFCorrect.t_spk_times;
SDFCorrect.t_sdf_trigger = SDFCorrect.t_sdf;
SDFCorrect.sdf_trigger   = SDFCorrect.sdf;
% trim whole trial sdf to time warp duration
n_correct = height(SDFCorrect);
for i = 1:n_correct
    i_centin  = SDFCorrect.t_centin(i);
    i_trigger = SDFCorrect.t_trigger(i);
    i_centout = SDFCorrect.t_centout(i);
    i_choice  = SDFCorrect.t_choice(i);
    i_mt = i_choice - i_centout;

%%%%% For Cent-In around
    total_dur = pre_centin*1000 + post_centin*1000; % the unit is ms, [pre, cent-in, post]
    this_spk_train = spktimes(spktimes>=i_centin-pre_centin*1000 & spktimes<=i_centin+post_centin*1000);

    % map spike times to spike matrix, then calculate spike density function
    tsdf = (0:total_dur-1) - pre_centin*1000; % in ms
    spkmat = zeros(1, total_dur); % each entry represents for 1 ms
    if ~isempty(this_spk_train)
        % align to cent-in
        this_spk_train = this_spk_train - i_centin;
        [~, ind_spikes] = intersect(round(tsdf), round(this_spk_train));
        spkmat(ind_spikes) = 1;
    end
    % convert the spike train to sdf
    sdfout = sdf(tsdf/1000, spkmat, 20); % spkout = sdf(tspk, spkin, kernel_width)

    SDFCorrect.t_spk_centin(i) = {this_spk_train};
    SDFCorrect.t_sdf_centin(i) = {tsdf};
    SDFCorrect.sdf_centin(i) = {sdfout'};

%%%%% For Trigger around
    total_dur = pre_trigger*1000 + round(i_choice-i_trigger) + post_choice*1000; % the unit is ms, [pre, cent-in, post]
    this_spk_train = spktimes(spktimes>=i_trigger-pre_trigger*1000 & spktimes<=i_choice+post_choice*1000);

    % map spike times to spike matrix, then calculate spike density function
    tsdf = (0:total_dur-1) - pre_trigger*1000; % in ms
    spkmat = zeros(1, total_dur); % each entry represents for 1 ms
    if ~isempty(this_spk_train)
        % align to cent-in
        this_spk_train = this_spk_train - i_trigger;
        [~, ind_spikes] = intersect(round(tsdf), round(this_spk_train));
        spkmat(ind_spikes) = 1;
    end
    % convert the spike train to sdf
    sdfout = sdf(tsdf/1000, spkmat, 20); % spkout = sdf(tspk, spkin, kernel_width)

    SDFCorrect.t_spk_trigger(i) = {this_spk_train};
    SDFCorrect.t_sdf_trigger(i) = {tsdf};
    SDFCorrect.sdf_trigger(i) = {sdfout'};
end

% Warp for trigger around
fprintf('\n** Performing warp for pooled foreperiods **\n');
warped_time_points = cell(1, NumPorts);
ioc_sdf_warped = cell(1, NumPorts);
t_warped       = cell(1, NumPorts);
for j = 1:NumPorts
    j_port = Ports(j);
    id_j  = SDFCorrect.Port==j_port;
    SDF_j = SDFCorrect(id_j, :);

    median_reaction_time = median(SDF_j.t_centout - SDF_j.t_trigger); % median, in ms
    median_movement_time = median(SDF_j.t_choice - SDF_j.t_centout); % median, in ms

    % jt_template has the following structures: end of FP, median hold
    % duration, median movement time
    jt_template = [0 median_reaction_time median_movement_time+median_reaction_time]; % [trigger, cent-out, choice] in ms
    warped_time_points{1,j} = jt_template;
    dt = 1; % 1 ms
    jt_target_time = (0:dt:median_movement_time+median_reaction_time);
    sprintf('Critical time points are %2.2f\n', jt_template(1), jt_template(2), jt_template(3));

    t_warptarget_first  = jt_target_time(jt_target_time>=0 & jt_target_time<jt_template(2)); % from trigger to cent-out
    t_warptarget_second = jt_target_time(jt_target_time>=jt_template(2) & jt_target_time<jt_template(3)); % from cent-out to choice

    ioc_sdf_warped{1,j} = cell(1,2);
    t_warped{1,j} = cell(1,2);
    for k = 1:sum(id_j) % for each trial
        % jt is [cent-in, cent-out, choice]
        k_centin  = SDF_j.t_centin(k);
        k_centout = SDF_j.t_centout(k);
        k_choice  = SDF_j.t_choice(k);

        jt = SDF_j(k, ["t_trigger", "t_centout", "t_choice"]); % normalize so that the first point is time 0
        jt = table2array(jt);
        jt = jt - jt(1);

        tsdf = SDF_j.t_sdf_trigger{k}; % this is the time of sdf, defined previously as: tsdf = (0:total_dur-1)-pre_*1000; % in ms
        jsdf = SDF_j.sdf_trigger{k}; % this is the sdf

        not_warped = jsdf(tsdf<0);

        towarp_first   = jsdf(tsdf>=0 & tsdf<jt(2)); % from trigger to cent-out
        t_towarp_first = tsdf(tsdf>=0 & tsdf<jt(2)); % from trigger to cent-out
        towarp_second   = jsdf(tsdf>=jt(2) & tsdf<jt(3)); % from cent-out to choice
        t_towarp_second = tsdf(tsdf>=jt(2) & tsdf<jt(3)); % from cent-out to choice

        not_warped2 = jsdf(tsdf>=jt(3));
        not_warped2 = not_warped2(1:post_keep*1000); % max 2.5 sec after choice poke

        dt = 1; % 1 ms
        if ~isempty(t_towarp_first) && ~isempty(t_towarp_second)
            sdf_warped_first  = SpikesGPS.SRT.warp_sdf(t_towarp_first, towarp_first, t_warptarget_first); % input is V, X, and duration to warp
            sdf_warped_second = SpikesGPS.SRT.warp_sdf(t_towarp_second, towarp_second, t_warptarget_second);
            new_sdf = [not_warped, sdf_warped_first, sdf_warped_second, not_warped2];
            ioc_sdf_warped{1,j}{2} = [ioc_sdf_warped{1,j}{2}; new_sdf];

            % dont warp cent-in around
            ioc_sdf_warped{1,j}{1} = [ioc_sdf_warped{1,j}{1}; SDF_j.sdf_centin{k}];
        end
    end
    t_warped{1,j}{2} = (-pre_trigger*1000:dt:length(new_sdf)-pre_trigger*1000-dt);
    t_warped{1,j}{1} = SDF_j.t_sdf_centin{1};
end

% Calculate mean and CI of warped sdf
sdf_warped_ci = cell(1, NumPorts);
sdf_warped_mean = cell(1, NumPorts);
for j = 1:NumPorts
    for k = 1:2
        if size(ioc_sdf_warped{j}{k}, 1)>10
            sdf_warped_ci{j}{k} = bootci(1000, @mean, ioc_sdf_warped{j}{k});
        end
        sdf_warped_mean{j}{k} = mean(ioc_sdf_warped{j}{k}, 1);
    end
end

SDFWarpPool.t_points = warped_time_points;
SDFWarpPool.t_points_mean = ["Trigger", "CentOut", "Choice"];
SDFWarpPool.pre_post_durations = [pre_, post_, post_keep];

SDFWarpPool.t_warp   = t_warped;
SDFWarpPool.sdf_warp = ioc_sdf_warped;
SDFWarpPool.sdf_mean = sdf_warped_mean;
SDFWarpPool.sdf_ci   = sdf_warped_ci;

SDFOut.SDFWarpPool = SDFWarpPool;

SDFOut.SDFCorrect = table2struct(SDFCorrect, 'ToScalar', true);

end % ComputeSDFWarped
