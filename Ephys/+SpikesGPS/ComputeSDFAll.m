function [SDFAll, SDFTable] = ComputeSDFAll(r, id)
%
% figure out which unit it is
spk_note = r.Units.SpikeNotes;
switch length(id)
    case 1
        ind_unit = id;
        id = spk_note(ind_unit, [1 2]);
    case 2
        ind_unit = find(spk_note(:, 1)==id(1) & spk_note(:, 2)==id(2));
end
switch spk_note(ind_unit, 3)
    case 1
        spk_quality = 'SU';
    case 2
        spk_quality = 'MU';
end

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

%% Gather all
e_tbl = r.EphysTable;

% event times in original order
t_initin  = e_tbl.tInitIn;
t_initout = e_tbl.tInitOut;
t_centin  = e_tbl.tCentIn;
t_trigger = e_tbl.tTrigger;
t_centout = e_tbl.tCentOut;
t_choice  = e_tbl.tChoice;

% trial information
trial_id = e_tbl.Trials;
outcome  = e_tbl.Outcome;
FP       = e_tbl.FP;
port     = e_tbl.PortCorrect;
cued     = e_tbl.Cued;
FP(outcome=="Probe") = Inf;

n_trials = length(t_centin);

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

win_pre = zeros(n_trials, 1);
for i = 1:n_trials
    win_pre(i) = min([ST(i) + pre_*1000, pre_max*1000]);
    win_pre(i) = max([win_pre(i), pre_min*1000]);
end

win_post = t_choice - t_centin + post_*1000;
win_post(isnan(win_post)) = t_centout(isnan(win_post)) - t_centin(isnan(win_post)) + post_*1000;

t_spk_times = getSpikeTimingsWithin(spktimes, t_centin, [-win_pre win_post]);

t_sdf_trial = cell(n_trials, 1);
sdf_trial   = cell(n_trials, 1);

for i = 1:n_trials
    [sdf_trial{i}, t_sdf_trial{i}] = sdf26(t_spk_times{i}, [-win_pre(i) win_post(i)], 20, 1);
end

%% save whole trial sdf to strcut SDFOut
SDFAll.trial_id = trial_id;

SDFAll.t_initin  = t_initin;
SDFAll.t_initout = t_initout;
SDFAll.t_centin  = t_centin;
SDFAll.t_trigger = t_trigger;
SDFAll.t_centout = t_centout;
SDFAll.t_choice  = t_choice;

SDFAll.FP = FP;
SDFAll.Port = port;
SDFAll.Cued = cued;
SDFAll.Outcome = outcome;
SDFAll.ST = ST;
SDFAll.RT = RT;
SDFAll.HD = HD;
SDFAll.MT = MT;

SDFAll.t_spk_times = t_spk_times;
SDFAll.t_sdf = t_sdf_trial;
SDFAll.sdf = sdf_trial;

SDFTable = struct2table(SDFAll);

% basal information
this_unit = r.BehaviorClass.Subject + '|Session' + r.BehaviorClass.Session + '|Ch' + num2str(id(1)) + '|Unit' + num2str(id(2));
fprintf('\n*** %s ***\n', this_unit);

SDFAll.meta.unit_name = this_unit;
SDFAll.meta.subject   = r.BehaviorClass.Subject;
SDFAll.meta.session   = r.BehaviorClass.Session;
SDFAll.meta.Channel   = num2str(id(1));
SDFAll.meta.Unit      = num2str(id(2));
SDFAll.meta.Quality   = spk_quality;

SDFAll.spike_waveform.tspk     = tspk;
SDFAll.spike_waveform.spk_mean = spk_mean;
SDFAll.spike_waveform.spk_std  = spk_std;

end % ComputeSDFWarped
