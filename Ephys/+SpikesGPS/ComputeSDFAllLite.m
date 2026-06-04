function SDFAll = ComputeSDFAllLite(spktimes, e_tbl, sigma, dt)
% 
if nargin<3
    sigma = 20;
    dt    = 1;
elseif nargin<4
    dt = 1;
end

%% Gather behavior information
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

win_ext = [-50 50]; % ms
win_pre = zeros(n_trials, 1);
for i = 1:n_trials
    win_pre(i) = min([ST(i) + pre_*1000, pre_max*1000]);
    win_pre(i) = max([win_pre(i), pre_min*1000]);
end

win_post = t_choice - t_centin + post_*1000;
win_post(isnan(win_post)) = median(HD, 'omitnan') + median(MT, 'omitnan') + post_*1000;

t_spk_times = getSpikeTimingsWithin(spktimes, t_centin, [-win_pre win_post] + win_ext);

t_sdf_trial = cell(n_trials, 1);
sdf_trial   = cell(n_trials, 1);

for i = 1:n_trials
    win_i = [-win_pre(i) win_post(i)];
    [sdf_i, t_sdf_i] = sdf26(t_spk_times{i}, win_i + win_ext, sigma, dt);
    ind_i = t_sdf_i>=win_i(1) & t_sdf_i<win_i(2);
    sdf_trial{i}   = sdf_i(ind_i);
    t_sdf_trial{i} = t_sdf_i(ind_i);
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
SDFAll.sdf   = sdf_trial;

end % ComputeSDFWarpedLite