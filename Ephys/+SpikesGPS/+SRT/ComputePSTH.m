function PSTH = ComputePSTH(r, PSTHOut, ku, varargin)

% Jianing Yu 5/8/2023
% For plotting PSTHs under SRT condition.
% Extracted from SRTSpikes

% Modified by Yue Huang on 6/26/2023
% Change the way of making raster plots to run faster

% 5/21/2025
% Revised from ComputePlotPSTH.m
% Remove plot part

% close all;
PSTH.UnitID = ku;
if nargin>2
    for i=1:2:size(varargin,2)
        switch varargin{i}
            %             case 'FRrange'
            %                 FRrange = varargin{i+1};
            case 'CentInTimeDomain'
                CentInTimeDomain = varargin{i+1}; % PSTH time domain
            case 'CentOutTimeDomain'
                CentOutTimeDomain = varargin{i+1}; % PSTH time domain
            case 'ChoiceTimeDomain'
                ChoiceTimeDomain = varargin{i+1};
            case 'TriggerTimeDomain'
                TriggerTimeDomain = varargin{i+1};
            case 'InitInTimeDomain'
                InitInTimeDomain = varargin{i+1};
            case 'InitOutTimeDomain'
                InitOutTimeDomain = varargin{i+1};
            otherwise
                errordlg('unknown argument')
        end
    end
end

FPs    = unique(PSTHOut.CentIn.FP{1});
NumFPs = length(FPs);

Ports    = unique(PSTHOut.CentIn.Port{1});
NumPorts = length(Ports);

%% PSTHs for CentIn and CentOut
params_cent_in.pre      = 5000; % take a longer pre-press activity so we can compute z score easily later.
params_cent_in.post     = CentInTimeDomain(2);
params_cent_in.binwidth = 20;

t_cent_in = PSTHOut.CentIn.Time{end};
[psth_cent_in_all, ts_cent_in_all, trialspxmat_cent_in_all, tspkmat_cent_in_all, t_correct_cent_in_all, ~] = jpsth(r.Units.SpikeTimes(ku).timings, t_cent_in, params_cent_in);

psth_cent_in_all = smoothdata(psth_cent_in_all, 'gaussian', 5);
PSTH.CentInAll   = {psth_cent_in_all, ts_cent_in_all, trialspxmat_cent_in_all, tspkmat_cent_in_all,  t_correct_cent_in_all};

% CentIn PSTH (corrected, sorted)
params_cent_in.pre      = CentInTimeDomain(1);
params_cent_in.post     = CentInTimeDomain(2);
params_cent_in.binwidth = 20;

t_cent_in_correct    = cell(NumFPs, NumPorts);
psth_cent_in_correct = cell(NumFPs, NumPorts);
ts_cent_in           = cell(NumFPs, NumPorts);
trialspxmat_cent_in  = cell(NumFPs, NumPorts);
tspkmat_cent_in      = cell(NumFPs, NumPorts);
rt_cent_in_sorted    = cell(NumFPs, NumPorts);

ind_correct = PSTHOut.CentIn.Labels=="Correct";
for i = 1:NumFPs
    for j = 1:NumPorts
        t_ij = PSTHOut.CentIn.Time{ind_correct}{i,j};
        [psth_cent_in_correct{i,j}, ts_cent_in{i,j}, trialspxmat_cent_in{i,j}, tspkmat_cent_in{i,j}, t_cent_in_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_cent_in);
        psth_cent_in_correct{i,j} = smoothdata(psth_cent_in_correct{i,j}, 'gaussian', 5);

        rt_cent_in_sorted{i,j} = PSTHOut.CentIn.RT_Correct{i,j};
        rt_cent_in_sorted{i,j} = rt_cent_in_sorted{i,j}(ind);

        PSTH.CentIn{i,j} = {psth_cent_in_correct{i,j}, ts_cent_in{i,j}, trialspxmat_cent_in{i,j}, tspkmat_cent_in{i,j}, t_cent_in_correct{i,j}, rt_cent_in_sorted{i,j}};
    end
end
PSTH.CentInLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT'};

% CentOut PSTH (corrected, sorted)
params_cent_out.pre      = CentOutTimeDomain(1);
params_cent_out.post     = CentOutTimeDomain(2);
params_cent_out.binwidth = 20;

t_cent_out_correct    = cell(NumFPs, NumPorts);
psth_cent_out_correct = cell(NumFPs, NumPorts);
ts_cent_out           = cell(NumFPs, NumPorts);
trialspxmat_cent_out  = cell(NumFPs, NumPorts);
tspkmat_cent_out      = cell(NumFPs, NumPorts);
rt_cent_out_sorted    = cell(NumFPs, NumPorts);

ind_correct = PSTHOut.CentOut.Labels=="Correct";
for i = 1:NumFPs
    for j = 1:NumPorts
        t_ij = PSTHOut.CentOut.Time{ind_correct}{i,j};
        [psth_cent_out_correct{i,j}, ts_cent_out{i,j}, trialspxmat_cent_out{i,j}, tspkmat_cent_out{i,j}, t_cent_out_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_cent_out);
        psth_cent_out_correct{i,j} = smoothdata(psth_cent_out_correct{i,j}, 'gaussian', 5);

        rt_cent_out_sorted{i,j} = PSTHOut.CentIn.RT_Correct{i,j};
        rt_cent_out_sorted{i,j} = rt_cent_out_sorted{i,j}(ind);

        PSTH.CentOut{i,j} = {psth_cent_out_correct{i,j}, ts_cent_out{i,j}, trialspxmat_cent_out{i,j}, tspkmat_cent_out{i,j}, t_cent_out_correct{i,j}, rt_cent_out_sorted{i,j}};
    end
end
PSTH.CentOutLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT'};

% premature cent_in PSTH
t_premature_cent_in           = cell(1, NumPorts);
psth_premature_cent_in        = cell(1, NumPorts);
ts_premature_cent_in          = cell(1, NumPorts);
trialspxmat_premature_cent_in = cell(1, NumPorts);
tspkmat_premature_cent_in     = cell(1, NumPorts);
hd_premature_cent_in          = cell(1, NumPorts);
FP_premature_cent_in          = cell(1, NumPorts);

ind_premature = PSTHOut.CentIn.Labels=="Premature";
for j = 1:NumPorts
    t_j = PSTHOut.CentIn.Time{ind_premature}{j};
    [psth_premature_cent_in{j}, ts_premature_cent_in{j}, trialspxmat_premature_cent_in{j}, tspkmat_premature_cent_in{j}, t_premature_cent_in{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_in);
    psth_premature_cent_in{j} = smoothdata(psth_premature_cent_in{j}, 'gaussian', 5);

    hd_premature_cent_in{j} = PSTHOut.CentIn.HoldDur.Premature{j};
    hd_premature_cent_in{j} = hd_premature_cent_in{j}(ind);

    FP_premature_cent_in{j} = PSTHOut.CentIn.FP{ind_premature}{j};
    FP_premature_cent_in{j} = FP_premature_cent_in{j}(ind);

    PSTH.PrematureCentIn{j} = {psth_premature_cent_in{j}, ts_premature_cent_in{j}, trialspxmat_premature_cent_in{j}, tspkmat_premature_cent_in{j}, t_premature_cent_in{j}, hd_premature_cent_in{j}, FP_premature_cent_in{j}};
end
PSTH.PrematureLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP'};

% premature cent_out PSTH
t_premature_cent_out           = cell(1, NumPorts);
psth_premature_cent_out        = cell(1, NumPorts);
ts_premature_cent_out          = cell(1, NumPorts);
trialspxmat_premature_cent_out = cell(1, NumPorts);
tspkmat_premature_cent_out     = cell(1, NumPorts);
hd_premature_cent_out          = cell(1, NumPorts);
FP_premature_cent_out          = cell(1, NumPorts);

ind_premature = PSTHOut.CentOut.Labels=="Premature";
for j = 1:NumPorts
    t_j = PSTHOut.CentOut.Time{ind_premature}{j};
    [psth_premature_cent_out{j}, ts_premature_cent_out{j}, trialspxmat_premature_cent_out{j}, tspkmat_premature_cent_out{j}, t_premature_cent_out{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_out);
    psth_premature_cent_out{j} = smoothdata(psth_premature_cent_out{j}, 'gaussian', 5);

    hd_premature_cent_out{j} = PSTHOut.CentOut.HoldDur.Premature{j};
    hd_premature_cent_out{j} = hd_premature_cent_out{j}(ind);

    FP_premature_cent_out{j} = PSTHOut.CentOut.FP{ind_premature}{j};
    FP_premature_cent_out{j} = FP_premature_cent_out{j}(ind);

    PSTH.PrematureCentOut{j} = {psth_premature_cent_out{j}, ts_premature_cent_out{j}, trialspxmat_premature_cent_out{j}, tspkmat_premature_cent_out{j}, t_premature_cent_out{j}, hd_premature_cent_out{j}, FP_premature_cent_out{j}};
end

% late cent_in PSTH
t_late_cent_in           = cell(1, NumPorts);
psth_late_cent_in        = cell(1, NumPorts);
ts_late_cent_in          = cell(1, NumPorts);
trialspxmat_late_cent_in = cell(1, NumPorts);
tspkmat_late_cent_in     = cell(1, NumPorts);
hd_late_cent_in          = cell(1, NumPorts);
FP_late_cent_in          = cell(1, NumPorts);

ind_late = PSTHOut.CentIn.Labels=="Late";
for j = 1:NumPorts
    t_j = PSTHOut.CentIn.Time{ind_late}{j};
    [psth_late_cent_in{j}, ts_late_cent_in{j}, trialspxmat_late_cent_in{j}, tspkmat_late_cent_in{j}, t_late_cent_in{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_in);
    psth_late_cent_in{j} = smoothdata(psth_late_cent_in{j}, 'gaussian', 5);

    hd_late_cent_in{j} = PSTHOut.CentIn.HoldDur.Late{j};
    hd_late_cent_in{j} = hd_late_cent_in{j}(ind);

    FP_late_cent_in{j} = PSTHOut.CentIn.FP{ind_late}{j};
    FP_late_cent_in{j} = FP_late_cent_in{j}(ind);

    PSTH.LateCentIn{j} = {psth_late_cent_in{j}, ts_late_cent_in{j}, trialspxmat_late_cent_in{j}, tspkmat_late_cent_in{j}, t_late_cent_in{j}, hd_late_cent_in{j}, FP_late_cent_in{j}};
end
PSTH.LateLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP'};

% late cent_out PSTH
t_late_cent_out           = cell(1, NumPorts);
psth_late_cent_out        = cell(1, NumPorts);
ts_late_cent_out          = cell(1, NumPorts);
trialspxmat_late_cent_out = cell(1, NumPorts);
tspkmat_late_cent_out     = cell(1, NumPorts);
hd_late_cent_out          = cell(1, NumPorts);
FP_late_cent_out          = cell(1, NumPorts);

ind_late = PSTHOut.CentOut.Labels=="Late";
for j = 1:NumPorts
    t_j = PSTHOut.CentOut.Time{ind_late}{j};
    [psth_late_cent_out{j}, ts_late_cent_out{j}, trialspxmat_late_cent_out{j}, tspkmat_late_cent_out{j}, t_late_cent_out{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_out);
    psth_late_cent_out{j} = smoothdata(psth_late_cent_out{j}, 'gaussian', 5);

    hd_late_cent_out{j} = PSTHOut.CentOut.HoldDur.Late{j};
    hd_late_cent_out{j} = hd_late_cent_out{j}(ind);

    FP_late_cent_out{j} = PSTHOut.CentOut.FP{ind_late}{j};
    FP_late_cent_out{j} = FP_late_cent_out{j}(ind);

    PSTH.LateCentOut{j} = {psth_late_cent_out{j}, ts_late_cent_out{j}, trialspxmat_late_cent_out{j}, tspkmat_late_cent_out{j}, t_late_cent_out{j}, hd_late_cent_out{j}, FP_late_cent_out{j}};
end

% probe cent_in PSTH
t_probe_cent_in           = cell(1, NumPorts);
psth_probe_cent_in        = cell(1, NumPorts);
ts_probe_cent_in          = cell(1, NumPorts);
trialspxmat_probe_cent_in = cell(1, NumPorts);
tspkmat_probe_cent_in     = cell(1, NumPorts);
hd_probe_cent_in          = cell(1, NumPorts);

ind_probe = PSTHOut.CentIn.Labels=="Probe";
for j = 1:NumPorts
    t_j = PSTHOut.CentIn.Time{ind_probe}{j};
    [psth_probe_cent_in{j}, ts_probe_cent_in{j}, trialspxmat_probe_cent_in{j}, tspkmat_probe_cent_in{j}, t_probe_cent_in{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_in);
    psth_probe_cent_in{j} = smoothdata(psth_probe_cent_in{j}, 'gaussian', 5);

    hd_probe_cent_in{j} = PSTHOut.CentIn.HoldDur.Probe{j};
    hd_probe_cent_in{j} = hd_probe_cent_in{j}(ind);

    PSTH.ProbeCentIn{j} = {psth_probe_cent_in{j}, ts_probe_cent_in{j}, trialspxmat_probe_cent_in{j}, tspkmat_probe_cent_in{j}, t_probe_cent_in{j}, hd_probe_cent_in{j}};
end

% late cent_out PSTH
t_probe_cent_out           = cell(1, NumPorts);
psth_probe_cent_out        = cell(1, NumPorts);
ts_probe_cent_out          = cell(1, NumPorts);
trialspxmat_probe_cent_out = cell(1, NumPorts);
tspkmat_probe_cent_out     = cell(1, NumPorts);
hd_probe_cent_out          = cell(1, NumPorts);

ind_probe = PSTHOut.CentOut.Labels=="Probe";
for j = 1:NumPorts
    t_j = PSTHOut.CentOut.Time{ind_probe}{j};
    [psth_probe_cent_out{j}, ts_probe_cent_out{j}, trialspxmat_probe_cent_out{j}, tspkmat_probe_cent_out{j}, t_probe_cent_out{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_out);
    psth_probe_cent_out{j} = smoothdata(psth_probe_cent_out{j}, 'gaussian', 5);

    hd_probe_cent_out{j} = PSTHOut.CentOut.HoldDur.Probe{j};
    hd_probe_cent_out{j} = hd_probe_cent_out{j}(ind);

    PSTH.ProbeCentOut{j} = {psth_probe_cent_out{j}, ts_probe_cent_out{j}, trialspxmat_probe_cent_out{j}, tspkmat_probe_cent_out{j}, t_probe_cent_out{j}, hd_probe_cent_out{j}};
end
PSTH.ProbeLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration'};

%% PSTH for choice poke
% use t_reward_poke and move_time to construct reward_poke PSTH
% reward PSTH
params_choice.pre  = ChoiceTimeDomain(1);
params_choice.post = ChoiceTimeDomain(2);
params_choice.binwidth = 20;

t_reward_choice           = cell(NumFPs, NumPorts);
psth_reward_choice        = cell(NumFPs, NumPorts);
ts_reward_choice          = cell(NumFPs, NumPorts);
trialspxmat_reward_choice = cell(NumFPs, NumPorts);
tspkmat_reward_choice     = cell(NumFPs, NumPorts);
mt_reward_choice          = cell(NumFPs, NumPorts);
for i = 1:NumFPs
    for j = 1:NumPorts
        t_ij = PSTHOut.ChoiceIn.RewardPoke.Time{i,j};
        [psth_reward_choice{i,j}, ts_reward_choice{i,j}, trialspxmat_reward_choice{i,j}, tspkmat_reward_choice{i,j}, t_reward_choice{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_choice);
        psth_reward_choice{i,j} = smoothdata(psth_reward_choice{i,j}, 'gaussian', 5);

        mt_reward_choice{i,j} = PSTHOut.ChoiceIn.RewardPoke.MoveTime{i,j}(ind);
        PSTH.RewardChoice{i,j} = {psth_reward_choice{i,j}, ts_reward_choice{i,j}, trialspxmat_reward_choice{i,j}, tspkmat_reward_choice{i,j}, t_reward_choice{i,j}, mt_reward_choice{i,j}};
    end
end
PSTH.ChoiceLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'MoveTime'};

% nonreward choice PSTH
t_nonreward_choice           = cell(1, NumPorts);
psth_nonreward_choice        = cell(1, NumPorts);
ts_nonreward_choice          = cell(1, NumPorts);
trialspxmat_nonreward_choice = cell(1, NumPorts);
tspkmat_nonreward_choice     = cell(1, NumPorts);
mt_nonreward_choice          = cell(1, NumPorts);
for j = 1:NumPorts
    t_j = PSTHOut.ChoiceIn.NonrewardPoke.Time{j};
    [psth_nonreward_choice{j}, ts_nonreward_choice{j}, trialspxmat_nonreward_choice{j}, tspkmat_nonreward_choice{j}, t_nonreward_choice{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_choice);
    psth_nonreward_choice{j} = smoothdata(psth_nonreward_choice{j}, 'gaussian', 5);

    mt_nonreward_choice{j}  = PSTHOut.ChoiceIn.NonrewardPoke.MoveTime{j}(ind);
    PSTH.NonrewardChoice{j} = {psth_nonreward_choice{j}, ts_nonreward_choice{j}, trialspxmat_nonreward_choice{j}, tspkmat_nonreward_choice{j}, t_nonreward_choice{j}, mt_nonreward_choice{j}};
end

%% PSTH for trigger
% trigger PSTH
params_trigger.pre  = TriggerTimeDomain(1);
params_trigger.post = TriggerTimeDomain(2);
params_trigger.binwidth = 20;

% correct response
t_trigger_correct           = cell(NumFPs, NumPorts);
psth_trigger_correct        = cell(NumFPs, NumPorts);
ts_trigger_correct          = cell(NumFPs, NumPorts);
trialspxmat_trigger_correct = cell(NumFPs, NumPorts);
tspkmat_trigger_correct     = cell(NumFPs, NumPorts);
rt_trigger_correct          = cell(NumFPs, NumPorts);

ind_correct = PSTHOut.Triggers.Labels=="Correct";
for i = 1:NumFPs
    for j = 1:NumPorts
        t_ij = PSTHOut.Triggers.Time{ind_correct}{i,j};
        [psth_trigger_correct{i,j}, ts_trigger_correct{i,j}, trialspxmat_trigger_correct{i,j}, tspkmat_trigger_correct{i,j}, t_trigger_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_trigger);
        psth_trigger_correct{i,j} = smoothdata(psth_trigger_correct{i,j}, 'gaussian', 5);

        rt_trigger_correct{i,j} = PSTHOut.Triggers.RT{ind_correct}{i,j};
        rt_trigger_correct{i,j} = rt_trigger_correct{i,j}(ind);
        PSTH.Triggers{i,j} = {psth_trigger_correct{i,j}, ts_trigger_correct{i,j}, trialspxmat_trigger_correct{i,j}, tspkmat_trigger_correct{i,j}, t_trigger_correct{i,j}, rt_trigger_correct{i,j}, FPs(i), Ports(j)};
    end
end
PSTH.TriggerLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT', 'FP', 'Port'};

% late response
t_trigger_late           = cell(1, NumPorts);
psth_trigger_late        = cell(1, NumPorts);
ts_trigger_late          = cell(1, NumPorts);
trialspxmat_trigger_late = cell(1, NumPorts);
tspkmat_trigger_late     = cell(1, NumPorts);
rt_trigger_late          = cell(1, NumPorts);
FP_trigger_late          = cell(1, NumPorts);

ind_late = PSTHOut.Triggers.Labels=="Late";
for j = 1:NumPorts
    t_j = PSTHOut.Triggers.Time{ind_late}{j};
    [psth_trigger_late{j}, ts_trigger_late{j}, trialspxmat_trigger_late{j}, tspkmat_trigger_late{j}, t_trigger_late{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_trigger);
    psth_trigger_late{j} = smoothdata(psth_trigger_late{j}, 'gaussian', 5);

    rt_trigger_late{j} = PSTHOut.Triggers.RT{ind_late}{j};
    rt_trigger_late{j} = rt_trigger_late{j}(ind);

    FP_trigger_late{j} = PSTHOut.Triggers.FP{ind_late}{j};
    FP_trigger_late{j} = FP_trigger_late{j}(ind);

    PSTH.TriggerLate{j} = {psth_trigger_late{j}, ts_trigger_late{j}, trialspxmat_trigger_late{j}, tspkmat_trigger_late{j}, t_trigger_late{j}, rt_trigger_late{j}, FP_trigger_late{j}, Ports(j)};
end

%% PSTH for init-in and init-out
% Init-In
params_initin.pre  = InitInTimeDomain(1);
params_initin.post = InitInTimeDomain(2);
params_initin.binwidth = 20;

% Init-In pre-correct
t_cor = PSTHOut.InitIn.PreCor.Time;
[psth_cor_initin, ts_cor_initin, trialspxmat_cor_initin, tspkmat_cor_initin, t_cor_initin, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_cor, params_initin);
psth_cor_initin = smoothdata(psth_cor_initin, 'gaussian', 5);

dur_cor_initin = PSTHOut.InitIn.PreCor.InitDur(ind);
PSTH.PreCorrectInitIn = {psth_cor_initin, ts_cor_initin, trialspxmat_cor_initin, tspkmat_cor_initin, t_cor_initin, dur_cor_initin};

% Init-In pre-error
t_err = PSTHOut.InitIn.PreErr.Time;
[psth_err_initin, ts_err_initin, trialspxmat_err_initin, tspkmat_err_initin, t_err_initin, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_err, params_initin);
psth_err_initin = smoothdata(psth_err_initin, 'gaussian', 5);

dur_err_initin = PSTHOut.InitIn.PreErr.InitDur(ind);
PSTH.PreErrorInitIn = {psth_err_initin, ts_err_initin, trialspxmat_err_initin, tspkmat_err_initin, t_err_initin, dur_err_initin};

% Init-Out
params_initout.pre  = InitOutTimeDomain(1);
params_initout.post = InitOutTimeDomain(2);
params_initout.binwidth = 20;

% Init-Out pre-correct
t_cor = PSTHOut.InitOut.PreCor.Time;
[psth_cor_initout, ts_cor_initout, trialspxmat_cor_initout, tspkmat_cor_initout, t_cor_initout, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_cor, params_initout);
psth_cor_initout = smoothdata(psth_cor_initout, 'gaussian', 5);

dur_cor_initout = PSTHOut.InitOut.PreCor.InitDur(ind);
PSTH.PreCorrectInitOut = {psth_cor_initout, ts_cor_initout, trialspxmat_cor_initout, tspkmat_cor_initout, t_cor_initout, dur_cor_initout};

% Init-Out pre-error
t_err = PSTHOut.InitOut.PreErr.Time;
[psth_err_initout, ts_err_initout, trialspxmat_err_initout, tspkmat_err_initout, t_err_initout, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_err, params_initout);
psth_err_initout = smoothdata(psth_err_initout, 'gaussian', 5);

dur_err_initout = PSTHOut.InitOut.PreErr.InitDur(ind);
PSTH.PreErrorInitOut = {psth_err_initout, ts_err_initout, trialspxmat_err_initout, tspkmat_err_initout, t_err_initout, dur_err_initout};

end