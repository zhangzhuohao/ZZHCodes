function PSTH = ComputePlotPSTH(r, PSTHOut, ku, varargin)

% Jianing Yu 5/8/2023
% For plotting PSTHs under SRT condition.
% Extracted from SRTSpikes

% Modified by Yue Huang on 6/26/2023
% Change the way of making raster plots to run faster

% close all;
PSTH.UnitID = ku;
ToSave = 'on';
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
            case 'ToSave'
                ToSave = varargin{i+1};
            otherwise
                errordlg('unknown argument')
        end
    end
end

c = GPSColor();
p_c = .9;
c_port = {p_c*c.Contra+(1-p_c)*[1 1 1], p_c*c.Ipsi+(1-p_c)*[1 1 1]};

% For PSTH and raster plots
c_cent_in  = [5 191 219] / 255;
c_trigger  = [247 182 45] / 255;
c_cent_out = [238 5 219] / 255;
c_reward   = [164 208 164] / 255;

TargetFPs = unique(PSTHOut.CentIn.FP{1});
nFPs = length(TargetFPs);
if nFPs == 2
    c_FP = [192, 127, 0; 76, 61, 61]/255;
    lw = [1.5 1];
    ls = [":", "-"];
else
    c_FP = [255, 217, 90; 192, 127, 0; 76, 61, 61]/255;
    lw = [1.5 1.25 1];
    ls = [":", "-.", "-"];
end

Ports = unique(PSTHOut.CentIn.Port{1});
nPorts = length(Ports);

nSort = nFPs * nPorts;

c_premature = [0.9 0.4 0.1];
c_late      = [0.6 0.6 0.6];
c_probe     = [.7 .5 .3];
printsize   = [2 1 17 25];

%% PSTHs for CentIn and CentOut
params_cent_in.pre      = 5000; % take a longer pre-press activity so we can compute z score easily later.
params_cent_in.post     = CentInTimeDomain(2);
params_cent_in.binwidth = 20;

t_cent_in = PSTHOut.CentIn.Time{end};
[psth_cent_in_all, ts_cent_in_all, trialspxmat_cent_in_all, tspkmat_cent_in_all, t_correct_cent_in_all, ~] = jpsth(r.Units.SpikeTimes(ku).timings, t_cent_in, params_cent_in);

psth_cent_in_all = smoothdata (psth_cent_in_all, 'gaussian', 5);
PSTH.CentInAll   = {psth_cent_in_all, ts_cent_in_all, trialspxmat_cent_in_all, tspkmat_cent_in_all,  t_correct_cent_in_all};

% CentIn PSTH (corrected, sorted)
params_cent_in.pre      = CentInTimeDomain(1);
params_cent_in.post     = CentInTimeDomain(2);
params_cent_in.binwidth = 20;

t_cent_in_correct    = cell(nFPs, nPorts);
psth_cent_in_correct = cell(nFPs, nPorts);
ts_cent_in           = cell(nFPs, nPorts);
trialspxmat_cent_in  = cell(nFPs, nPorts);
tspkmat_cent_in      = cell(nFPs, nPorts);
rt_cent_in_sorted    = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        ind_ij = nFPs*(j-1)+i;
        t_ij = PSTHOut.CentIn.Time{ind_ij};
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

t_cent_out_correct    = cell(nFPs, nPorts);
psth_cent_out_correct = cell(nFPs, nPorts);
ts_cent_out           = cell(nFPs, nPorts);
trialspxmat_cent_out  = cell(nFPs, nPorts);
tspkmat_cent_out      = cell(nFPs, nPorts);
rt_cent_out_sorted    = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        ind_ij = nFPs*(j-1)+i;
        t_ij = PSTHOut.CentOut.Time{ind_ij};
        [psth_cent_out_correct{i,j}, ts_cent_out{i,j}, trialspxmat_cent_out{i,j}, tspkmat_cent_out{i,j}, t_cent_out_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_cent_out);
        psth_cent_out_correct{i,j} = smoothdata (psth_cent_out_correct{i,j}, 'gaussian', 5);

        rt_cent_out_sorted{i,j} = PSTHOut.CentIn.RT_Correct{i,j};
        rt_cent_out_sorted{i,j} = rt_cent_out_sorted{i,j}(ind);

        PSTH.CentOut{i,j} = {psth_cent_out_correct{i,j}, ts_cent_out{i,j}, trialspxmat_cent_out{i,j}, tspkmat_cent_out{i,j}, t_cent_out_correct{i,j}, rt_cent_out_sorted{i,j}};
    end
end
PSTH.CentOutLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT'};

% premature cent_in PSTH
t_premature_cent_in           = cell(1, nPorts);
psth_premature_cent_in        = cell(1, nPorts);
ts_premature_cent_in          = cell(1, nPorts);
trialspxmat_premature_cent_in = cell(1, nPorts);
tspkmat_premature_cent_in     = cell(1, nPorts);
hd_premature_cent_in          = cell(1, nPorts);
FP_premature_cent_in          = cell(1, nPorts);
ind_premature = find(strcmp(PSTHOut.CentIn.Labels, 'Premature'));
for j = 1:nPorts
    ind_j = ind_premature(j);
    t_j = PSTHOut.CentIn.Time{ind_j};
    [psth_premature_cent_in{j}, ts_premature_cent_in{j}, trialspxmat_premature_cent_in{j}, tspkmat_premature_cent_in{j}, t_premature_cent_in{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_in);
    psth_premature_cent_in{j} = smoothdata(psth_premature_cent_in{j}, 'gaussian', 5);

    hd_premature_cent_in{j} = PSTHOut.CentIn.HoldDur.Premature{j};
    hd_premature_cent_in{j} = hd_premature_cent_in{j}(ind);

    FP_premature_cent_in{j} = PSTHOut.CentIn.FP{2}{j};
    FP_premature_cent_in{j} = FP_premature_cent_in{j}(ind);

    PSTH.PrematureCentIn{j} = {psth_premature_cent_in{j}, ts_premature_cent_in{j}, trialspxmat_premature_cent_in{j}, tspkmat_premature_cent_in{j}, t_premature_cent_in{j}, hd_premature_cent_in{j}, FP_premature_cent_in{j}};
end
PSTH.PrematureLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP'};

% premature cent_out PSTH
t_premature_cent_out           = cell(1, nPorts);
psth_premature_cent_out        = cell(1, nPorts);
ts_premature_cent_out          = cell(1, nPorts);
trialspxmat_premature_cent_out = cell(1, nPorts);
tspkmat_premature_cent_out     = cell(1, nPorts);
hd_premature_cent_out          = cell(1, nPorts);
FP_premature_cent_out          = cell(1, nPorts);
ind_premature = find(strcmp(PSTHOut.CentOut.Labels, 'Premature'));
for j = 1:nPorts
    ind_j = ind_premature(j);
    t_j = PSTHOut.CentOut.Time{ind_j};
    [psth_premature_cent_out{j}, ts_premature_cent_out{j}, trialspxmat_premature_cent_out{j}, tspkmat_premature_cent_out{j}, t_premature_cent_out{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_out);
    psth_premature_cent_out{j} = smoothdata(psth_premature_cent_out{j}, 'gaussian', 5);

    hd_premature_cent_out{j} = PSTHOut.CentOut.HoldDur.Premature{j};
    hd_premature_cent_out{j} = hd_premature_cent_out{j}(ind);

    FP_premature_cent_out{j} = PSTHOut.CentOut.FP{2}{j};
    FP_premature_cent_out{j} = FP_premature_cent_out{j}(ind);

    PSTH.PrematureCentOut{j} = {psth_premature_cent_out{j}, ts_premature_cent_out{j}, trialspxmat_premature_cent_out{j}, tspkmat_premature_cent_out{j}, t_premature_cent_out{j}, hd_premature_cent_out{j}, FP_premature_cent_out{j}};
end

% late cent_in PSTH
t_late_cent_in           = cell(1, nPorts);
psth_late_cent_in        = cell(1, nPorts);
ts_late_cent_in          = cell(1, nPorts);
trialspxmat_late_cent_in = cell(1, nPorts);
tspkmat_late_cent_in     = cell(1, nPorts);
hd_late_cent_in          = cell(1, nPorts);
FP_late_cent_in          = cell(1, nPorts);
ind_late = find(strcmp(PSTHOut.CentIn.Labels, 'Late'));
for j = 1:nPorts
    ind_j = ind_late(j);
    t_j = PSTHOut.CentIn.Time{ind_j};
    [psth_late_cent_in{j}, ts_late_cent_in{j}, trialspxmat_late_cent_in{j}, tspkmat_late_cent_in{j}, t_late_cent_in{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_in);
    psth_late_cent_in{j} = smoothdata(psth_late_cent_in{j}, 'gaussian', 5);

    hd_late_cent_in{j} = PSTHOut.CentIn.HoldDur.Late{j};
    hd_late_cent_in{j} = hd_late_cent_in{j}(ind);

    FP_late_cent_in{j} = PSTHOut.CentIn.FP{3}{j};
    FP_late_cent_in{j} = FP_late_cent_in{j}(ind);

    PSTH.LateCentIn{j} = {psth_late_cent_in{j}, ts_late_cent_in{j}, trialspxmat_late_cent_in{j}, tspkmat_late_cent_in{j}, t_late_cent_in{j}, hd_late_cent_in{j}, FP_late_cent_in{j}};
end
PSTH.LateLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP'};

% late cent_out PSTH
t_late_cent_out           = cell(1, nPorts);
psth_late_cent_out        = cell(1, nPorts);
ts_late_cent_out          = cell(1, nPorts);
trialspxmat_late_cent_out = cell(1, nPorts);
tspkmat_late_cent_out     = cell(1, nPorts);
hd_late_cent_out          = cell(1, nPorts);
FP_late_cent_out          = cell(1, nPorts);
ind_late = find(strcmp(PSTHOut.CentOut.Labels, 'Late'));
for j = 1:nPorts
    ind_j = ind_late(j);
    t_j = PSTHOut.CentOut.Time{ind_j};
    [psth_late_cent_out{j}, ts_late_cent_out{j}, trialspxmat_late_cent_out{j}, tspkmat_late_cent_out{j}, t_late_cent_out{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_out);
    psth_late_cent_out{j} = smoothdata(psth_late_cent_out{j}, 'gaussian', 5);

    hd_late_cent_out{j} = PSTHOut.CentOut.HoldDur.Late{j};
    hd_late_cent_out{j} = hd_late_cent_out{j}(ind);

    FP_late_cent_out{j} = PSTHOut.CentOut.FP{3}{j};
    FP_late_cent_out{j} = FP_late_cent_out{j}(ind);

    PSTH.LateCentOut{j} = {psth_late_cent_out{j}, ts_late_cent_out{j}, trialspxmat_late_cent_out{j}, tspkmat_late_cent_out{j}, t_late_cent_out{j}, hd_late_cent_out{j}, FP_late_cent_out{j}};
end

% probe cent_in PSTH
t_probe_cent_in           = cell(1, nPorts);
psth_probe_cent_in        = cell(1, nPorts);
ts_probe_cent_in          = cell(1, nPorts);
trialspxmat_probe_cent_in = cell(1, nPorts);
tspkmat_probe_cent_in     = cell(1, nPorts);
hd_probe_cent_in          = cell(1, nPorts);
ind_probe = find(strcmp(PSTHOut.CentIn.Labels, 'Probe'));
for j = 1:nPorts
    ind_j = ind_probe(j);
    t_j = PSTHOut.CentIn.Time{ind_j};
    [psth_probe_cent_in{j}, ts_probe_cent_in{j}, trialspxmat_probe_cent_in{j}, tspkmat_probe_cent_in{j}, t_probe_cent_in{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_in);
    psth_probe_cent_in{j} = smoothdata(psth_probe_cent_in{j}, 'gaussian', 5);

    hd_probe_cent_in{j} = PSTHOut.CentIn.HoldDur.Probe{j};
    hd_probe_cent_in{j} = hd_probe_cent_in{j}(ind);

    PSTH.ProbeCentIn{j} = {psth_probe_cent_in{j}, ts_probe_cent_in{j}, trialspxmat_probe_cent_in{j}, tspkmat_probe_cent_in{j}, t_probe_cent_in{j}, hd_probe_cent_in{j}};
end

% late cent_out PSTH
t_probe_cent_out           = cell(1, nPorts);
psth_probe_cent_out        = cell(1, nPorts);
ts_probe_cent_out          = cell(1, nPorts);
trialspxmat_probe_cent_out = cell(1, nPorts);
tspkmat_probe_cent_out     = cell(1, nPorts);
hd_probe_cent_out          = cell(1, nPorts);
ind_probe = find(strcmp(PSTHOut.CentIn.Labels, 'Probe'));
for j = 1:nPorts
    ind_j = ind_probe(j);
    t_j = PSTHOut.CentOut.Time{ind_j};
    [psth_probe_cent_out{j}, ts_probe_cent_out{j}, trialspxmat_probe_cent_out{j}, tspkmat_probe_cent_out{j}, t_probe_cent_out{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_cent_out);
    psth_probe_cent_out{j} = smoothdata(psth_probe_cent_out{j}, 'gaussian', 5);

    hd_probe_cent_out{j} = PSTHOut.CentOut.HoldDur.Probe{j};
    hd_probe_cent_out{j} = hd_probe_cent_out{j}(ind);

    PSTH.ProbeCentOut{j} = {psth_probe_cent_out{j}, ts_probe_cent_out{j}, trialspxmat_probe_cent_out{j}, tspkmat_probe_cent_out{j}, t_probe_cent_out{j}, hd_probe_cent_out{j}};
end
PSTH.ProbeLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration'};

% use t_reward_poke and move_time to construct reward_poke PSTH
% reward PSTH
params_choice.pre  = ChoiceTimeDomain(1);
params_choice.post = ChoiceTimeDomain(2);
params_choice.binwidth = 20;

t_reward_choice           = cell(nFPs, nPorts);
psth_reward_choice        = cell(nFPs, nPorts);
ts_reward_choice          = cell(nFPs, nPorts);
trialspxmat_reward_choice = cell(nFPs, nPorts);
tspkmat_reward_choice     = cell(nFPs, nPorts);
mt_reward_choice          = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        ind_ij = nFPs*(j-1)+i;
        t_ij = PSTHOut.ChoiceIn.RewardPoke.Time{ind_ij};
        [psth_reward_choice{i,j}, ts_reward_choice{i,j}, trialspxmat_reward_choice{i,j}, tspkmat_reward_choice{i,j}, t_reward_choice{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_choice);
        psth_reward_choice{i,j} = smoothdata(psth_reward_choice{i,j}, 'gaussian', 5);

        mt_reward_choice{i,j} = PSTHOut.ChoiceIn.RewardPoke.Move_Time{i,j}(ind);
        PSTH.RewardChoice{i,j} = {psth_reward_choice{i,j}, ts_reward_choice{i,j}, trialspxmat_reward_choice{i,j}, tspkmat_reward_choice{i,j}, t_reward_choice{i,j}, mt_reward_choice{i,j}};
    end
end
PSTH.ChoiceLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'MoveTime'};

% nonreward choice PSTH
t_nonreward_choice           = cell(1, nPorts);
psth_nonreward_choice        = cell(1, nPorts);
ts_nonreward_choice          = cell(1, nPorts);
trialspxmat_nonreward_choice = cell(1, nPorts);
tspkmat_nonreward_choice     = cell(1, nPorts);
mt_nonreward_choice          = cell(1, nPorts);
for j = 1:nPorts
    t_j = PSTHOut.ChoiceIn.NonrewardPoke.Time{j};
    [psth_nonreward_choice{j}, ts_nonreward_choice{j}, trialspxmat_nonreward_choice{j}, tspkmat_nonreward_choice{j}, t_nonreward_choice{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_choice);
    psth_nonreward_choice{j} = smoothdata(psth_nonreward_choice{j}, 'gaussian', 5);

    mt_nonreward_choice{j} = PSTHOut.ChoiceIn.NonrewardPoke.Move_Time{j}(ind);
    PSTH.NonrewardChoice{j} = {psth_nonreward_choice{j}, ts_nonreward_choice{j}, trialspxmat_nonreward_choice{j}, tspkmat_nonreward_choice{j}, t_nonreward_choice{j}, mt_nonreward_choice{j}};
end

% trigger PSTH
params_trigger.pre  = TriggerTimeDomain(1);
params_trigger.post = TriggerTimeDomain(2);
params_trigger.binwidth = 20;

% late response
% t_trigger_late    = PSTHOut.Triggers.Time{end};
% RT_trigger_late   = PSTHOut.Triggers.RT{end};
% FP_trigger_late   = PSTHOut.Triggers.FP{end};
% Port_trigger_late = PSTHOut.Triggers.Port{end};
% 
% [psth_late_trigger, ts_late_trigger, trialspxmat_late_trigger, tspkmat_late_trigger, t_trigger_late, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_trigger_late, params_trigger);
% psth_late_trigger = smoothdata (psth_late_trigger, 'gaussian', 5);
% 
% RT_trigger_late   = RT_trigger_late(ind);
% FP_trigger_late   = FP_trigger_late(ind);
% Port_trigger_late = Port_trigger_late(ind);
% 
% PSTH.TriggerLate   = {psth_late_trigger, ts_late_trigger, trialspxmat_late_trigger, tspkmat_late_trigger, t_trigger_late, RT_trigger_late, FP_trigger_late, Port_trigger_late};

% correct response
t_trigger_correct           = cell(nFPs, nPorts);
psth_trigger_correct        = cell(nFPs, nPorts);
ts_trigger_correct          = cell(nFPs, nPorts);
trialspxmat_trigger_correct = cell(nFPs, nPorts);
tspkmat_trigger_correct     = cell(nFPs, nPorts);
rt_trigger_correct          = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        ind_ij = nFPs*(j-1)+i;
        t_ij = PSTHOut.Triggers.Time{ind_ij};
        [psth_trigger_correct{i,j}, ts_trigger_correct{i,j}, trialspxmat_trigger_correct{i,j}, tspkmat_trigger_correct{i,j}, t_trigger_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params_trigger);
        psth_trigger_correct{i,j} = smoothdata (psth_trigger_correct{i,j}, 'gaussian', 5);

        rt_trigger_correct{i,j} = PSTHOut.Triggers.RT{ind_ij};
        rt_trigger_correct{i,j} = rt_trigger_correct{i,j}(ind);
        PSTH.Triggers{i,j} = {psth_trigger_correct{i,j}, ts_trigger_correct{i,j}, trialspxmat_trigger_correct{i,j}, tspkmat_trigger_correct{i,j}, t_trigger_correct{i,j}, rt_trigger_correct{i,j}, TargetFPs(i), Ports(j)};
    end
end
PSTH.TriggerLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT', 'FP', 'Port'};

% late response
t_trigger_late           = cell(1, nPorts);
psth_trigger_late        = cell(1, nPorts);
ts_trigger_late          = cell(1, nPorts);
trialspxmat_trigger_late = cell(1, nPorts);
tspkmat_trigger_late     = cell(1, nPorts);
rt_trigger_late          = cell(1, nPorts);
FP_trigger_late          = cell(1, nPorts);
for j = 1:nPorts
    t_j = PSTHOut.Triggers.Time{end+j-nPorts};
    [psth_trigger_late{j}, ts_trigger_late{j}, trialspxmat_trigger_late{j}, tspkmat_trigger_late{j}, t_trigger_late{j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_j, params_trigger);
    psth_trigger_late{j} = smoothdata (psth_trigger_late{j}, 'gaussian', 5);

    rt_trigger_late{j} = PSTHOut.Triggers.RT{end+j-nPorts};
    rt_trigger_late{j} = rt_trigger_late{j}(ind);

    FP_trigger_late{j} = PSTHOut.Triggers.FP{end}{j};
    FP_trigger_late{j} = FP_trigger_late{j}(ind);

    PSTH.TriggerLate{j} = {psth_trigger_late{j}, ts_trigger_late{j}, trialspxmat_trigger_late{j}, tspkmat_trigger_late{j}, t_trigger_late{j}, rt_trigger_late{j}, FP_trigger_late{j}, Ports(j)};
end

%% Plot raster and spks

h_psth = 1;
fig = figure(40); clf(fig);
set(fig, 'unit', 'centimeters', 'position', printsize, 'paperpositionmode', 'auto' ,'color', 'w', 'Visible', 'on')

%% Align to CentIn

w_psth = sum(CentInTimeDomain) / 1000;
% Correct trials

% PSTH of correct trials
yshift_row1 = 1;
FRMax = 3;

ha_cent_in_psth = cell(1,nFPs);
for i = 1:nFPs
    yshift = yshift_row1 + (h_psth+0.1)*(nFPs-i);
    ha_cent_in_psth{i} = axes(fig, 'unit', 'centimeters', 'position', [1.25 yshift w_psth h_psth], 'nextplot', 'add', 'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'FontSize', 7, 'TickDir', 'Out');
    xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
    xline(TargetFPs(i)*1000, 'Color', c_trigger, 'LineWidth', 1, 'Alpha', 1);
    for j = 1:nPorts
        plot(ts_cent_in{i,j}, psth_cent_in_correct{i,j}, 'color', c_port{j}, 'linewidth', 1);
        FRMax = max([FRMax max(psth_cent_in_correct{i,j})]);
        %     disp(FRMax)
    end
    if i==nFPs
        xlabel('Time from cent-in (ms)')
        ylabel('Spks per s')
    else
        xticklabels([]);
        yticklabels([]);
    end
end
axis 'auto y'

% make raster plot  750 ms FP
if num2str(length(t_cent_in))>200
    rasterheight = 0.01;
elseif num2str(length(t_cent_in))>100
    rasterheight = 0.02;
else
    rasterheight = 0.04;
end

yshift_row2 = yshift_row1 + (h_psth+.1)*nFPs;
% Plot spike raster of correct trials (all FPs)
ntrials_cent_in = 0;
nFP_ij = zeros(nFPs, nPorts);
t_choice = PSTHOut.ChoiceIn.Time;
for i = 1:nFPs
    for j = 1:nPorts
        nFP_ij(i,j) = size(trialspxmat_cent_in{i,j}, 2);
        ntrials_cent_in = ntrials_cent_in + nFP_ij(i,j);
    end
end
axes('unit', 'centimeters', 'position', [1.25 yshift_row2 w_psth ntrials_cent_in*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrials_cent_in 1], 'box', 'on');
% Paint the foreperiod
k=0;
for m = 1:nFPs
    for n = 1:nPorts
        ap_mat = trialspxmat_cent_in{m,n};
        t_mat = tspkmat_cent_in{m,n};
        rt = rt_cent_in_sorted{m,n};
        xx_all = [];
        yy_all = [];
        xxrt_all = [];
        yyrt_all = [];
        x_portin = [];
        y_portin = [];
        for i = 1:nFP_ij(m,n)
            irt = rt(i); % time from foreperiod to release
            xx = t_mat(ap_mat(:,i)==1);
            yy1 = [0 0.8]-k;
            yy2 = [0 1]-k;
            xxrt = irt+TargetFPs(m)*1000;
%             fill(1000*[0 TargetFPs(m) TargetFPs(m) 0],[-k -k 1-k 1-k], 'r', 'FaceColor', c_trigger, 'FaceAlpha', 1, 'EdgeColor', 'none');
            line(1000*[TargetFPs(m) TargetFPs(m)], [-k 1-k], 'Color', c_trigger, 'linewidth', 1);

            if isempty(find(isnan(ap_mat(:, i)), 1))
                for i_xx = 1:length(xx)
                    xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
                    yy_all = [yy_all, yy1, NaN];
                end
                xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
                yyrt_all = [yyrt_all, yy2, NaN];
            end
            % port access time
            it_cent_in = t_cent_in_correct{m,n}(i);
            it_choice = t_choice - it_cent_in;
            it_choice = it_choice(it_choice>=-CentInTimeDomain(1) & it_choice<=CentInTimeDomain(2));
            if ~isempty(it_choice)
                it_choice = reshape(it_choice,1,[]);
                x_portin = [x_portin, it_choice];
                y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
            end
            k = k+1;
        end
        line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2);
        line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1);
        scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
    end
end

xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
title('Correct', 'fontsize', 7, 'fontweight','bold');
axis off

% Probe trials

% PSTH of probe trials
yshift_row3 = yshift_row2+ntrials_cent_in*rasterheight+0.5;
ha_cent_in_psth_probe = axes('unit', 'centimeters', 'position', [1.25 yshift_row3 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-CentInTimeDomain(1) CentInTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
% plot probe trials as well
for j = 1:nPorts
    if size(trialspxmat_probe_cent_in{j}, 2)>3
        plot(ts_probe_cent_in{j}, psth_probe_cent_in{j}, 'color', c_port{j}, 'linewidth',1);
%         FRMax = max([FRMax max(psth_premature_cent_in)]);
        %      disp(FRMax)
    end
end

axis 'auto y'
yshift_row4 = yshift_row3 + (h_psth+.1);
% Probe cent-in raster plot
ntrial_probe = sum(cellfun(@(x) size(x, 2), trialspxmat_probe_cent_in)); % number of trials
axes('unit', 'centimeters', 'position', [1.25 yshift_row4 w_psth ntrial_probe*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrial_probe 1], 'box', 'on');
k =0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_probe_cent_in{n};
    t_mat             =     tspkmat_probe_cent_in{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = hd_probe_cent_in{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = ipredur;
        % plot trigger stimulus FPs_premature_cent_in
%         itrigger = FP_probe_cent_in{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_probe_cent_in{n}(i);
        it_choice = it_choice(it_choice>=-CentInTimeDomain(1) & it_choice<=CentInTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
title('Probe', 'fontsize', 7, 'fontweight','bold')
axis off

% Late trials

% PSTH of late trials
yshift_row5 = yshift_row4 + 0.5 + ntrial_probe*rasterheight;

ha_cent_in_psth_late =  axes('unit', 'centimeters', 'position', [1.25 yshift_row5 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-CentInTimeDomain(1) CentInTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if  size(trialspxmat_late_cent_in{j}, 2)>3
        plot(ts_late_cent_in{j}, psth_late_cent_in{j}, 'color', c_port{j}, 'linewidth', 1)
        %     FRMax = max([FRMax max(psth_late_cent_in)]);
        %     disp(FRMax)
    end
end
axis 'auto y'

hline_cent_in_error = line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1);

% Late response raster plot
yshift_row6 = yshift_row5 + (h_psth+.1);
ntrial_late = sum(cellfun(@(x) size(x, 2), trialspxmat_late_cent_in)); % number of trials
axes('unit', 'centimeters', 'position', [1.25 yshift_row6 w_psth ntrial_late*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrial_late 1], 'box', 'on');
k = 0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_late_cent_in{n};
    t_mat             =     tspkmat_late_cent_in{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = hd_late_cent_in{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = ipredur;
        % plot trigger stimulus FPs_premature_cent_in
%         itrigger = FP_late_cent_in{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_late_cent_in{n}(i);
        it_choice = it_choice(it_choice>=-CentInTimeDomain(1) & it_choice<=CentInTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
title('Late', 'fontsize', 7, 'fontweight','bold')
axis off

% Premature trials

% PSTH
yshift_row7             =      yshift_row6 + 0.5 + ntrial_late*rasterheight;
ha_cent_in_psth_premature =  axes('unit', 'centimeters', 'position', [1.25 yshift_row7 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-CentInTimeDomain(1) CentInTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if size(trialspxmat_premature_cent_in{j}, 2)>3
        plot(ts_premature_cent_in{j}, psth_premature_cent_in{j}, 'color', c_port{j}, 'linewidth',1);
        %     FRMax = max([FRMax max(psth_premature_cent_in)]);
        %      disp(FRMax)
    end
end

% Premature press raster plot
yshift_row8 = yshift_row7 + (h_psth+.1);
ntrial_premature = sum(cellfun(@(x) size(x, 2), trialspxmat_premature_cent_in)); % number of trials
axes('unit', 'centimeters', 'position', [1.25 yshift_row8 w_psth ntrial_premature*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrial_premature 1], 'box', 'on');
yshift_row9    =      yshift_row8 + 0.5 + ntrial_premature*rasterheight;
k = 0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_premature_cent_in{n};
    t_mat             =     tspkmat_premature_cent_in{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = hd_premature_cent_in{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = ipredur;
        % plot trigger stimulus FPs_premature_cent_in
%         itrigger = FP_late_cent_in{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_premature_cent_in{n}(i);
        it_choice = it_choice(it_choice>=-CentInTimeDomain(1) & it_choice<=CentInTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_cent_in, 'LineWidth', 1, 'Alpha', 1);
title('Premature', 'fontsize', 7, 'fontweight','bold')
axis off

% this is the position of last panel
% Add information
uicontrol('Style','text','Units','centimeters','Position',[0.75 yshift_row9  6 1],...
    'string', 'A. CentIn-related', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 8,'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','Left');

yshift_row10=yshift_row9+1.25;
ch = r.Units.SpikeNotes(ku, 1);
unit_no = r.Units.SpikeNotes(ku, 2);

if size(r.Units.SpikeNotes, 2) == 4
    cluster_id = r.Units.SpikeNotes(ku, 4);
    uicontrol('style', 'text', 'units', 'centimeters', 'position', [1 yshift_row10 6 1.2],...
        'string', (['Unit #' num2str(ku) ' (Ch ' num2str(ch) ' | UnitOnCh ' num2str(unit_no) ' | ' 'Kilosort cluster ' num2str(cluster_id) ')']),...
        'BackgroundColor','w', 'fontsize', 8, 'fontweight','bold',  'FontName','Dejavu Sans')
else
    cluster_id = [];
    uicontrol('style', 'text', 'units', 'centimeters', 'position', [1 yshift_row10 6 1.2],...
        'string', (['Unit #' num2str(ku) ' (' num2str(ch) ' | ' num2str(unit_no) ')']),...
        'BackgroundColor','w', 'fontsize', 8, 'fontweight','bold',  'FontName','Dejavu Sans')
end
uicontrol('style', 'text', 'units', 'centimeters', 'position', [1 yshift_row10+1.2 4 0.5],...
    'string', ([r.Meta(1).Subject ' ' r.Meta(1).DateTime(1:11)]), 'BackgroundColor','w',...
    'fontsize', 8, 'fontweight', 'bold',  'FontName','Dejavu Sans')

fig_height = yshift_row10+2;

%% Align to CentOut

x_pos = 6;
w_psth = sum(CentOutTimeDomain) / 1000;
% Correct trials

% PSTH of correct trials
yshift_row1 = 1;
ha_cent_out_psth = cell(1,nFPs);
for i = 1:nFPs
    yshift = yshift_row1 + (h_psth+0.1)*(nFPs-i);
    ha_cent_out_psth{i} = axes(fig, 'unit', 'centimeters', 'position', [x_pos yshift w_psth h_psth], 'nextplot', 'add', 'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'FontSize', 7, 'TickDir', 'Out');
    xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
    for j = 1:nPorts
        plot(ts_cent_out{i,j}, psth_cent_out_correct{i,j}, 'color', c_port{j}, 'linewidth', 1);
        FRMax = max([FRMax max(psth_cent_out_correct{i,j})]);
        %     disp(FRMax)
    end
    if i==nFPs
        xlabel('Time from cent-out (ms)')
        ylabel('Spks per s')
    else
        xticklabels([]);
        yticklabels([]);
    end
end
axis 'auto y'

% make raster plot  750 ms FP
if num2str(length(t_cent_out_correct))>200
    rasterheight = 0.01;
elseif num2str(length(t_cent_out_correct))>100
    rasterheight = 0.02;
else
    rasterheight = 0.04;
end

yshift_row2 = yshift_row1 + (h_psth+.1)*nFPs;
% Plot spike raster of correct trials (all FPs)
ntrials_cent_out = 0;
nFP_ij = zeros(nFPs, nPorts);
t_choice = PSTHOut.ChoiceIn.Time;
for i = 1:nFPs
    for j = 1:nPorts
        nFP_ij(i,j) = size(trialspxmat_cent_out{i,j}, 2);
        ntrials_cent_out = ntrials_cent_out + nFP_ij(i,j);
    end
end
axes('unit', 'centimeters', 'position', [x_pos yshift_row2 w_psth ntrials_cent_out*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrials_cent_out 1], 'box', 'on');
% Paint the foreperiod
k=0;
for m = 1:nFPs
    for n = 1:nPorts
        ap_mat = trialspxmat_cent_out{m,n};
        t_mat = tspkmat_cent_out{m,n};
        rt = rt_cent_out_sorted{m,n};
        xx_all = [];
        yy_all = [];
        xxrt_all = [];
        yyrt_all = [];
        x_portin = [];
        y_portin = [];
        for i = 1:nFP_ij(m,n)
            irt = rt(i); % time from foreperiod to release
            xx = t_mat(ap_mat(:,i)==1);
            yy1 = [0 0.8]-k;
            yy2 = [0 1]-k;
            xxrt = -irt;
%             fill(1000*[0 TargetFPs(m) TargetFPs(m) 0],[-k -k 1-k 1-k], 'r', 'FaceColor', c_trigger, 'FaceAlpha', .4, 'EdgeColor', 'none');

            if isempty(find(isnan(ap_mat(:, i)), 1))
                for i_xx = 1:length(xx)
                    xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
                    yy_all = [yy_all, yy1, NaN];
                end
                xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
                yyrt_all = [yyrt_all, yy2, NaN];
            end
            % port access time
            it_cent_out = t_cent_out_correct{m,n}(i);
            it_choice = t_choice - it_cent_out;
            it_choice = it_choice(it_choice>=-CentOutTimeDomain(1) & it_choice<=CentOutTimeDomain(2));
            if ~isempty(it_choice)
                it_choice = reshape(it_choice,1,[]);
                x_portin = [x_portin, it_choice];
                y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
            end
            k = k+1;
        end
        line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2);
        line(xxrt_all, yyrt_all, 'color', c_trigger, 'linewidth', 1);
        scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
    end
end

xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
title('Correct', 'fontsize', 7, 'fontweight','bold');
axis off

% Probe trials

% PSTH of probe trials
yshift_row3 = yshift_row2+ntrials_cent_out*rasterheight+0.5;
ha_cent_out_psth_probe = axes('unit', 'centimeters', 'position', [x_pos yshift_row3 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
% plot probe trials as well
xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if size(trialspxmat_probe_cent_out{j}, 2)>3
        plot(ts_probe_cent_out{j}, psth_probe_cent_out{j}, 'color', c_port{j}, 'linewidth',1);
        %     FRMax = max([FRMax max(psth_premature_cent_out)]);
        %      disp(FRMax)
    end
end

axis 'auto y'
yshift_row4 = yshift_row3 + (h_psth+.1);
% Probe cent-in raster plot
ntrial_probe = sum(cellfun(@(x) size(x, 2), trialspxmat_probe_cent_out)); % number of trials
axes('unit', 'centimeters', 'position', [x_pos yshift_row4 w_psth ntrial_probe*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_probe 1], 'box', 'on');
k =0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_probe_cent_out{n};
    t_mat             =     tspkmat_probe_cent_out{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = hd_probe_cent_out{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = ipredur;
        % plot trigger stimulus FPs_premature_cent_out
%         itrigger = FP_probe_cent_out{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_probe_cent_out{n}(i);
        it_choice = it_choice(it_choice>=-CentOutTimeDomain(1) & it_choice<=CentOutTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
title('Probe', 'fontsize', 7, 'fontweight','bold')
axis off

% Late trials

% PSTH of late trials
yshift_row5 = yshift_row4 + 0.5 + ntrial_probe*rasterheight;

ha_cent_out_psth_late =  axes('unit', 'centimeters', 'position', [x_pos yshift_row5 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if  size(trialspxmat_late_cent_out{j}, 2)>3
        plot(ts_late_cent_out{j}, psth_late_cent_out{j}, 'color', c_port{j}, 'linewidth', 1)
        %     FRMax = max([FRMax max(psth_late_cent_out)]);
        %     disp(FRMax)
    end
end
axis 'auto y'

% Late response raster plot
yshift_row6 = yshift_row5 + (h_psth+.1);
ntrial_late = sum(cellfun(@(x) size(x, 2), trialspxmat_late_cent_out)); % number of trials
axes('unit', 'centimeters', 'position', [x_pos yshift_row6 w_psth ntrial_late*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_late 1], 'box', 'on');
k = 0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_late_cent_out{n};
    t_mat             =     tspkmat_late_cent_out{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = hd_late_cent_out{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = ipredur;
        % plot trigger stimulus FPs_premature_cent_out
%         itrigger = FP_late_cent_out{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_late_cent_out{n}(i);
        it_choice = it_choice(it_choice>=-CentOutTimeDomain(1) & it_choice<=CentOutTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
title('Late', 'fontsize', 7, 'fontweight','bold')
axis off

% Premature trials

% PSTH
yshift_row7             =      yshift_row6 + 0.5 + ntrial_late*rasterheight;
ha_cent_out_psth_premature =  axes('unit', 'centimeters', 'position', [x_pos yshift_row7 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if size(trialspxmat_premature_cent_out{j}, 2)>3
        plot(ts_premature_cent_out{j}, psth_premature_cent_out{j}, 'color', c_port{j}, 'linewidth',1);
        %     FRMax = max([FRMax max(psth_premature_cent_out)]);
        %      disp(FRMax)
    end
end

% Premature press raster plot
yshift_row8 = yshift_row7 + (h_psth+.1);
ntrial_premature = sum(cellfun(@(x) size(x, 2), trialspxmat_premature_cent_out)); % number of trials
axes('unit', 'centimeters', 'position', [x_pos yshift_row8 w_psth ntrial_premature*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_premature 1], 'box', 'on');
yshift_row9    =      yshift_row8 + 0.5 + ntrial_premature*rasterheight;
k = 0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_premature_cent_out{n};
    t_mat             =     tspkmat_premature_cent_out{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = hd_premature_cent_out{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = -ipredur;
        % plot trigger stimulus FPs_premature_cent_out
%         itrigger = FP_late_cent_out{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_premature_cent_out{n}(i);
        it_choice = it_choice(it_choice>=-CentOutTimeDomain(1) & it_choice<=CentOutTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_in, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_cent_out, 'LineWidth', 1, 'Alpha', 1);
title('Premature', 'fontsize', 7, 'fontweight','bold')
axis off

% this is the position of last panel
% Add information
uicontrol('Style','text','Units','centimeters','Position',[x_pos-.25 yshift_row9 4 1],...
    'string', 'B. CentOut-related', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 8,'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','Left');

%% Align to Trigger

x_pos = 9.75;
w_psth = 2 * sum(TriggerTimeDomain) / 1000;
% Correct trials

% PSTH of correct trials
yshift_row1 = 1;
ha_trigger_psth = cell(1,nFPs);
for i = 1:nFPs
    yshift = yshift_row1 + (h_psth+0.1)*(nFPs-i);
    ha_trigger_psth{i} = axes(fig, 'unit', 'centimeters', 'position', [x_pos yshift w_psth h_psth], 'nextplot', 'add', 'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'FontSize', 7, 'TickDir', 'Out');
    xline(0, 'Color', c_trigger, 'LineWidth', 1, 'Alpha', 1);
    for j = 1:nPorts
        plot(ts_trigger_correct{i,j}, psth_trigger_correct{i,j}, 'color', c_port{j}, 'linewidth', 1);
        FRMax = max([FRMax max(psth_trigger_correct{i,j})]);
        %     disp(FRMax)
    end
    if i==nFPs
        xlabel('Time from trigger (ms)')
        ylabel('Spks per s')
    else
        xticklabels([]);
        yticklabels([]);
    end
end
axis 'auto y'

% make raster plot  750 ms FP
if num2str(length(t_trigger_correct))>200
    rasterheight = 0.01;
elseif num2str(length(t_trigger_correct))>100
    rasterheight = 0.02;
else
    rasterheight = 0.04;
end

yshift_row2 = yshift_row1 + (h_psth+.1)*nFPs;
% Plot spike raster of correct trials (all FPs)
ntrials_trigger = 0;
nFP_ij = zeros(nFPs, nPorts);
t_choice = PSTHOut.ChoiceIn.Time;
for i = 1:nFPs
    for j = 1:nPorts
        nFP_ij(i,j) = size(trialspxmat_trigger_correct{i,j}, 2);
        ntrials_trigger = ntrials_trigger + nFP_ij(i,j);
    end
end
axes('unit', 'centimeters', 'position', [x_pos yshift_row2 w_psth ntrials_trigger*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'ylim', [-ntrials_trigger 1], 'box', 'on');
% Paint the foreperiod
k=0;
for m = 1:nFPs
    for n = 1:nPorts
        ap_mat = trialspxmat_trigger_correct{m,n};
        t_mat = tspkmat_trigger_correct{m,n};
        rt = rt_trigger_correct{m,n};
        xx_all = [];
        yy_all = [];
        xxrt_all = [];
        yyrt_all = [];
        x_portin = [];
        y_portin = [];
        for i = 1:nFP_ij(m,n)
            irt = rt(i); % time from foreperiod to release
            xx = t_mat(ap_mat(:,i)==1);
            yy1 = [0 0.8]-k;
            yy2 = [0 1]-k;
            xxrt = irt;
%             fill(1000*[0 TargetFPs(m) TargetFPs(m) 0],[-k -k 1-k 1-k], 'r', 'FaceColor', c_trigger, 'FaceAlpha', .4, 'EdgeColor', 'none');

            if isempty(find(isnan(ap_mat(:, i)), 1))
                for i_xx = 1:length(xx)
                    xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
                    yy_all = [yy_all, yy1, NaN];
                end
                xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
                yyrt_all = [yyrt_all, yy2, NaN];
            end
            % port access time
            it_trigger = t_trigger_correct{m,n}(i);
            it_choice = t_choice - it_trigger;
            it_choice = it_choice(it_choice>=-TriggerTimeDomain(1) & it_choice<=TriggerTimeDomain(2));
            if ~isempty(it_choice)
                it_choice = reshape(it_choice,1,[]);
                x_portin = [x_portin, it_choice];
                y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
            end
            k = k+1;
        end
        line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2);
        line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1);
        scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
    end
end

xline(0, 'Color', c_trigger, 'LineWidth', 1, 'Alpha', 1);
title('Correct', 'fontsize', 7, 'fontweight','bold');
axis off
% Late trials

% PSTH of late trials
yshift_row3 = yshift_row2 + 0.5 + ntrials_trigger*rasterheight;

ha_trigger_psth_late =  axes('unit', 'centimeters', 'position', [x_pos yshift_row3 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_trigger, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if  size(trialspxmat_trigger_late{j}, 2)>3
        plot(ts_trigger_late{j}, psth_trigger_late{j}, 'color', c_port{j}, 'linewidth', 1)
        %     FRMax = max([FRMax max(psth_late_cent_out)]);
        %     disp(FRMax)
    end
end
axis 'auto y'

% Late response raster plot
yshift_row4 = yshift_row3 + (h_psth+.1);
ntrial_trigger_late = sum(cellfun(@(x) size(x, 2), trialspxmat_trigger_late)); % number of trials
axes('unit', 'centimeters', 'position', [x_pos yshift_row4 w_psth ntrial_trigger_late*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'ylim', [-ntrial_trigger_late 1], 'box', 'on');
k = 0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_trigger_late{n};
    t_mat             =     tspkmat_trigger_late{n};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        ipredur = rt_trigger_late{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = ipredur;
        % plot trigger stimulus FPs_premature_cent_out
%         itrigger = FP_late_cent_out{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % plot port poke time
        it_choice = t_choice - t_trigger_late{n}(i);
        it_choice = it_choice(it_choice>=-TriggerTimeDomain(1) & it_choice<=TriggerTimeDomain(2));
        if ~isempty(it_choice)
            it_choice = reshape(it_choice,1,[]);
            x_portin = [x_portin, it_choice];
            y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
        end
        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end

xline(0, 'Color', c_trigger, 'LineWidth', 1, 'Alpha', 1);
title('Late', 'fontsize', 7, 'fontweight','bold')
axis off

yshift_row5    =      yshift_row4 + 0.5 + ntrial_trigger_late*rasterheight;
% Add information
uicontrol('Style','text','Units','centimeters','Position',[x_pos-.25 yshift_row5 4 1],...
    'string', 'C. Trigger-related', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 8,'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','Left');

%% Align to ChoiceIn

x_pos = 13;
w_psth = sum(ChoiceTimeDomain) / 1000;
% Correct trials

% PSTH of correct trials
yshift_row1 = 1;
ha_choice_psth = cell(1,nFPs);
for i = 1:nFPs
    yshift = yshift_row1 + (h_psth+0.1)*(nFPs-i);
    ha_choice_psth{i} = axes(fig, 'unit', 'centimeters', 'position', [x_pos yshift w_psth h_psth], 'nextplot', 'add', 'xlim', [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], 'FontSize', 7, 'TickDir', 'Out');
    xline(0, 'Color', c_reward, 'LineWidth', 1, 'Alpha', 1);
    for j = 1:nPorts
        plot(ts_reward_choice{i,j}, psth_reward_choice{i,j}, 'color', c_port{j}, 'linewidth', 1);
        FRMax = max([FRMax max(psth_reward_choice{i,j})]);
        %     disp(FRMax)
    end
    if i==nFPs
        xlabel('Time from Choice (ms)')
        ylabel('Spks per s')
    else
        xticklabels([]);
        yticklabels([]);
    end
end
axis 'auto y'

% make raster plot  750 ms FP
if num2str(length(t_reward_choice))>200
    rasterheight = 0.01;
elseif num2str(length(t_reward_choice))>100
    rasterheight = 0.02;
else
    rasterheight = 0.04;
end

yshift_row2 = yshift_row1 + (h_psth+.1)*nFPs;
% Plot spike raster of correct trials (all FPs)
ntrials_reward_choice = 0;
nFP_ij = zeros(nFPs, nPorts);
t_choice = PSTHOut.ChoiceIn.Time;
for i = 1:nFPs
    for j = 1:nPorts
        nFP_ij(i,j) = size(trialspxmat_reward_choice{i,j}, 2);
        ntrials_reward_choice = ntrials_reward_choice + nFP_ij(i,j);
    end
end
axes('unit', 'centimeters', 'position', [x_pos yshift_row2 w_psth ntrials_reward_choice*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], 'ylim', [-ntrials_reward_choice 1], 'box', 'on');
% Paint the foreperiod
k=0;
for m = 1:nFPs
    for n = 1:nPorts
        ap_mat = trialspxmat_reward_choice{m,n};
        t_mat = tspkmat_reward_choice{m,n};
        mt = mt_reward_choice{m,n};
        xx_all = [];
        yy_all = [];
        xxmt_all = [];
        yymt_all = [];
        x_portin = [];
        y_portin = [];
        for i = 1:nFP_ij(m,n)
            imt = mt(i); % time from foreperiod to release
            xx = t_mat(ap_mat(:,i)==1);
            yy1 = [0 0.8]-k;
            yy2 = [0 1]-k;
            xxmt = -imt;
%             fill(1000*[0 TargetFPs(m) TargetFPs(m) 0],[-k -k 1-k 1-k], 'r', 'FaceColor', c_trigger, 'FaceAlpha', .4, 'EdgeColor', 'none');

            if isempty(find(isnan(ap_mat(:, i)), 1))
                for i_xx = 1:length(xx)
                    xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
                    yy_all = [yy_all, yy1, NaN];
                end
                xxmt_all = [xxmt_all, xxmt, xxmt, NaN];
                yymt_all = [yymt_all, yy2, NaN];
            end
            k = k+1;
        end
        line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2);
        line(xxmt_all, yymt_all, 'color', c_cent_out, 'linewidth', 1);
    end
end

xline(0, 'Color', c_reward, 'LineWidth', 1, 'Alpha', 1);
title('Reward', 'fontsize', 7, 'fontweight','bold');
axis off
% Non-reward trials

% PSTH of non-reward trials
yshift_row3 = yshift_row2 + 0.5 + ntrials_reward_choice*rasterheight;

ha_nonreward_choice_psth =  axes('unit', 'centimeters', 'position', [x_pos yshift_row3 w_psth h_psth], 'nextplot', 'add',...
    'xlim',  [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], 'xticklabel', [], 'FontSize', 7, 'TickDir', 'Out');
xline(0, 'Color', c_reward, 'LineWidth', 1, 'Alpha', 1);
for j = 1:nPorts
    if  size(trialspxmat_nonreward_choice{j}, 2)>3
        plot(ts_nonreward_choice{j}, psth_nonreward_choice{j}, 'color', c_port{j}, 'linewidth', 1)
        %     FRMax = max([FRMax max(psth_late_cent_out)]);
        %     disp(FRMax)
    end
end
axis 'auto y'

% Late response raster plot
yshift_row4 = yshift_row3 + (h_psth+.1);
ntrial_nonreward_choice = sum(cellfun(@(x) size(x, 2), trialspxmat_nonreward_choice)); % number of trials
axes('unit', 'centimeters', 'position', [x_pos yshift_row4 w_psth ntrial_nonreward_choice*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], 'ylim', [-ntrial_nonreward_choice 1], 'box', 'on');
k = 0;
for n = 1:nPorts
    ap_mat          =     trialspxmat_nonreward_choice{n};
    t_mat             =     tspkmat_nonreward_choice{n};
    xx_all = [];
    yy_all = [];
    xxmt_all = [];
    yymt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:size(ap_mat, 2)
        imt = mt_nonreward_choice{n}(i);
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxmt = -imt;
        % plot trigger stimulus FPs_premature_cent_out
%         itrigger = FP_late_cent_out{n}(i);
        %     plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxmt_all = [xxmt_all, xxmt, xxmt, NaN];
        yymt_all = [yymt_all, yy2, NaN];

        k = k+1;
    end

    line(xx_all, yy_all, 'color', c_port{n}, 'linewidth', .2)
    line(xxmt_all, yymt_all, 'color', c_cent_out, 'linewidth', 1)
end

xline(0, 'Color', c_reward, 'LineWidth', 1, 'Alpha', 1);
title('Non-reward', 'fontsize', 7, 'fontweight','bold')
axis off

yshift_row5    =      yshift_row4 + 0.5 + ntrial_nonreward_choice*rasterheight;
% Add information
uicontrol('Style','text','Units','centimeters','Position',[x_pos-.25 yshift_row5 4 1],...
    'string', 'D. Reward-related', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 8,'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','Left');


%% Align to InitIn






%% Align to InitOut







% % % % % % % % % % % % 
% % % % % % % % % % % % %% Release PSTHs
% % % % % % % % % % % % % Release-related PSTHs
% % % % % % % % % % % % width = 6*sum(CentOutTimeDomain)/sum(CentInTimeDomain);
% % % % % % % % % % % % yshift_row1 = 1;
% % % % % % % % % % % % ha_cent_out_psth =  axes('unit', 'centimeters', 'position', [9.25 yshift_row1 width 2], 'nextplot', 'add', ...
% % % % % % % % % % % %     'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)]);
% % % % % % % % % % % % yshift_row2 = yshift_row1+2+0.25;
% % % % % % % % % % % % 
% % % % % % % % % % % % for i =1:nFPs
% % % % % % % % % % % %     hplot_cent_out(i) = plot(ts_cent_out{i}, psth_cent_out_correct{i}, 'color', c_FP(i, :),  'linewidth', 1.5);
% % % % % % % % % % % %     FRMax = max([FRMax max(psth_cent_out_correct{i})]);
% % % % % % % % % % % % %     disp(FRMax)
% % % % % % % % % % % % end
% % % % % % % % % % % % axis 'auto y'
% % % % % % % % % % % % hline_cent_out = line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1);
% % % % % % % % % % % % xlabel('Time from release (ms)')
% % % % % % % % % % % % ylabel ('Spks per s')
% % % % % % % % % % % % 
% % % % % % % % % % % % % error PSTHs
% % % % % % % % % % % % ha_cent_out_psth_error =  axes('unit', 'centimeters', 'position', [9.25 yshift_row2 width 2], 'nextplot', 'add', ...
% % % % % % % % % % % %     'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'xticklabel',[]);
% % % % % % % % % % % % yshift_row3 = yshift_row2 +2+0.25;
% % % % % % % % % % % % if  size(trialspxmat_premature_cent_out, 2)>3
% % % % % % % % % % % %     plot(ts_premature_cent_out, psth_premature_cent_out, 'color', c_premature, 'linewidth', 1.5)
% % % % % % % % % % % % %     FRMax = max([FRMax max(psth_premature_cent_out)]);
% % % % % % % % % % % % %     disp(FRMax)
% % % % % % % % % % % % end
% % % % % % % % % % % % if  size(trialspxmat_late_cent_out, 2)>3
% % % % % % % % % % % %     plot(ts_late_cent_out, psth_late_cent_out, 'color', c_late, 'linewidth', 1.5)
% % % % % % % % % % % % %     FRMax = max([FRMax max(psth_late_cent_out)]);
% % % % % % % % % % % % %     disp(FRMax)
% % % % % % % % % % % % end
% % % % % % % % % % % % axis 'auto y'
% % % % % % % % % % % % hline_cent_out_error =line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1);
% % % % % % % % % % % % 
% % % % % % % % % % % % % Make raster plot
% % % % % % % % % % % % % Plot spike raster of correct trials (all FPs)
% % % % % % % % % % % % ntrials_cent_out = 0;
% % % % % % % % % % % % nFP_ij = zeros(1, nFPs);
% % % % % % % % % % % % for i =1:nFPs
% % % % % % % % % % % %     nFP_ij(i) = size(trialspxmat_cent_out{i}, 2);
% % % % % % % % % % % %     ntrials_cent_out = ntrials_cent_out + nFP_ij(i);
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [9.25 yshift_row3 width ntrials_cent_out*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add',...
% % % % % % % % % % % %     'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrials_cent_out 1], 'box', 'on');
% % % % % % % % % % % % yshift_row4 = yshift_row3+ntrials_cent_out*rasterheight+0.5;
% % % % % % % % % % % % % Paint the foreperiod
% % % % % % % % % % % % n_start = 1;
% % % % % % % % % % % % k=0;
% % % % % % % % % % % % for m =1:nFPs
% % % % % % % % % % % %     ap_mat = trialspxmat_cent_out{m};
% % % % % % % % % % % %     t_mat = tspkmat_cent_out{m};
% % % % % % % % % % % %     rt = rt_cent_out_sorted{m};
% % % % % % % % % % % %     mFP = TargetFPs(m);
% % % % % % % % % % % %     xx_all = [];
% % % % % % % % % % % %     yy_all = [];
% % % % % % % % % % % %     xxrt_all = [];
% % % % % % % % % % % %     yyrt_all = [];
% % % % % % % % % % % %     x_portin = [];
% % % % % % % % % % % %     y_portin = [];
% % % % % % % % % % % %     for i =1:nFP_ij(m)
% % % % % % % % % % % %         irt = rt(i); % time from foreperiod to release
% % % % % % % % % % % %         xx = t_mat(ap_mat(:, i)==1);
% % % % % % % % % % % %         yy1 = [0 0.8]-k;
% % % % % % % % % % % %         yy2 = [0 1]-k;
% % % % % % % % % % % % 
% % % % % % % % % % % %         % paint foreperiod
% % % % % % % % % % % %         plotshaded([-irt-mFP -irt]-n_start, [-k -k; 1-k 1-k], c_trigger);
% % % % % % % % % % % % 
% % % % % % % % % % % %         for i_xx = 1:length(xx)
% % % % % % % % % % % %             xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %             yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %         end
% % % % % % % % % % % %         xxrt_all = [xxrt_all, -irt-mFP, -irt-mFP, NaN];
% % % % % % % % % % % %         yyrt_all = [yyrt_all, yy2, NaN];
% % % % % % % % % % % % 
% % % % % % % % % % % %         % port access time
% % % % % % % % % % % %         itrelease =t_cent_out_correct{m}(i);
% % % % % % % % % % % %         it_choice = t_choice - itrelease;
% % % % % % % % % % % %         it_choice = it_choice(it_choice>=-CentOutTimeDomain(1) & it_choice<=CentOutTimeDomain(2));
% % % % % % % % % % % %         if ~isempty(it_choice)
% % % % % % % % % % % %             it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %             x_portin = [x_portin, it_choice];
% % % % % % % % % % % %             y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];            
% % % % % % % % % % % %         end
% % % % % % % % % % % %         k = k+1;
% % % % % % % % % % % %     end
% % % % % % % % % % % %     n_start = n_start - nFP_ij(m);
% % % % % % % % % % % % 
% % % % % % % % % % % %     line(xx_all, yy_all, 'color', c_FP(m, :), 'linewidth', 1);
% % % % % % % % % % % %     line(xxrt_all, yyrt_all, 'color', c_cent_in, 'linewidth', 1.5);
% % % % % % % % % % % %     scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
% % % % % % % % % % % % end
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1);
% % % % % % % % % % % % title('Correct', 'fontsize', 7);
% % % % % % % % % % % % axis off
% % % % % % % % % % % % 
% % % % % % % % % % % % % Premature release raster plot
% % % % % % % % % % % % ntrial_premature = size(trialspxmat_premature_cent_out, 2); % number of trials
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [9.25 yshift_row4 width ntrial_premature*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add',...
% % % % % % % % % % % %     'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_premature 1], 'box', 'on');
% % % % % % % % % % % % yshift_row5    =      yshift_row4 + 0.5 + ntrial_premature*rasterheight;
% % % % % % % % % % % % ap_mat          =     trialspxmat_premature_cent_out;
% % % % % % % % % % % % t_mat             =     tspkmat_premature_cent_out;
% % % % % % % % % % % % k =0;
% % % % % % % % % % % % xx_all = [];
% % % % % % % % % % % % yy_all = [];
% % % % % % % % % % % % x_predur_all = [];
% % % % % % % % % % % % y_predur_all = [];
% % % % % % % % % % % % x_portin = [];
% % % % % % % % % % % % y_portin = [];
% % % % % % % % % % % % for i =1:size(ap_mat, 2)
% % % % % % % % % % % %     ipredur = premature_duration_cent_out(i);
% % % % % % % % % % % %     iFP = FPs_premature_cent_out(i);
% % % % % % % % % % % %     xx =  t_mat(ap_mat(:, i)==1);
% % % % % % % % % % % %     yy1 = [0 0.8]-k;
% % % % % % % % % % % %     yy = [0 1]-k;
% % % % % % % % % % % % 
% % % % % % % % % % % %     x_predur_all = [x_predur_all, -ipredur, -ipredur, NaN];
% % % % % % % % % % % %     y_predur_all = [y_predur_all, yy, NaN];
% % % % % % % % % % % % 
% % % % % % % % % % % %     % paint foreperiod
% % % % % % % % % % % %     plotshaded([-ipredur -ipredur+iFP], [-k -k; 1-k 1-k], c_trigger);
% % % % % % % % % % % % 
% % % % % % % % % % % %     for i_xx = 1:length(xx)
% % % % % % % % % % % %         xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %         yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %     end    
% % % % % % % % % % % % 
% % % % % % % % % % % %     % plot port poke time
% % % % % % % % % % % %     it_choice = t_choice - t_premature_cent_out(i);
% % % % % % % % % % % %     it_choice = it_choice(it_choice>=-CentOutTimeDomain(1) & it_choice<=CentOutTimeDomain(2));
% % % % % % % % % % % %     if ~isempty(it_choice)
% % % % % % % % % % % %         it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %         x_portin = [x_portin, it_choice];
% % % % % % % % % % % %         y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];         
% % % % % % % % % % % %     end
% % % % % % % % % % % %     k = k+1;
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % line(x_predur_all, y_predur_all, 'color', c_cent_in, 'linewidth', 1.5);
% % % % % % % % % % % % line(xx_all, yy_all, 'color', c_premature, 'linewidth', 1)
% % % % % % % % % % % % scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
% % % % % % % % % % % % 
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1)
% % % % % % % % % % % % title('Premature', 'fontsize', 7)
% % % % % % % % % % % % axis off
% % % % % % % % % % % % 
% % % % % % % % % % % % % Late response raster plot
% % % % % % % % % % % % ntrial_late = size(trialspxmat_late_cent_out, 2); % number of trials
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [9.25 yshift_row5  width ntrial_late*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add',...
% % % % % % % % % % % %     'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_late 1], 'box', 'on');
% % % % % % % % % % % % yshift_row6    =      yshift_row5 + 0.5 + ntrial_late*rasterheight;
% % % % % % % % % % % % ap_mat          =     trialspxmat_late_cent_out;
% % % % % % % % % % % % t_mat             =     tspkmat_late_cent_out;
% % % % % % % % % % % % k =0;
% % % % % % % % % % % % xx_all = [];
% % % % % % % % % % % % yy_all = [];
% % % % % % % % % % % % x_latedur_all = [];
% % % % % % % % % % % % y_latedur_all = [];
% % % % % % % % % % % % x_portin = [];
% % % % % % % % % % % % y_portin = [];
% % % % % % % % % % % % for i =1:size(ap_mat, 2)
% % % % % % % % % % % %     ilatedur =late_duration_cent_out(i);
% % % % % % % % % % % %     iFP = FPs_late_cent_out(i);
% % % % % % % % % % % %     xx =  t_mat(ap_mat(:, i)==1);
% % % % % % % % % % % %     yy1 = [0 0.8]-k;
% % % % % % % % % % % %     yy = [0 1]-k;
% % % % % % % % % % % %     % paint foreperiod
% % % % % % % % % % % %     plotshaded([-ilatedur -ilatedur+iFP], [-k -k; 1-k 1-k], c_trigger);
% % % % % % % % % % % %     x_latedur_all = [x_latedur_all, -ilatedur, -ilatedur, NaN];
% % % % % % % % % % % %     y_latedur_all = [y_latedur_all, yy, NaN];
% % % % % % % % % % % % 
% % % % % % % % % % % %     for i_xx = 1:length(xx)
% % % % % % % % % % % %         xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %         yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %     end  
% % % % % % % % % % % % 
% % % % % % % % % % % %     % plot port poke time
% % % % % % % % % % % %     it_choice = t_choice - t_late_cent_out(i);
% % % % % % % % % % % %     it_choice = it_choice(it_choice>=-CentInTimeDomain(1) & it_choice<=CentInTimeDomain(2));
% % % % % % % % % % % %     if ~isempty(it_choice)
% % % % % % % % % % % %         it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %         x_portin = [x_portin, it_choice];
% % % % % % % % % % % %         y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];  
% % % % % % % % % % % %     end
% % % % % % % % % % % %     k = k+1;
% % % % % % % % % % % % end
% % % % % % % % % % % % line(x_latedur_all, y_latedur_all, 'color', c_cent_in, 'linewidth', 1.5);
% % % % % % % % % % % % line(xx_all, yy_all, 'color', c_late, 'linewidth', 1)
% % % % % % % % % % % % scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
% % % % % % % % % % % % 
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1)
% % % % % % % % % % % % title('Late', 'fontsize', 7)
% % % % % % % % % % % % axis off
% % % % % % % % % % % % 
% % % % % % % % % % % % % Add information
% % % % % % % % % % % % uicontrol('Style','text','Units','centimeters','Position',[4.75 yshift_row6 width+1 1],...
% % % % % % % % % % % %     'string', 'B. CentOut-related', ...
% % % % % % % % % % % %     'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],...
% % % % % % % % % % % %     'HorizontalAlignment','Left');
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % %% Reward
% % % % % % % % % % % % col3 = 13;
% % % % % % % % % % % % width = 6*sum(ChoiceTimeDomain)/sum(CentInTimeDomain);
% % % % % % % % % % % % ha_poke =  axes('unit', 'centimeters', 'position', [col3 yshift_row1 width 2], 'nextplot', 'add', ...
% % % % % % % % % % % %     'xlim', [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)]);
% % % % % % % % % % % % for i =1:nFPs
% % % % % % % % % % % %     plot(ts_reward_choice{i}, psth_reward_choice{i}, 'color', c_FP(i, :), 'linewidth', 1.5);
% % % % % % % % % % % %     FRMax = max([FRMax max(psth_reward_choice{i})]);
% % % % % % % % % % % % %     disp(FRMax)
% % % % % % % % % % % % end
% % % % % % % % % % % % % also add non-rewarded pokes
% % % % % % % % % % % % plot(ts_nonreward_choice, psth_nonreward_choice, 'color', [0.6 0.6 0.6], 'linewidth', .5);
% % % % % % % % % % % % xlabel('Time from rewarded/nonrewarded poke (ms)')
% % % % % % % % % % % % ylabel ('Spks per s')
% % % % % % % % % % % % axis 'auto y'
% % % % % % % % % % % % hline_poke = line([0 0], get(gca, 'ylim'), 'color', c_reward, 'linewidth', 1);
% % % % % % % % % % % % % FRMax = max([FRMax max(psth_nonreward_pokes)]);
% % % % % % % % % % % % % disp(FRMax)
% % % % % % % % % % % % % Raster plot
% % % % % % % % % % % % 
% % % % % % % % % % % % % Make raster plot
% % % % % % % % % % % % % Plot spike raster of correct trials (all FPs)
% % % % % % % % % % % % ntrials_rewardpoke = 0;
% % % % % % % % % % % % nFP_ij = zeros(1, nFPs);
% % % % % % % % % % % % for i =1:nFPs
% % % % % % % % % % % %     nFP_ij(i) = size(trialspxmat_reward_choice{i}, 2);
% % % % % % % % % % % %     ntrials_rewardpoke = ntrials_rewardpoke + nFP_ij(i);
% % % % % % % % % % % % end;
% % % % % % % % % % % % 
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [col3 yshift_row2 width ntrials_rewardpoke*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add', 'xlim',  [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], ...
% % % % % % % % % % % %     'ylim', [-ntrials_rewardpoke 1], 'box', 'on', 'xticklabel', []);
% % % % % % % % % % % % yshift_row3new = yshift_row2 + ntrials_rewardpoke*rasterheight + 0.5;
% % % % % % % % % % % % 
% % % % % % % % % % % % % Paint the foreperiod
% % % % % % % % % % % % k = 0;
% % % % % % % % % % % % for m =1:nFPs
% % % % % % % % % % % %     ap_mat = trialspxmat_reward_choice{m};
% % % % % % % % % % % %     t_mat = tspkmat_reward_choice{m};
% % % % % % % % % % % %     move=move_time{m};
% % % % % % % % % % % %     mFP = TargetFPs(m);
% % % % % % % % % % % %     xx_all = [];
% % % % % % % % % % % %     yy_all = [];
% % % % % % % % % % % %     x_move_all = [];
% % % % % % % % % % % %     y_move_all = [];
% % % % % % % % % % % %     x_portin = [];
% % % % % % % % % % % %     y_portin = [];
% % % % % % % % % % % %     for i =1:nFP_ij(m)
% % % % % % % % % % % %         xx =  t_mat(ap_mat(:, i)==1);
% % % % % % % % % % % %         yy1 = [0 0.8]-k;
% % % % % % % % % % % %         yy2 = [0 1]-k;
% % % % % % % % % % % %         imov = -move(i);
% % % % % % % % % % % % 
% % % % % % % % % % % %         for i_xx = 1:length(xx)
% % % % % % % % % % % %             xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %             yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %         end 
% % % % % % % % % % % % 
% % % % % % % % % % % %         x_move_all = [x_move_all, imov, imov, NaN];
% % % % % % % % % % % %         y_move_all = [y_move_all, yy2, NaN];
% % % % % % % % % % % % 
% % % % % % % % % % % %         % plot port poke time
% % % % % % % % % % % %         itreward =t_reward_choice{m}(i);
% % % % % % % % % % % %         it_choice = t_choice-itreward;
% % % % % % % % % % % %         it_choice = it_choice(it_choice>=-ChoiceTimeDomain(1) & it_choice<=ChoiceTimeDomain(2));
% % % % % % % % % % % %         if ~isempty(it_choice)
% % % % % % % % % % % %             it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %             x_portin = [x_portin, it_choice];
% % % % % % % % % % % %             y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];            
% % % % % % % % % % % %         end
% % % % % % % % % % % %         k = k+1;
% % % % % % % % % % % %     end
% % % % % % % % % % % %     line(xx_all, yy_all, 'color', c_FP(m,:), 'linewidth', 1)
% % % % % % % % % % % %     line(x_move_all, y_move_all, 'color', c_cent_out, 'linewidth', 1)
% % % % % % % % % % % %     scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
% % % % % % % % % % % % end
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1);
% % % % % % % % % % % % title('Correct', 'fontsize', 7);
% % % % % % % % % % % % axis off
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % Raster plot for unrewarded pokes
% % % % % % % % % % % % % use at most 50 events
% % % % % % % % % % % % if size(trialspxmat_nonreward_choice, 2)>50
% % % % % % % % % % % %     plot_ind = sort(randperm(size(trialspxmat_nonreward_choice, 2), 50));
% % % % % % % % % % % %     trialspxmat_nonreward_pokes_plot = trialspxmat_nonreward_choice(:, plot_ind);
% % % % % % % % % % % % else
% % % % % % % % % % % %     trialspxmat_nonreward_pokes_plot = trialspxmat_nonreward_choice;
% % % % % % % % % % % %     plot_ind = 1:size(trialspxmat_nonreward_choice, 2);
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % ntrial_nonrewardpoke = size(trialspxmat_nonreward_pokes_plot, 2);
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [col3 yshift_row3new width ntrial_nonrewardpoke*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add', 'xlim',  [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], ...
% % % % % % % % % % % %     'ylim', [-ntrial_nonrewardpoke 1], 'box', 'on', 'xticklabel', []);
% % % % % % % % % % % % yshift_row4new = yshift_row3new+0.5+ntrial_nonrewardpoke*rasterheight;
% % % % % % % % % % % % k =0;
% % % % % % % % % % % % move_time_nonreward_plot             =                          mt_nonreward_choice(plot_ind);
% % % % % % % % % % % % t_nonreward_pokes_plot                  =                          t_nonreward_choice(plot_ind);
% % % % % % % % % % % % 
% % % % % % % % % % % % xx_all = [];
% % % % % % % % % % % % yy_all = [];
% % % % % % % % % % % % x_move_all = [];
% % % % % % % % % % % % y_move_all = [];
% % % % % % % % % % % % x_portin = [];
% % % % % % % % % % % % y_portin = [];
% % % % % % % % % % % % for i =1:ntrial_nonrewardpoke
% % % % % % % % % % % %     if isempty(find(isnan(trialspxmat_nonreward_pokes_plot(:, i)), 1))
% % % % % % % % % % % %         xx =  tspkmat_nonreward_choice(trialspxmat_nonreward_pokes_plot(:, i)==1);
% % % % % % % % % % % %         yy1 = [0 0.8]-k;
% % % % % % % % % % % %         yy2 = [0 1]-k;
% % % % % % % % % % % %         imov = -move_time_nonreward_plot(i);
% % % % % % % % % % % % 
% % % % % % % % % % % %         for i_xx = 1:length(xx)
% % % % % % % % % % % %             xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %             yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %         end 
% % % % % % % % % % % % 
% % % % % % % % % % % %         x_move_all = [x_move_all, imov, imov, NaN];
% % % % % % % % % % % %         y_move_all = [y_move_all, yy2, NaN];
% % % % % % % % % % % % 
% % % % % % % % % % % %         % plot port poke time
% % % % % % % % % % % %         itreward =t_nonreward_pokes_plot(i);
% % % % % % % % % % % %         it_choice = t_choice-itreward;
% % % % % % % % % % % %         it_choice = it_choice(it_choice>=-ChoiceTimeDomain(1) & it_choice<=ChoiceTimeDomain(2));
% % % % % % % % % % % %         if ~isempty(it_choice)
% % % % % % % % % % % %             it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %             x_portin = [x_portin, it_choice];
% % % % % % % % % % % %             y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];             
% % % % % % % % % % % %         end
% % % % % % % % % % % %     end
% % % % % % % % % % % %     k = k+1;
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % line(xx_all, yy_all, 'color', [0.6 0.6 0.6], 'linewidth', 1)
% % % % % % % % % % % % line(x_move_all, y_move_all, 'color', c_cent_out, 'linewidth', 1)
% % % % % % % % % % % % scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
% % % % % % % % % % % % 
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color', c_reward, 'linewidth', 1)
% % % % % % % % % % % % axis off
% % % % % % % % % % % % title('Nonrewarded pokes', 'fontname', 'dejavu sans', 'fontsize', 7)
% % % % % % % % % % % % 
% % % % % % % % % % % % % Add information  13.5 3+0.5 6 ntrial4*rasterheight
% % % % % % % % % % % % uicontrol('Style','text','Units','centimeters','Position',[col3-0.5 yshift_row4new 5 1.75],...
% % % % % % % % % % % %     'string', 'C. Rewarded/Nonrewarded poke-related activity', ...
% % % % % % % % % % % %     'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
% % % % % % % % % % % %     'HorizontalAlignment','Left');
% % % % % % % % % % % % 
% % % % % % % % % % % % yshift_row5new = yshift_row4new+1.75+1.5;
% % % % % % % % % % % % 
% % % % % % % % % % % % %% plot pre-press activity vs trial num or time
% % % % % % % % % % % % 
% % % % % % % % % % % % ha10=axes('unit', 'centimeters', 'position', [col3 yshift_row5new 5 2.5], ...
% % % % % % % % % % % %     'nextplot', 'add', 'xlim', [min(t_correct_cent_in_all/1000) max(t_correct_cent_in_all/1000)]);
% % % % % % % % % % % % 
% % % % % % % % % % % % ind_prepress = find(tspkmat_cent_in_all<0);
% % % % % % % % % % % % spkmat_prepress =  trialspxmat_cent_in_all(ind_prepress, :);
% % % % % % % % % % % % dur_prepress = abs(tspkmat_cent_in_all(ind_prepress(1)))/1000; % total time
% % % % % % % % % % % % 
% % % % % % % % % % % % rate_prepress = sum(spkmat_prepress, 1)/dur_prepress; % spk rate across time
% % % % % % % % % % % % plot(ha10, t_correct_cent_in_all/1000, rate_prepress, 'k', 'marker', 'o', 'markersize', 3, 'linestyle', 'none');
% % % % % % % % % % % % 
% % % % % % % % % % % % % linear regression
% % % % % % % % % % % % Pfit = polyfit(t_correct_cent_in_all/1000,rate_prepress,1);
% % % % % % % % % % % % yfit = Pfit(1)*t_correct_cent_in_all/1000+Pfit(2);
% % % % % % % % % % % % plot(t_correct_cent_in_all/1000,yfit,'r:', 'linewidth', 1.5);
% % % % % % % % % % % % 
% % % % % % % % % % % % xlabel('Time (s)')
% % % % % % % % % % % % ylabel('Spk rate (Hz)')
% % % % % % % % % % % % 
% % % % % % % % % % % % yshift_row6new = yshift_row5new+3;
% % % % % % % % % % % % % Add information  13.5 3+0.5 6 ntrial4*rasterheight
% % % % % % % % % % % % uicontrol('Style','text','Units','centimeters','Position',[col3-0.5 yshift_row6new 4 0.5],...
% % % % % % % % % % % %     'string', 'D.  Activity vs time', ...
% % % % % % % % % % % %     'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
% % % % % % % % % % % %     'HorizontalAlignment','Left');
% % % % % % % % % % % % 
% % % % % % % % % % % % fig_height = max([fig_height, yshift_row6new+1]);
% % % % % % % % % % % % 
% % % % % % % % % % % % %% plot trigger-related activity
% % % % % % % % % % % % col4 = 20;
% % % % % % % % % % % % width = 6*sum(TriggerTimeDomain)/sum(CentInTimeDomain);
% % % % % % % % % % % % 
% % % % % % % % % % % % ha_trigger =  axes('unit', 'centimeters', 'position', [col4 yshift_row1 width 2], 'nextplot', 'add', ...
% % % % % % % % % % % %     'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)]);
% % % % % % % % % % % % for i=1:nFPs
% % % % % % % % % % % %     plot(ts_trigger_correct{i}, psth_trigger_correct{i}, 'color', c_FP(i, :), 'linewidth', 1.5);
% % % % % % % % % % % %     FRMax = max([FRMax max(psth_trigger_correct{i})]);
% % % % % % % % % % % % %     disp(FRMax)
% % % % % % % % % % % % end
% % % % % % % % % % % % plot(ts_late_trigger, psth_late_trigger, 'color', c_late, 'linewidth', 1.5)
% % % % % % % % % % % % xlabel('Time from trigger stimulus (ms)')
% % % % % % % % % % % % ylabel ('Spks per s')
% % % % % % % % % % % % 
% % % % % % % % % % % % % FRMax = max([FRMax max(psth_late_trigger)]);
% % % % % % % % % % % % % disp(FRMax)
% % % % % % % % % % % % xlim = max(get(gca, 'xlim'));
% % % % % % % % % % % % axis 'auto y'
% % % % % % % % % % % % 
% % % % % % % % % % % % % raster plot of trigger-related activity
% % % % % % % % % % % % ntrials_trigger = 0;
% % % % % % % % % % % % nFP_ij = zeros(1, nFPs);
% % % % % % % % % % % % for i =1:nFPs
% % % % % % % % % % % %     nFP_ij(i) = size(trialspxmat_trigger_correct{i}, 2);
% % % % % % % % % % % %     ntrials_trigger = ntrials_trigger + nFP_ij(i);
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [col4 yshift_row2 width ntrials_trigger*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add',...
% % % % % % % % % % % %     'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'ylim', [-ntrials_trigger 1], 'box', 'on');
% % % % % % % % % % % % yshift_row3 = yshift_row2+ntrials_trigger*rasterheight+0.5;
% % % % % % % % % % % % % Paint the foreperiod
% % % % % % % % % % % % k=0;
% % % % % % % % % % % % for m =1:nFPs
% % % % % % % % % % % %     ap_mat = trialspxmat_trigger_correct{m};
% % % % % % % % % % % %     t_mat = tspkmat_trigger_correct{m};
% % % % % % % % % % % %     rt = rt_trigger_correct{m};
% % % % % % % % % % % %     mFP = TargetFPs(m);
% % % % % % % % % % % %     xx_all = [];
% % % % % % % % % % % %     yy_all = [];
% % % % % % % % % % % %     xxrt_all = [];
% % % % % % % % % % % %     yyrt_all = [];
% % % % % % % % % % % %     x_portin = [];
% % % % % % % % % % % %     y_portin = [];
% % % % % % % % % % % %     for i =1:nFP_ij(m)
% % % % % % % % % % % %         irt = rt(i); % time from trigger to release
% % % % % % % % % % % %         xx = t_mat(ap_mat(:, i)==1);
% % % % % % % % % % % %         yy1 = [0 0.8]-k;
% % % % % % % % % % % %         yy = [0 1]-k;
% % % % % % % % % % % %         % paint foreperiod
% % % % % % % % % % % %         plotshaded([-mFP 0], [-k -k; 1-k 1-k], c_trigger);
% % % % % % % % % % % % 
% % % % % % % % % % % %         for i_xx = 1:length(xx)
% % % % % % % % % % % %             xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %             yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %         end 
% % % % % % % % % % % % 
% % % % % % % % % % % %         % plot release time
% % % % % % % % % % % %         xxrt_all = [xxrt_all, irt, irt, NaN];
% % % % % % % % % % % %         yyrt_all = [yyrt_all, yy, NaN];
% % % % % % % % % % % % 
% % % % % % % % % % % %         % port access time
% % % % % % % % % % % %         it_trigger =t_trigger_correct{m}(i);
% % % % % % % % % % % %         it_choice = t_choice - it_trigger;
% % % % % % % % % % % %         it_choice = it_choice(it_choice>=-TriggerTimeDomain(1) & it_choice<=TriggerTimeDomain(2));
% % % % % % % % % % % %         if ~isempty(it_choice)
% % % % % % % % % % % %             it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %             x_portin = [x_portin, it_choice];
% % % % % % % % % % % %             y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
% % % % % % % % % % % %         end
% % % % % % % % % % % %         k = k+1;
% % % % % % % % % % % %     end
% % % % % % % % % % % %     line(xx_all, yy_all, 'color', c_FP(m, :), 'linewidth', 1);
% % % % % % % % % % % %     line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5);
% % % % % % % % % % % %     scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
% % % % % % % % % % % % end
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color', c_trigger, 'linewidth', 1);
% % % % % % % % % % % % title('Correct', 'fontsize', 7);
% % % % % % % % % % % % axis off
% % % % % % % % % % % % 
% % % % % % % % % % % % % trigger following late FP
% % % % % % % % % % % % ntrials_trigger_late = size(trialspxmat_late_trigger, 2);
% % % % % % % % % % % % axes('unit', 'centimeters', 'position', [col4 yshift_row3 width ntrials_trigger_late*rasterheight],...
% % % % % % % % % % % %     'nextplot', 'add', 'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'ylim', [-ntrials_trigger_late 1], ...
% % % % % % % % % % % %     'box', 'on', 'xticklabel', []);
% % % % % % % % % % % % yshift_row4 = yshift_row3+ntrials_trigger_late*rasterheight+0.5;
% % % % % % % % % % % % k =0;
% % % % % % % % % % % % xx_all = [];
% % % % % % % % % % % % yy_all = [];
% % % % % % % % % % % % xxrt_all = [];
% % % % % % % % % % % % yyrt_all = [];
% % % % % % % % % % % % x_portin = [];
% % % % % % % % % % % % y_portin = [];
% % % % % % % % % % % % for i =1:ntrials_trigger_late
% % % % % % % % % % % %     xx =  tspkmat_late_trigger(trialspxmat_late_trigger(:, i)==1);
% % % % % % % % % % % %     iFP = FP_trigger_late(i);
% % % % % % % % % % % %     yy1 = [0 0.8]-k;
% % % % % % % % % % % %     yy2 = [0 1]-k;
% % % % % % % % % % % %     % paint foreperiod
% % % % % % % % % % % %     plotshaded([-iFP 0], [-k -k; 1-k 1-k], c_trigger);
% % % % % % % % % % % % 
% % % % % % % % % % % %     for i_xx = 1:length(xx)
% % % % % % % % % % % %         xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
% % % % % % % % % % % %         yy_all = [yy_all, yy1, NaN];
% % % % % % % % % % % %     end 
% % % % % % % % % % % % 
% % % % % % % % % % % %     % plot release time
% % % % % % % % % % % %     irt = RT_trigger_late(i);
% % % % % % % % % % % %     xxrt_all = [xxrt_all, irt, irt, NaN];
% % % % % % % % % % % %     yyrt_all = [yyrt_all, yy2, NaN];
% % % % % % % % % % % %     % plot port poke time
% % % % % % % % % % % %     itrigger = t_trigger_late(i);
% % % % % % % % % % % %     it_choice = t_choice-itrigger;
% % % % % % % % % % % %     it_choice = it_choice(it_choice>=-TriggerTimeDomain(1) & it_choice<=TriggerTimeDomain(2));
% % % % % % % % % % % %     if ~isempty(it_choice)
% % % % % % % % % % % %         it_choice = reshape(it_choice,1,[]);
% % % % % % % % % % % %         x_portin = [x_portin, it_choice];
% % % % % % % % % % % %         y_portin = [y_portin, (0.4-k)*ones(1,length(it_choice))];
% % % % % % % % % % % %     end
% % % % % % % % % % % %     k = k+1;
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % line(xx_all, yy_all, 'color', c_late, 'linewidth', 1);
% % % % % % % % % % % % line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5);
% % % % % % % % % % % % scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
% % % % % % % % % % % % 
% % % % % % % % % % % % line([0 0], get(gca, 'ylim'), 'color',c_trigger, 'linewidth', 1)
% % % % % % % % % % % % title('late', 'fontname', 'dejavu sans', 'fontsize', 7)
% % % % % % % % % % % % axis off
% % % % % % % % % % % % 
% % % % % % % % % % % % % Add information
% % % % % % % % % % % % uicontrol('Style','text','Units','centimeters','Position',[col4-0.5  yshift_row4 5 1.2],...
% % % % % % % % % % % %     'string', 'E. Trigger-related activity', ...
% % % % % % % % % % % %     'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
% % % % % % % % % % % %     'HorizontalAlignment','Left');
% % % % % % % % % % % % yshift_row5=yshift_row4+1.2+1.2;
% % % % % % % % % % % % 

FRrange = [0 FRMax*1.1];
cellfun(@(x) set(x, 'ylim', FRrange), ha_cent_in_psth);
set(ha_cent_in_psth_late, 'ylim', FRrange);
set(ha_cent_in_psth_premature, 'ylim', FRrange);
set(ha_cent_in_psth_probe, 'ylim', FRrange);
cellfun(@(x) set(x, 'ylim', FRrange), ha_cent_out_psth);
set(ha_cent_out_psth_late, 'ylim', FRrange);
set(ha_cent_out_psth_premature, 'ylim', FRrange);
set(ha_cent_out_psth_probe, 'ylim', FRrange);
cellfun(@(x) set(x, 'ylim', FRrange), ha_choice_psth);
set(ha_nonreward_choice_psth, 'ylim', FRrange);
cellfun(@(x) set(x, 'ylim', FRrange), ha_trigger_psth);
set(ha_trigger_psth_late, 'ylim', FRrange);

%% plot spks

col5=10.5;
thiscolor = [0 0 0];
Lspk = size(r.Units.SpikeTimes(ku).wave, 2);
ha0=axes('unit', 'centimeters', 'position', [col5 yshift_row7+.5 1.5 1.5], ...
    'nextplot', 'add', 'xlim', [0 Lspk], 'ytick', -500:100:200, 'xticklabel', []);
set(ha0, 'nextplot', 'add');
ylabel('uV')
allwaves = r.Units.SpikeTimes(ku).wave/4;
if size(allwaves, 1)>100
    nplot = randperm(size(allwaves, 1), 100);
else
    nplot=1:size(allwaves, 1);
end
wave2plot = allwaves(nplot, :);
plot(1:Lspk, wave2plot, 'color', [0.8 .8 0.8]);
plot(1:Lspk, mean(allwaves, 1), 'color', thiscolor, 'linewidth', 2)
axis([0 Lspk min(wave2plot(:)) max(wave2plot(:))])
set (gca, 'ylim', [min(mean(allwaves, 1))*1.25 max(mean(allwaves, 1))*1.25])
axis tight
line([30 60], min(get(gca, 'ylim')), 'color', 'k', 'linewidth', 2.5)
PSTH.SpikeWave = mean(allwaves, 1);
% plot autocorrelation
kutime = round(r.Units.SpikeTimes(ku).timings);
kutime = kutime(kutime>0);
kutime2 = zeros(1, max(kutime));
kutime2(kutime)=1;
[c, lags] = xcorr(kutime2, 100); % max lag 100 ms
c(lags==0)=0;

ha00= axes('unit', 'centimeters', 'position', [col5+1.5+1 yshift_row7+.5 2 1.5], 'nextplot', 'add', 'xlim', [-25 25]);
if median(c)>1
    set(ha00, 'nextplot', 'add', 'xtick', -50:10:50, 'ytick', [0 median(c)]);
else
    set(ha00, 'nextplot', 'add', 'xtick', -50:10:50, 'ytick', [0 1], 'ylim', [0 1]);
end

switch r.Units.SpikeNotes(ku, 3)
    case 1
        title(['#' num2str(ku) '(Ch ' num2str(r.Units.SpikeNotes(ku, 1)) ' | unit' num2str(r.Units.SpikeNotes(ku, 2))  ' | SU'], 'fontsize', 7);
    case 2
        title(['#' num2str(ku) '(Ch ' num2str(r.Units.SpikeNotes(ku, 1))  ' | unit' num2str(r.Units.SpikeNotes(ku, 2))  ' | MU'], 'fontsize', 7);
    otherwise
end

PSTH.AutoCorrelation = {lags, c};

hbar = bar(lags, c);
set(hbar, 'facecolor', 'k');
xlabel('Lag(ms)')

yshift_row8 = yshift_row7+3;
%% Plot all waveforms if it is a polytrode

if isfield(r.Units.SpikeTimes(ku), 'wave_mean')
    ha_wave_poly = axes('unit', 'centimeters', 'position', [col5 yshift_row8 4 1.5], 'nextplot', 'add');
    wave_form = r.Units.SpikeTimes(ku).wave_mean/4;
    PSTH.SpikeWaveMean = wave_form;
    n_chs = size(wave_form, 1); % number of channels
    ch_selected = 1:n_chs;
    if n_chs > 32
        n_chs = 32;
        ch_largest = r.Units.SpikeNotes(ku,1);
        if ch_largest < n_chs/2
            ch_selected = 1:n_chs;
        elseif ch_largest > size(wave_form, 1) - n_chs/2
            ch_selected = size(wave_form, 1)-n_chs+1:size(wave_form, 1);
        else
            ch_selected = ch_largest-15:ch_largest+16;
        end
    end
    n_sample = size(wave_form, 2); % sample size per spike
    n_cols = 8;
    n_rows = n_chs/n_cols;
    max_x = 0;
    colors = [25, 167, 206]/255;
    if n_rows<1
        n_rows=1;
    end
    v_sep = 100;

    t_wave_all = [];
    wave_all = [];
    for i = 1:n_rows
        for j = 1:n_cols
            k = j+(i-1)*n_cols;
            wave_k = wave_form(ch_selected(k), :)+v_sep*(i-1);
            t_wave = (1:n_sample)+n_sample*(j-1)+4;
            t_wave_all = [t_wave_all, t_wave, NaN];
            wave_all = [wave_all, wave_k, NaN];
            max_x = max([max_x, max(t_wave)]);
        end
    end
    plot(ha_wave_poly, t_wave_all, wave_all, 'linewidth', 1, 'color', colors);

    set(ha_wave_poly, 'xlim', [0 max_x], 'ylim', [-400  v_sep*(n_rows-1)+200]);
    axis off
    axis tight

    yshift_row9 = yshift_row8+1.5;
else
    yshift_row9 = yshift_row8;
end

uicontrol('Style','text','Units','centimeters','Position',[col5 yshift_row9 4 1.5],...
    'string', 'E. Spike waveform and autocorrelation', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 8,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
    'HorizontalAlignment','Left');
fig_height = max([fig_height, yshift_row9+2]);
% change the height of the figure
set(gcf, 'position', [2 2 17 fig_height])

toc;
if strcmpi(ToSave,'on')
    % save to a folder
    anm_name        =     r.BehaviorClass.Subject;
    session              =     r.BehaviorClass.Session;
    
    PSTH.ANM_Session = {anm_name, session};
    thisFolder = fullfile(pwd, 'Fig');
    if ~exist(thisFolder, 'dir')
        mkdir(thisFolder)
    end
    tosavename2= fullfile(thisFolder, anm_name+"_"+session+"_Ch"+num2str(ch)+"_Unit"+num2str(unit_no));
    print (gcf,'-dpng', tosavename2, '-r300')
%     exportgraphics(gcf, tosavename2+".png");
    
    % save PSTH as well save(psth_new_name, 'PSTHOut');
    save(tosavename2+".mat", 'PSTH')
    
    try
        tic
        % C:\Users\jiani\OneDrive\00_Work\03_Projects
        thisFolder = fullfile(findonedrive, '00_Work', '03_Projects', '05_Physiology', 'Data', 'UnitsCollection', anm_name, session);
        if ~exist(thisFolder, 'dir')
            mkdir(thisFolder)
        end
        copyfile([tosavename2 '.png'], thisFolder)
        copyfile([tosavename2 '.mat'], thisFolder)
    
        toc
    end
end