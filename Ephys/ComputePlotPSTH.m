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
c_port     = {c.PortL, c.PortR};

% For PSTH and raster plots
c_cent_in  = [5 191 219] / 255;
c_trigger  = [242 182 250] / 255;
c_cent_out = [87 108 188] / 255;
c_reward   = [164 208 164] / 255;

TargetFPs = unique(PSTHOut.CentIn.FP{1});
nFPs = length(TargetFPs);
if nFPs == 2
    c_FP = [192, 127, 0; 76, 61, 61]/255;
else
    c_FP = [255, 217, 90; 192, 127, 0; 76, 61, 61]/255;
end

Ports = unique(PSTHOut.CentIn.Port{1});
nPorts = length(Ports);

nSort = nFPs * nPorts;

c_premature = [0.9 0.4 0.1];
c_late      = [0.6 0.6 0.6];
printsize   = [2 2 25 25];

%% PSTHs for CentIn and CentOut
params_cent_in.pre      = 5000; % take a longer pre-press activity so we can compute z score easily later.
params_cent_in.post     = CentInTimeDomain(2);
params_cent_in.binwidth = 20;

t_cent_in = PSTHOut.CentIn.Time{end};
[psth_cent_in_all, ts_cent_in_all, trialspxmat_cent_in_all, tspkmat_cent_in_all, t_correct_cent_in_all, ~] = jpsth(r.Units.SpikeTimes(ku).timings, t_cent_in, params_cent_in);

psth_cent_in_all = smoothdata (psth_cent_in_all, 'gaussian', 5);
PSTH.CentInAll   = {psth_cent_in_all, ts_cent_in_all, trialspxmat_cent_in_all, tspkmat_cent_in_all,  t_correct_cent_in_all};

% CentIn PSTH (corrected, sorted)
params.pre      = CentInTimeDomain(1);
params.post     = CentInTimeDomain(2);
params.binwidth = 20;

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
        [psth_cent_in_correct{i,j}, ts_cent_in{i,j}, trialspxmat_cent_in{i,j}, tspkmat_cent_in{i,j}, t_cent_in_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params);
        psth_cent_in_correct{i,j} = smoothdata(psth_cent_in_correct{i,j}, 'gaussian', 5);

        rt_cent_in_sorted{i,j} = PSTHOut.CentIn.RT_Correct{i,j};
        rt_cent_in_sorted{i,j} = rt_cent_in_sorted{i,j}(ind);

        PSTH.CentIn{i,j} = {psth_cent_in_correct{i,j}, ts_cent_in{i,j}, trialspxmat_cent_in{i,j}, tspkmat_cent_in{i,j}, t_cent_in_correct{i,j}, rt_cent_in_sorted{i,j}};
    end
end
PSTH.CentInLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT'};

% CentOut PSTH (corrected, sorted)
params.pre      = CentOutTimeDomain(1);
params.post     = CentOutTimeDomain(2);
params.binwidth = 20;

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
        [psth_cent_out_correct{i,j}, ts_cent_out{i,j}, trialspxmat_cent_out{i,j}, tspkmat_cent_out{i,j}, t_cent_out_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params);
        psth_cent_out_correct{i,j} = smoothdata (psth_cent_out_correct{i,j}, 'gaussian', 5);

        rt_cent_out_sorted{i,j} = PSTHOut.CentIn.RT_Correct{i,j};
        rt_cent_out_sorted{i,j} = rt_cent_out_sorted{i,j}(ind);

        PSTH.CentOut{i,j} = {psth_cent_out_correct{i,j}, ts_cent_out{i,j}, trialspxmat_cent_out{i,j}, tspkmat_cent_out{i,j}, t_cent_out_correct{i,j}, rt_cent_out_sorted{i,j}};
    end
end
PSTH.CentOutLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT'};

% premature cent_in PSTH
t_premature_cent_in = PSTHOut.CentIn.Time{nSort+1};
[psth_premature_cent_in, ts_premature_cent_in, trialspxmat_premature_cent_in, tspkmat_premature_cent_in, t_premature_cent_in, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_premature_cent_in, params);
psth_premature_cent_in = smoothdata (psth_premature_cent_in, 'gaussian', 5);

FPs_premature_cent_in = PSTHOut.CentIn.FP{2};
FPs_premature_cent_in = FPs_premature_cent_in(ind);
Ports_premature_cent_in = PSTHOut.CentIn.Port{2};
Ports_premature_cent_in = Ports_premature_cent_in(ind);

HD_premature_cent_in = PSTHOut.CentIn.HoldDur.Premature;
HD_premature_cent_in = HD_premature_cent_in(ind);
PSTH.PrematureCentIn = {psth_premature_cent_in, ts_premature_cent_in, trialspxmat_premature_cent_in, tspkmat_premature_cent_in, t_premature_cent_in, HD_premature_cent_in, FPs_premature_cent_in, Ports_premature_cent_in};
PSTH.PrematureLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP', 'Port'};

% premature cent_out PSTH
t_premature_cent_out = PSTHOut.CentOut.Time{nSort+1};
[psth_premature_cent_out, ts_premature_cent_out, trialspxmat_premature_cent_out, tspkmat_premature_cent_out, t_premature_cent_out, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_premature_cent_out, params);
psth_premature_cent_out = smoothdata (psth_premature_cent_out, 'gaussian', 5);

FPs_premature_cent_out = PSTHOut.CentOut.FP{2};
FPs_premature_cent_out = FPs_premature_cent_out(ind);
Ports_premature_cent_out = PSTHOut.CentOut.Port{2};
Ports_premature_cent_out = Ports_premature_cent_out(ind);
HD_premature_cent_out = PSTHOut.CentOut.HoldDur.Premature;
HD_premature_cent_out = HD_premature_cent_out(ind);

PSTH.PrematureCentOut = {psth_premature_cent_out, ts_premature_cent_out, trialspxmat_premature_cent_out, tspkmat_premature_cent_out, t_premature_cent_out, HD_premature_cent_out, FPs_premature_cent_out, Ports_premature_cent_out};
PSTH.PrematureLabels  = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP', 'Port'};

% late cent_in PSTH
t_late_cent_in = PSTHOut.CentIn.Time{nSort+2};
[psth_late_cent_in, ts_late_cent_in, trialspxmat_late_cent_in, tspkmat_late_cent_in, t_late_cent_in, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_late_cent_in, params);
psth_late_cent_in  = smoothdata (psth_late_cent_in, 'gaussian', 5);

FPs_late_cent_in   = PSTHOut.CentIn.FP{3};
FPs_late_cent_in   = FPs_late_cent_in(ind);
Ports_late_cent_in = PSTHOut.CentIn.Port{3};
Ports_late_cent_in = Ports_late_cent_in(ind);
HD_late_cent_in    = PSTHOut.CentIn.HoldDur.Late;
HD_late_cent_in    = HD_late_cent_in(ind);

PSTH.LateCentIn = {psth_late_cent_in, ts_late_cent_in, trialspxmat_late_cent_in, tspkmat_late_cent_in, t_late_cent_in, HD_late_cent_in, FPs_late_cent_in, Ports_late_cent_in};
PSTH.LateLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP', 'Port'};

% late cent_out PSTH
t_late_cent_out = PSTHOut.CentOut.Time{nSort+2};
[psth_late_cent_out, ts_late_cent_out, trialspxmat_late_cent_out, tspkmat_late_cent_out, t_late_cent_out, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_late_cent_out, params);
psth_late_cent_out  = smoothdata (psth_late_cent_out, 'gaussian', 5);

FPs_late_cent_out   = PSTHOut.CentOut.FP{3};
FPs_late_cent_out   = FPs_late_cent_out(ind);
Ports_late_cent_out = PSTHOut.CentOut.Port{3};
Ports_late_cent_out = Ports_late_cent_out(ind);
HD_late_cent_out    = PSTHOut.CentOut.HoldDur.Late;
HD_late_cent_out    = HD_late_cent_out(ind);

PSTH.LateCentOut = {psth_late_cent_out, ts_late_cent_out, trialspxmat_late_cent_out, tspkmat_late_cent_out, t_late_cent_out, HD_late_cent_out, FPs_late_cent_out, Ports_late_cent_out};
PSTH.LateLabels  = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP', 'Port'};

% probe cent_in PSTH
t_probe_cent_in = PSTHOut.CentIn.Time{nSort+2};
[psth_probe_cent_in, ts_probe_cent_in, trialspxmat_probe_cent_in, tspkmat_probe_cent_in, t_probe_cent_in, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_probe_cent_in, params);
psth_probe_cent_in  = smoothdata (psth_probe_cent_in, 'gaussian', 5);

FPs_probe_cent_in   = PSTHOut.CentIn.FP{3};
FPs_probe_cent_in   = FPs_probe_cent_in(ind);
Ports_probe_cent_in = PSTHOut.CentIn.Port{3};
Ports_probe_cent_in = Ports_probe_cent_in(ind);
HD_probe_cent_in    = PSTHOut.CentIn.HoldDur.Probe;
HD_probe_cent_in    = HD_probe_cent_in(ind);

PSTH.ProbeCentIn = {psth_probe_cent_in, ts_probe_cent_in, trialspxmat_probe_cent_in, tspkmat_probe_cent_in, t_probe_cent_in, HD_probe_cent_in, FPs_probe_cent_in, Ports_probe_cent_in};
PSTH.ProbeLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP', 'Port'};

% probe cent_out PSTH
t_probe_cent_out = PSTHOut.CentOut.Time{nSort+2};
[psth_probe_cent_out, ts_probe_cent_out, trialspxmat_probe_cent_out, tspkmat_probe_cent_out, t_probe_cent_out, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_probe_cent_out, params);
psth_probe_cent_out  = smoothdata (psth_probe_cent_out, 'gaussian', 5);

FPs_probe_cent_out   = PSTHOut.CentOut.FP{3};
FPs_probe_cent_out   = FPs_probe_cent_out(ind);
Ports_probe_cent_out = PSTHOut.CentOut.Port{3};
Ports_probe_cent_out = Ports_probe_cent_out(ind);
HD_probe_cent_out    = PSTHOut.CentOut.HoldDur.Probe;
HD_probe_cent_out    = HD_probe_cent_out(ind);

PSTH.ProbeCentOut = {psth_probe_cent_out, ts_probe_cent_out, trialspxmat_probe_cent_out, tspkmat_probe_cent_out, t_probe_cent_out, HD_probe_cent_out, FPs_probe_cent_out, Ports_probe_cent_out};
PSTH.ProbeLabels  = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'HoldDuration', 'FP', 'Port'};

% use t_reward_poke and move_time to construct reward_poke PSTH
% reward PSTH
params.pre  = ChoiceTimeDomain(1);
params.post = ChoiceTimeDomain(2);

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
        [psth_reward_choice{i,j}, ts_reward_choice{i,j}, trialspxmat_reward_choice{i,j}, tspkmat_reward_choice{i,j}, t_reward_choice{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params);
        psth_reward_choice{i,j} = smoothdata (psth_reward_choice{i,j}, 'gaussian', 5);

        mt_reward_choice{i,j} = PSTHOut.ChoiceIn.RewardPoke.Move_Time{i,j}(ind);
        PSTH.RewardChoice{i,j} = {psth_reward_choice{i,j}, ts_reward_choice{i,j}, trialspxmat_reward_choice{i,j}, tspkmat_reward_choice{i,j}, t_reward_choice{i,j}, mt_reward_choice{i,j}};
    end
end
PSTH.ChoiceLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'MoveTime'};

% nonreward choice PSTH
t_nonreward_choice    = PSTHOut.ChoiceIn.NonrewardPoke.Time;
mt_nonreward_choice   = PSTHOut.ChoiceIn.NonrewardPoke.Move_Time;
port_nonreward_choice = PSTHOut.ChoiceIn.NonrewardPoke.PortChosen;
[psth_nonreward_choice, ts_nonreward_choice, trialspxmat_nonreward_choice, tspkmat_nonreward_choice, t_nonreward_choice, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_nonreward_choice, params);
psth_nonreward_choice = smoothdata (psth_nonreward_choice, 'gaussian', 5);
mt_nonreward_choice   = mt_nonreward_choice(ind);
port_nonreward_choice = port_nonreward_choice(ind);
PSTH.NonrewardChoice  = {psth_nonreward_choice, ts_nonreward_choice, trialspxmat_nonreward_choice, tspkmat_nonreward_choice, t_nonreward_choice, mt_nonreward_choice, port_nonreward_choice};

% trigger PSTH
params.pre  = TriggerTimeDomain(1);
params.post = TriggerTimeDomain(2);

% late response
t_trigger_late    = PSTHOut.Triggers.Time{end};
RT_trigger_late   = PSTHOut.Triggers.RT{end};
FP_trigger_late   = PSTHOut.Triggers.FP{end};
Port_trigger_late = PSTHOut.Triggers.Port{end};

[psth_late_trigger, ts_late_trigger, trialspxmat_late_trigger, tspkmat_late_trigger, t_trigger_late, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_trigger_late, params);
psth_late_trigger = smoothdata (psth_late_trigger, 'gaussian', 5);

RT_trigger_late   = RT_trigger_late(ind);
FP_trigger_late   = FP_trigger_late(ind);
Port_trigger_late = Port_trigger_late(ind);

PSTH.TriggerLate   = {psth_late_trigger, ts_late_trigger, trialspxmat_late_trigger, tspkmat_late_trigger, t_trigger_late, RT_trigger_late, FP_trigger_late, Port_trigger_late};
PSTH.TriggerLabels = {'PSTH', 'tPSTH', 'SpikeMat', 'tSpikeMat', 'tEvents', 'RT', 'FP', 'Port'};

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
        [psth_trigger_correct{i,j}, ts_trigger_correct{i,j}, trialspxmat_trigger_correct{i,j}, tspkmat_trigger_correct{i,j}, t_trigger_correct{i,j}, ind] = jpsth(r.Units.SpikeTimes(ku).timings, t_ij, params);
        psth_trigger_correct{i,j} = smoothdata (psth_trigger_correct{i,j}, 'gaussian', 5);

        rt_trigger_correct{i,j} = PSTHOut.Triggers.RT{ind_ij};
        rt_trigger_correct{i,j} = rt_trigger_correct{i,j}(ind);
        PSTH.Triggers{i,j} = {psth_trigger_correct{i,j}, ts_trigger_correct{i,j}, trialspxmat_trigger_correct{i,j}, tspkmat_trigger_correct{i,j}, t_trigger_correct{i,j}, rt_trigger_correct{i,j}, TargetFPs(i), Ports(j)};
    end
end

%% Plot raster and spks
figure();
set(gcf, 'unit', 'centimeters', 'position', printsize, 'paperpositionmode', 'auto' ,'color', 'w')
% PSTH of correct trials
yshift_row1 = 1;
ha_cent_in_psth =  axes('unit', 'centimeters', 'position', [1.25 yshift_row1 6 2], 'nextplot', 'add', 'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)]);
yshift_row2 = yshift_row1+2+0.25;
hplot_cent_in= zeros(1, nFPs);
FRMax = 3;
for i =1:nFPs
    hplot_cent_in(i) = plot(ts_cent_in{i}, psth_cent_in_correct{i}, 'color', c_FP(i, :),  'linewidth', 1.5);
    FRMax = max([FRMax max(psth_cent_in_correct{i})]);
%     disp(FRMax)
end
axis 'auto y'
xlabel('Time from press (ms)')
ylabel ('Spks per s')

% PSTH of error trials (premature and late)
ha_cent_in_psth_error =  axes('unit', 'centimeters', 'position', [1.25 yshift_row2 6 2], 'nextplot', 'add',...
    'xlim',  [-CentInTimeDomain(1) CentInTimeDomain(2)], 'xticklabel', []);
yshift_row3 = yshift_row2 +2+0.25;
% plot premature and late as well
if  size(trialspxmat_premature_cent_in, 2)>3
    plot(ts_premature_cent_in, psth_premature_cent_in, 'color', c_premature, 'linewidth',1.5);
%     FRMax = max([FRMax max(psth_premature_cent_in)]);
%      disp(FRMax)
end
if  size(trialspxmat_late_cent_in, 2)>3
    plot(ts_late_cent_in, psth_late_cent_in, 'color', c_late, 'linewidth', 1.5)
%     FRMax = max([FRMax max(psth_late_cent_in)]);
%     disp(FRMax)
end
axis 'auto y'
hline_cent_in_error = line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1);

% make raster plot  750 ms FP
if num2str(length(t_cent_in))>200
    rasterheight = 0.01;
elseif num2str(length(t_cent_in))>100
    rasterheight = 0.02;
else
    rasterheight = 0.04;
end

% Plot spike raster of correct trials (all FPs)
ntrials_cent_in = 0;
nFP_i = zeros(1, nFPs);
t_portin =  PSTHOut.Pokes.Time;
for i =1:nFPs
    nFP_i(i) = size(trialspxmat_cent_in{i}, 2);
    ntrials_cent_in = ntrials_cent_in + nFP_i(i);
end
axes('unit', 'centimeters', 'position', [1.25 yshift_row3 6 ntrials_cent_in*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrials_cent_in 1], 'box', 'on');
yshift_row4 = yshift_row3+ntrials_cent_in*rasterheight+0.5;
% Paint the foreperiod
k=0;
for m = 1:nFPs
    ap_mat = trialspxmat_cent_in{m};
    t_mat = tspkmat_cent_in{m};
    rt = rt_cent_in_sorted{m};
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:nFP_i(m)
        irt = rt(i); % time from foreperiod to release
        xx = t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        xxrt = irt+TargetFPs(m);
        plotshaded([0 TargetFPs(m)],[-k -k; 1-k 1-k], c_trigger);

        if isempty(find(isnan(ap_mat(:, i)), 1))
            for i_xx = 1:length(xx)
                xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
                yy_all = [yy_all, yy1, NaN];
            end
            xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
            yyrt_all = [yyrt_all, yy2, NaN];
        end
        % port access time
        itpress =t_cent_in_correct{m}(i);
        i_portin = t_portin - itpress;
        i_portin = i_portin(i_portin>=-CentInTimeDomain(1) & i_portin<=CentInTimeDomain(2));
        if ~isempty(i_portin)
            i_portin = reshape(i_portin,1,[]);
            x_portin = [x_portin, i_portin];
            y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];
        end
        k = k+1;
    end
    line(xx_all, yy_all, 'color', c_FP(m, :), 'linewidth', 1);
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5);
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
end

line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1);
title('Correct', 'fontsize', 7, 'fontweight','bold');
axis off

% Premature press raster plot
ntrial_premature = size(trialspxmat_premature_cent_in, 2); % number of trials
axes('unit', 'centimeters', 'position', [1.25 yshift_row4 6 ntrial_premature*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrial_premature 1], 'box', 'on');
yshift_row5    =      yshift_row4 + 0.5 + ntrial_premature*rasterheight;
ap_mat          =     trialspxmat_premature_cent_in;
t_mat             =     tspkmat_premature_cent_in;
k =0;
xx_all = [];
yy_all = [];
xxrt_all = [];
yyrt_all = [];
x_portin = [];
y_portin = [];
for i =1:size(ap_mat, 2)
    ipredur = HD_premature_cent_in(i);
    xx =  t_mat(ap_mat(:, i)==1);
    yy1 = [0 0.8]-k;
    yy2 = [0 1]-k;
    xxrt = ipredur;
    % plot trigger stimulus FPs_premature_cent_in
    itrigger = FPs_premature_cent_in(i);
    plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

    for i_xx = 1:length(xx)
        xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
        yy_all = [yy_all, yy1, NaN];
    end
    xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
    yyrt_all = [yyrt_all, yy2, NaN];

    % plot port poke time
    i_portin = t_portin - t_premature_cent_in(i);
    i_portin = i_portin(i_portin>=-CentInTimeDomain(1) & i_portin<=CentInTimeDomain(2));
    if ~isempty(i_portin)
        i_portin = reshape(i_portin,1,[]);
        x_portin = [x_portin, i_portin];
        y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];        
    end
    k = k+1;
end

line(xx_all, yy_all, 'color', c_premature, 'linewidth', 1)
line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5)
scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')

line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1)
title('Premature', 'fontsize', 7, 'fontweight','bold')
axis off

% Late response raster plot
ntrial_late = size(trialspxmat_late_cent_in, 2); % number of trials
axes('unit', 'centimeters', 'position', [1.25 yshift_row5  6 ntrial_late*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentInTimeDomain(1) CentInTimeDomain(2)], 'ylim', [-ntrial_late 1], 'box', 'on');
yshift_row6             =      yshift_row5 + 0.5 + ntrial_late*rasterheight;
ap_mat          =     trialspxmat_late_cent_in;
t_mat             =     tspkmat_late_cent_in;
k =0;
xx_all = [];
yy_all = [];
xxrt_all = [];
yyrt_all = [];
x_portin = [];
y_portin = [];
for i =1:size(ap_mat, 2)
    ilatedur =HD_late_cent_in(i);
    xx =  t_mat(ap_mat(:, i)==1);
    yy1 = [0 0.8]-k;
    yy2 = [0 1]-k;
    xxrt = ilatedur;
    % plot trigger stimulus FPs_premature_cent_in
    itrigger = FPs_late_cent_in(i);
    plotshaded([0 itrigger], [-k -k; 1-k 1-k], c_trigger);

    for i_xx = 1:length(xx)
        xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
        yy_all = [yy_all, yy1, NaN];
    end
    xxrt_all = [xxrt_all, xxrt, xxrt, NaN];
    yyrt_all = [yyrt_all, yy2, NaN];

    % plot port poke time
    i_portin = t_portin - t_late_cent_in(i);
    i_portin = i_portin(i_portin>=-CentInTimeDomain(1) & i_portin<=CentInTimeDomain(2));
    if ~isempty(i_portin)
        i_portin = reshape(i_portin,1,[]);
        x_portin = [x_portin, i_portin];
        y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];
    end
    k = k+1;
end

line(xx_all, yy_all, 'color', c_late, 'linewidth', 1)
line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5)
scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')

line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1)
title('Late', 'fontsize', 7, 'fontweight','bold')
axis off

% this is the position of last panel
% Add information
uicontrol('Style','text','Units','centimeters','Position',[0.5 yshift_row6  6 1],...
    'string', 'A. Press-related activity', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','Left');

yshift_row7=yshift_row6+1.25;
ch = r.Units.SpikeNotes(ku, 1);
unit_no = r.Units.SpikeNotes(ku, 2);

if size(r.Units.SpikeNotes, 2) == 4
    cluster_id = r.Units.SpikeNotes(ku, 4);
    uicontrol('style', 'text', 'units', 'centimeters', 'position', [1 yshift_row7 6 1.2],...
        'string', (['Unit #' num2str(ku) ' (Ch ' num2str(ch) ' | UnitOnCh ' num2str(unit_no) ' | ' 'Kilosort cluster ' num2str(cluster_id) ')']),...
        'BackgroundColor','w', 'fontsize', 10, 'fontweight','bold',  'FontName','Dejavu Sans')
else
    cluster_id = [];
    uicontrol('style', 'text', 'units', 'centimeters', 'position', [1 yshift_row7 6 1.2],...
        'string', (['Unit #' num2str(ku) ' (' num2str(ch) ' | ' num2str(unit_no) ')']),...
        'BackgroundColor','w', 'fontsize', 10, 'fontweight','bold',  'FontName','Dejavu Sans')
end
uicontrol('style', 'text', 'units', 'centimeters', 'position', [1 yshift_row7+1.2 4 0.5],...
    'string', ([r.Meta(1).Subject ' ' r.Meta(1).DateTime(1:11)]), 'BackgroundColor','w',...
    'fontsize', 10, 'fontweight', 'bold',  'FontName','Dejavu Sans')

fig_height = yshift_row7+2;

%% Release PSTHs
% Release-related PSTHs
width = 6*sum(CentOutTimeDomain)/sum(CentInTimeDomain);
yshift_row1 = 1;
ha_cent_out_psth =  axes('unit', 'centimeters', 'position', [8.25 yshift_row1 width 2], 'nextplot', 'add', ...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)]);
yshift_row2 = yshift_row1+2+0.25;

for i =1:nFPs
    hplot_cent_out(i) = plot(ts_cent_out{i}, psth_cent_out_correct{i}, 'color', c_FP(i, :),  'linewidth', 1.5);
    FRMax = max([FRMax max(psth_cent_out_correct{i})]);
%     disp(FRMax)
end
axis 'auto y'
hline_cent_out = line([0 0], get(gca, 'ylim'), 'color', c_cent_in, 'linewidth', 1);
xlabel('Time from release (ms)')
ylabel ('Spks per s')

% error PSTHs
ha_cent_out_psth_error =  axes('unit', 'centimeters', 'position', [8.25 yshift_row2 width 2], 'nextplot', 'add', ...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'xticklabel',[]);
yshift_row3 = yshift_row2 +2+0.25;
if  size(trialspxmat_premature_cent_out, 2)>3
    plot(ts_premature_cent_out, psth_premature_cent_out, 'color', c_premature, 'linewidth', 1.5)
%     FRMax = max([FRMax max(psth_premature_cent_out)]);
%     disp(FRMax)
end
if  size(trialspxmat_late_cent_out, 2)>3
    plot(ts_late_cent_out, psth_late_cent_out, 'color', c_late, 'linewidth', 1.5)
%     FRMax = max([FRMax max(psth_late_cent_out)]);
%     disp(FRMax)
end
axis 'auto y'
hline_cent_out_error =line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1);

% Make raster plot
% Plot spike raster of correct trials (all FPs)
ntrials_cent_out = 0;
nFP_i = zeros(1, nFPs);
for i =1:nFPs
    nFP_i(i) = size(trialspxmat_cent_out{i}, 2);
    ntrials_cent_out = ntrials_cent_out + nFP_i(i);
end

axes('unit', 'centimeters', 'position', [8.25 yshift_row3 width ntrials_cent_out*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrials_cent_out 1], 'box', 'on');
yshift_row4 = yshift_row3+ntrials_cent_out*rasterheight+0.5;
% Paint the foreperiod
n_start = 1;
k=0;
for m =1:nFPs
    ap_mat = trialspxmat_cent_out{m};
    t_mat = tspkmat_cent_out{m};
    rt = rt_cent_out_sorted{m};
    mFP = TargetFPs(m);
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:nFP_i(m)
        irt = rt(i); % time from foreperiod to release
        xx = t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;

        % paint foreperiod
        plotshaded([-irt-mFP -irt]-n_start, [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end
        xxrt_all = [xxrt_all, -irt-mFP, -irt-mFP, NaN];
        yyrt_all = [yyrt_all, yy2, NaN];

        % port access time
        itrelease =t_cent_out_correct{m}(i);
        i_portin = t_portin - itrelease;
        i_portin = i_portin(i_portin>=-CentOutTimeDomain(1) & i_portin<=CentOutTimeDomain(2));
        if ~isempty(i_portin)
            i_portin = reshape(i_portin,1,[]);
            x_portin = [x_portin, i_portin];
            y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];            
        end
        k = k+1;
    end
    n_start = n_start - nFP_i(m);

    line(xx_all, yy_all, 'color', c_FP(m, :), 'linewidth', 1);
    line(xxrt_all, yyrt_all, 'color', c_cent_in, 'linewidth', 1.5);
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
end
line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1);
title('Correct', 'fontsize', 7);
axis off

% Premature release raster plot
ntrial_premature = size(trialspxmat_premature_cent_out, 2); % number of trials
axes('unit', 'centimeters', 'position', [8.25 yshift_row4 width ntrial_premature*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_premature 1], 'box', 'on');
yshift_row5    =      yshift_row4 + 0.5 + ntrial_premature*rasterheight;
ap_mat          =     trialspxmat_premature_cent_out;
t_mat             =     tspkmat_premature_cent_out;
k =0;
xx_all = [];
yy_all = [];
x_predur_all = [];
y_predur_all = [];
x_portin = [];
y_portin = [];
for i =1:size(ap_mat, 2)
    ipredur = premature_duration_cent_out(i);
    iFP = FPs_premature_cent_out(i);
    xx =  t_mat(ap_mat(:, i)==1);
    yy1 = [0 0.8]-k;
    yy = [0 1]-k;

    x_predur_all = [x_predur_all, -ipredur, -ipredur, NaN];
    y_predur_all = [y_predur_all, yy, NaN];

    % paint foreperiod
    plotshaded([-ipredur -ipredur+iFP], [-k -k; 1-k 1-k], c_trigger);

    for i_xx = 1:length(xx)
        xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
        yy_all = [yy_all, yy1, NaN];
    end    

    % plot port poke time
    i_portin = t_portin - t_premature_cent_out(i);
    i_portin = i_portin(i_portin>=-CentOutTimeDomain(1) & i_portin<=CentOutTimeDomain(2));
    if ~isempty(i_portin)
        i_portin = reshape(i_portin,1,[]);
        x_portin = [x_portin, i_portin];
        y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];         
    end
    k = k+1;
end

line(x_predur_all, y_predur_all, 'color', c_cent_in, 'linewidth', 1.5);
line(xx_all, yy_all, 'color', c_premature, 'linewidth', 1)
scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')

line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1)
title('Premature', 'fontsize', 7)
axis off

% Late response raster plot
ntrial_late = size(trialspxmat_late_cent_out, 2); % number of trials
axes('unit', 'centimeters', 'position', [8.25 yshift_row5  width ntrial_late*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-CentOutTimeDomain(1) CentOutTimeDomain(2)], 'ylim', [-ntrial_late 1], 'box', 'on');
yshift_row6    =      yshift_row5 + 0.5 + ntrial_late*rasterheight;
ap_mat          =     trialspxmat_late_cent_out;
t_mat             =     tspkmat_late_cent_out;
k =0;
xx_all = [];
yy_all = [];
x_latedur_all = [];
y_latedur_all = [];
x_portin = [];
y_portin = [];
for i =1:size(ap_mat, 2)
    ilatedur =late_duration_cent_out(i);
    iFP = FPs_late_cent_out(i);
    xx =  t_mat(ap_mat(:, i)==1);
    yy1 = [0 0.8]-k;
    yy = [0 1]-k;
    % paint foreperiod
    plotshaded([-ilatedur -ilatedur+iFP], [-k -k; 1-k 1-k], c_trigger);
    x_latedur_all = [x_latedur_all, -ilatedur, -ilatedur, NaN];
    y_latedur_all = [y_latedur_all, yy, NaN];

    for i_xx = 1:length(xx)
        xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
        yy_all = [yy_all, yy1, NaN];
    end  

    % plot port poke time
    i_portin = t_portin - t_late_cent_out(i);
    i_portin = i_portin(i_portin>=-CentInTimeDomain(1) & i_portin<=CentInTimeDomain(2));
    if ~isempty(i_portin)
        i_portin = reshape(i_portin,1,[]);
        x_portin = [x_portin, i_portin];
        y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];  
    end
    k = k+1;
end
line(x_latedur_all, y_latedur_all, 'color', c_cent_in, 'linewidth', 1.5);
line(xx_all, yy_all, 'color', c_late, 'linewidth', 1)
scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')

line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1)
title('Late', 'fontsize', 7)
axis off

% Add information
uicontrol('Style','text','Units','centimeters','Position',[7.75 yshift_row6 width+1 1],...
    'string', 'B. Release-related', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','Left');

%% Reward
col3 = 13;
width = 6*sum(ChoiceTimeDomain)/sum(CentInTimeDomain);
ha_poke =  axes('unit', 'centimeters', 'position', [col3 yshift_row1 width 2], 'nextplot', 'add', ...
    'xlim', [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)]);
for i =1:nFPs
    plot(ts_reward_choice{i}, psth_reward_choice{i}, 'color', c_FP(i, :), 'linewidth', 1.5);
    FRMax = max([FRMax max(psth_reward_choice{i})]);
%     disp(FRMax)
end
% also add non-rewarded pokes
plot(ts_nonreward_choice, psth_nonreward_choice, 'color', [0.6 0.6 0.6], 'linewidth', .5);
xlabel('Time from rewarded/nonrewarded poke (ms)')
ylabel ('Spks per s')
axis 'auto y'
hline_poke = line([0 0], get(gca, 'ylim'), 'color', c_reward, 'linewidth', 1);
% FRMax = max([FRMax max(psth_nonreward_pokes)]);
% disp(FRMax)
% Raster plot

% Make raster plot
% Plot spike raster of correct trials (all FPs)
ntrials_rewardpoke = 0;
nFP_i = zeros(1, nFPs);
for i =1:nFPs
    nFP_i(i) = size(trialspxmat_reward_choice{i}, 2);
    ntrials_rewardpoke = ntrials_rewardpoke + nFP_i(i);
end;

axes('unit', 'centimeters', 'position', [col3 yshift_row2 width ntrials_rewardpoke*rasterheight],...
    'nextplot', 'add', 'xlim',  [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], ...
    'ylim', [-ntrials_rewardpoke 1], 'box', 'on', 'xticklabel', []);
yshift_row3new = yshift_row2 + ntrials_rewardpoke*rasterheight + 0.5;

% Paint the foreperiod
k = 0;
for m =1:nFPs
    ap_mat = trialspxmat_reward_choice{m};
    t_mat = tspkmat_reward_choice{m};
    move=move_time{m};
    mFP = TargetFPs(m);
    xx_all = [];
    yy_all = [];
    x_move_all = [];
    y_move_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:nFP_i(m)
        xx =  t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        imov = -move(i);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end 

        x_move_all = [x_move_all, imov, imov, NaN];
        y_move_all = [y_move_all, yy2, NaN];

        % plot port poke time
        itreward =t_reward_choice{m}(i);
        i_portin = t_portin-itreward;
        i_portin = i_portin(i_portin>=-ChoiceTimeDomain(1) & i_portin<=ChoiceTimeDomain(2));
        if ~isempty(i_portin)
            i_portin = reshape(i_portin,1,[]);
            x_portin = [x_portin, i_portin];
            y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];            
        end
        k = k+1;
    end
    line(xx_all, yy_all, 'color', c_FP(m,:), 'linewidth', 1)
    line(x_move_all, y_move_all, 'color', c_cent_out, 'linewidth', 1)
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')
end
line([0 0], get(gca, 'ylim'), 'color', c_cent_out, 'linewidth', 1);
title('Correct', 'fontsize', 7);
axis off


% Raster plot for unrewarded pokes
% use at most 50 events
if size(trialspxmat_nonreward_choice, 2)>50
    plot_ind = sort(randperm(size(trialspxmat_nonreward_choice, 2), 50));
    trialspxmat_nonreward_pokes_plot = trialspxmat_nonreward_choice(:, plot_ind);
else
    trialspxmat_nonreward_pokes_plot = trialspxmat_nonreward_choice;
    plot_ind = 1:size(trialspxmat_nonreward_choice, 2);
end

ntrial_nonrewardpoke = size(trialspxmat_nonreward_pokes_plot, 2);
axes('unit', 'centimeters', 'position', [col3 yshift_row3new width ntrial_nonrewardpoke*rasterheight],...
    'nextplot', 'add', 'xlim',  [-ChoiceTimeDomain(1) ChoiceTimeDomain(2)], ...
    'ylim', [-ntrial_nonrewardpoke 1], 'box', 'on', 'xticklabel', []);
yshift_row4new = yshift_row3new+0.5+ntrial_nonrewardpoke*rasterheight;
k =0;
move_time_nonreward_plot             =                          mt_nonreward_choice(plot_ind);
t_nonreward_pokes_plot                  =                          t_nonreward_choice(plot_ind);

xx_all = [];
yy_all = [];
x_move_all = [];
y_move_all = [];
x_portin = [];
y_portin = [];
for i =1:ntrial_nonrewardpoke
    if isempty(find(isnan(trialspxmat_nonreward_pokes_plot(:, i)), 1))
        xx =  tspkmat_nonreward_choice(trialspxmat_nonreward_pokes_plot(:, i)==1);
        yy1 = [0 0.8]-k;
        yy2 = [0 1]-k;
        imov = -move_time_nonreward_plot(i);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end 

        x_move_all = [x_move_all, imov, imov, NaN];
        y_move_all = [y_move_all, yy2, NaN];

        % plot port poke time
        itreward =t_nonreward_pokes_plot(i);
        i_portin = t_portin-itreward;
        i_portin = i_portin(i_portin>=-ChoiceTimeDomain(1) & i_portin<=ChoiceTimeDomain(2));
        if ~isempty(i_portin)
            i_portin = reshape(i_portin,1,[]);
            x_portin = [x_portin, i_portin];
            y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];             
        end
    end
    k = k+1;
end

line(xx_all, yy_all, 'color', [0.6 0.6 0.6], 'linewidth', 1)
line(x_move_all, y_move_all, 'color', c_cent_out, 'linewidth', 1)
scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')

line([0 0], get(gca, 'ylim'), 'color', c_reward, 'linewidth', 1)
axis off
title('Nonrewarded pokes', 'fontname', 'dejavu sans', 'fontsize', 7)

% Add information  13.5 3+0.5 6 ntrial4*rasterheight
uicontrol('Style','text','Units','centimeters','Position',[col3-0.5 yshift_row4new 5 1.75],...
    'string', 'C. Rewarded/Nonrewarded poke-related activity', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
    'HorizontalAlignment','Left');

yshift_row5new = yshift_row4new+1.75+1.5;

%% plot pre-press activity vs trial num or time

ha10=axes('unit', 'centimeters', 'position', [col3 yshift_row5new 5 2.5], ...
    'nextplot', 'add', 'xlim', [min(t_correct_cent_in_all/1000) max(t_correct_cent_in_all/1000)]);

ind_prepress = find(tspkmat_cent_in_all<0);
spkmat_prepress =  trialspxmat_cent_in_all(ind_prepress, :);
dur_prepress = abs(tspkmat_cent_in_all(ind_prepress(1)))/1000; % total time

rate_prepress = sum(spkmat_prepress, 1)/dur_prepress; % spk rate across time
plot(ha10, t_correct_cent_in_all/1000, rate_prepress, 'k', 'marker', 'o', 'markersize', 3, 'linestyle', 'none');

% linear regression
Pfit = polyfit(t_correct_cent_in_all/1000,rate_prepress,1);
yfit = Pfit(1)*t_correct_cent_in_all/1000+Pfit(2);
plot(t_correct_cent_in_all/1000,yfit,'r:', 'linewidth', 1.5);

xlabel('Time (s)')
ylabel('Spk rate (Hz)')

yshift_row6new = yshift_row5new+3;
% Add information  13.5 3+0.5 6 ntrial4*rasterheight
uicontrol('Style','text','Units','centimeters','Position',[col3-0.5 yshift_row6new 4 0.5],...
    'string', 'D.  Activity vs time', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
    'HorizontalAlignment','Left');

fig_height = max([fig_height, yshift_row6new+1]);

%% plot trigger-related activity
col4 = 20;
width = 6*sum(TriggerTimeDomain)/sum(CentInTimeDomain);

ha_trigger =  axes('unit', 'centimeters', 'position', [col4 yshift_row1 width 2], 'nextplot', 'add', ...
    'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)]);
for i=1:nFPs
    plot(ts_trigger_correct{i}, psth_trigger_correct{i}, 'color', c_FP(i, :), 'linewidth', 1.5);
    FRMax = max([FRMax max(psth_trigger_correct{i})]);
%     disp(FRMax)
end
plot(ts_late_trigger, psth_late_trigger, 'color', c_late, 'linewidth', 1.5)
xlabel('Time from trigger stimulus (ms)')
ylabel ('Spks per s')

% FRMax = max([FRMax max(psth_late_trigger)]);
% disp(FRMax)
xlim = max(get(gca, 'xlim'));
axis 'auto y'

% raster plot of trigger-related activity
ntrials_trigger = 0;
nFP_i = zeros(1, nFPs);
for i =1:nFPs
    nFP_i(i) = size(trialspxmat_trigger_correct{i}, 2);
    ntrials_trigger = ntrials_trigger + nFP_i(i);
end

axes('unit', 'centimeters', 'position', [col4 yshift_row2 width ntrials_trigger*rasterheight],...
    'nextplot', 'add',...
    'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'ylim', [-ntrials_trigger 1], 'box', 'on');
yshift_row3 = yshift_row2+ntrials_trigger*rasterheight+0.5;
% Paint the foreperiod
k=0;
for m =1:nFPs
    ap_mat = trialspxmat_trigger_correct{m};
    t_mat = tspkmat_trigger_correct{m};
    rt = rt_trigger_correct{m};
    mFP = TargetFPs(m);
    xx_all = [];
    yy_all = [];
    xxrt_all = [];
    yyrt_all = [];
    x_portin = [];
    y_portin = [];
    for i =1:nFP_i(m)
        irt = rt(i); % time from trigger to release
        xx = t_mat(ap_mat(:, i)==1);
        yy1 = [0 0.8]-k;
        yy = [0 1]-k;
        % paint foreperiod
        plotshaded([-mFP 0], [-k -k; 1-k 1-k], c_trigger);

        for i_xx = 1:length(xx)
            xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
            yy_all = [yy_all, yy1, NaN];
        end 

        % plot release time
        xxrt_all = [xxrt_all, irt, irt, NaN];
        yyrt_all = [yyrt_all, yy, NaN];

        % port access time
        it_trigger =t_trigger_correct{m}(i);
        i_portin = t_portin - it_trigger;
        i_portin = i_portin(i_portin>=-TriggerTimeDomain(1) & i_portin<=TriggerTimeDomain(2));
        if ~isempty(i_portin)
            i_portin = reshape(i_portin,1,[]);
            x_portin = [x_portin, i_portin];
            y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];
        end
        k = k+1;
    end
    line(xx_all, yy_all, 'color', c_FP(m, :), 'linewidth', 1);
    line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5);
    scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none');
end
line([0 0], get(gca, 'ylim'), 'color', c_trigger, 'linewidth', 1);
title('Correct', 'fontsize', 7);
axis off

% trigger following late FP
ntrials_trigger_late = size(trialspxmat_late_trigger, 2);
axes('unit', 'centimeters', 'position', [col4 yshift_row3 width ntrials_trigger_late*rasterheight],...
    'nextplot', 'add', 'xlim', [-TriggerTimeDomain(1) TriggerTimeDomain(2)], 'ylim', [-ntrials_trigger_late 1], ...
    'box', 'on', 'xticklabel', []);
yshift_row4 = yshift_row3+ntrials_trigger_late*rasterheight+0.5;
k =0;
xx_all = [];
yy_all = [];
xxrt_all = [];
yyrt_all = [];
x_portin = [];
y_portin = [];
for i =1:ntrials_trigger_late
    xx =  tspkmat_late_trigger(trialspxmat_late_trigger(:, i)==1);
    iFP = FP_trigger_late(i);
    yy1 = [0 0.8]-k;
    yy2 = [0 1]-k;
    % paint foreperiod
    plotshaded([-iFP 0], [-k -k; 1-k 1-k], c_trigger);

    for i_xx = 1:length(xx)
        xx_all = [xx_all, xx(i_xx), xx(i_xx), NaN];
        yy_all = [yy_all, yy1, NaN];
    end 

    % plot release time
    irt = RT_trigger_late(i);
    xxrt_all = [xxrt_all, irt, irt, NaN];
    yyrt_all = [yyrt_all, yy2, NaN];
    % plot port poke time
    itrigger = t_trigger_late(i);
    i_portin = t_portin-itrigger;
    i_portin = i_portin(i_portin>=-TriggerTimeDomain(1) & i_portin<=TriggerTimeDomain(2));
    if ~isempty(i_portin)
        i_portin = reshape(i_portin,1,[]);
        x_portin = [x_portin, i_portin];
        y_portin = [y_portin, (0.4-k)*ones(1,length(i_portin))];
    end
    k = k+1;
end

line(xx_all, yy_all, 'color', c_late, 'linewidth', 1);
line(xxrt_all, yyrt_all, 'color', c_cent_out, 'linewidth', 1.5);
scatter(x_portin, y_portin, 8, 'o', 'filled','MarkerFaceColor', c_reward,  'markerfacealpha', 0.5, 'MarkerEdgeColor','none')

line([0 0], get(gca, 'ylim'), 'color',c_trigger, 'linewidth', 1)
title('late', 'fontname', 'dejavu sans', 'fontsize', 7)
axis off

% Add information
uicontrol('Style','text','Units','centimeters','Position',[col4-0.5  yshift_row4 5 1.2],...
    'string', 'E. Trigger-related activity', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
    'HorizontalAlignment','Left');
yshift_row5=yshift_row4+1.2+1.2;

FRrange = [0 FRMax*1.1];
set(ha_cent_in_psth, 'ylim', FRrange);
line(ha_cent_in_psth, [0 0], FRrange, 'color', c_cent_in, 'linewidth', 1);

line(ha_cent_in_psth, [TargetFPs(1) TargetFPs(1)], get(gca, 'ylim'), 'color', c_trigger, 'linestyle', ':', 'linewidth', 1);
line(ha_cent_in_psth, [TargetFPs(2) TargetFPs(2)], get(gca, 'ylim'), 'color', c_trigger, 'linestyle', ':', 'linewidth', 1);

set(ha_cent_in_psth_error, 'ylim', FRrange);
line(ha_cent_in_psth_error, [0 0], FRrange, 'color', c_cent_in, 'linewidth', 1);
line(ha_cent_in_psth, [TargetFPs; TargetFPs], FRrange, 'color', c_trigger, 'linewidth', 1);
set(ha_cent_out_psth, 'ylim', FRrange);
line(ha_cent_out_psth, [0 0], FRrange, 'color', c_cent_out, 'linewidth', 1);
set(ha_cent_out_psth_error, 'ylim', FRrange);
line(ha_cent_out_psth_error, [0 0], FRrange, 'color', c_cent_out, 'linewidth', 1);
set(ha_poke, 'ylim', FRrange);
line(ha_poke, [0 0], FRrange, 'color', c_reward, 'linewidth', 1);
set(ha_trigger, 'ylim', FRrange);
line(ha_trigger, [0 0], FRrange, 'color', c_trigger, 'linewidth', 1);


%% plot spks
col5=19.5;
thiscolor = [0 0 0];
Lspk = size(r.Units.SpikeTimes(ku).wave, 2);
ha0=axes('unit', 'centimeters', 'position', [col5 yshift_row5 1.5 1.5], ...
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

ha00= axes('unit', 'centimeters', 'position', [19.5+1.5+1 yshift_row5 2 1.5], 'nextplot', 'add', 'xlim', [-25 25]);
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

yshift_row6 = yshift_row5+2;
% Plot all waveforms if it is a polytrode
if isfield(r.Units.SpikeTimes(ku), 'wave_mean')
    ha_wave_poly = axes('unit', 'centimeters', 'position', [col5 yshift_row6 4 3], 'nextplot', 'add');
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

    yshift_row7 = yshift_row6+3;
else
    yshift_row7 = yshift_row6;
end

uicontrol('Style','text','Units','centimeters','Position',[19 yshift_row7 5 1.5],...
    'string', 'F. Spike waveform and autocorrelation', ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 10,'BackgroundColor',[1 1 1],'ForegroundColor', 'k', ...
    'HorizontalAlignment','Left');
fig_height = max([fig_height, yshift_row7+2]);
% change the height of the figure
set(gcf, 'position', [2 2 25 fig_height])
toc;

if strcmpi(ToSave,'on')
    % save to a folder
    anm_name        =     r.BehaviorClass.Subject;
    session              =     r.BehaviorClass.Date;
    
    PSTH.ANM_Session = {anm_name, session};
    thisFolder = fullfile(pwd, 'Fig');
    if ~exist(thisFolder, 'dir')
        mkdir(thisFolder)
    end
    tosavename2= fullfile(thisFolder, [anm_name '_' session '_Ch'  num2str(ch) '_Unit' num2str(unit_no) ]);
    print (gcf,'-dpng', tosavename2)
    
    % save PSTH as well save(psth_new_name, 'PSTHOut');
    save([tosavename2 '.mat'], 'PSTH')
    
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