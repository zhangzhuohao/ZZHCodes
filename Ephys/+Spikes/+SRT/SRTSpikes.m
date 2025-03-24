function PSTHOut = SRTSpikes(r, ind, varargin)
% V3: added a few new parameters
% V4: 2/17/2021 regular style drinking port. No IR sensor in front of the
% port.
% V5: add poke events following an unsuccesful release

% SRTSpikes(r, 13, 'FRrange', [0 35])

% ind can be singular or a vector

% 8.9.2020
% sort out spikes trains according to reaction time

% 1. same foreperiod, randked by reaction time
% 2. PSTH, two different foreperiods

% 5/7/2023 revised to adopt new FP schedule (e.g., 500 1000 1500)

set_matlab_default;
takeall = 0;

if isempty(ind)
    ind = 1:length(r.Units.SpikeTimes);
     takeall = 1;
else
    if length(ind) == 2
        ind_unit = find(r.Units.SpikeNotes(:, 1)==ind(1) & r.Units.SpikeNotes(:, 2)==ind(2));
        ind = ind_unit;
    end
end
ku_all = ind; % ind is used in different places

tic
ComputeRange = [];  % this is the range where time is extracted. Event times outside of this range will be discarded. Empty means taking everything

CentInTimeDomain  = [1000 2500]; % default
CentOutTimeDomain = [1500 1000];
ChoiceTimeDomain  = [1000 2000];
TriggerTimeDomain = [500  500];
InitInTimeDomain  = [1000 2000];
InitOutTimeDomain = [500  2000];

c_reward = [237 43 42] / 255;
ToSave = 'on';

rb = r.Behavior;
% all FPs
if length(r.BehaviorClass) > 1
    r.BehaviorClass = r.BehaviorClass{1};
end

TargetFPs = r.BehaviorClass.TargetFP; % you have to use BuildR2023 or BuildR4Tetrodes2023 to have this included in r.
TargetFPs = TargetFPs(TargetFPs>0);
nFPs      = length(TargetFPs);

Ports  = r.BehaviorClass.LeftRight;
nPorts = length(Ports);

%% 1 Align to different events
%% 1.1 Align to cent-in and cent-out

ind_cent_in  = find(strcmp(rb.Labels, 'PokeCentIn'));
t_cent_in    = shape_it(rb.EventTimings(rb.EventMarkers==ind_cent_in));

ind_cent_out = find(strcmp(rb.Labels, 'PokeCentOut'));
t_cent_out   = shape_it(rb.EventTimings(rb.EventMarkers==ind_cent_out));

disp(['Number of Cent poke-in is ' num2str(length(t_cent_in))]);
% 1.1.1 Correct trials

% index and time of correct cent_in
% get correct response
ind_correct        = find(strcmp(r.Behavior.Outcome, 'Correct'));
t_correct_cent_in  = t_cent_in(ind_correct);
t_correct_cent_out = t_cent_out(ind_correct);
FPs_correct        = shape_it(rb.Foreperiods(ind_correct));
Ports_correct      = shape_it(rb.PortCorrect(ind_correct));

RT_correct         = t_correct_cent_out - t_correct_cent_in - FPs_correct*1000;

% initialize sorted cells, sorted by FPs and Ports, ordered by reaction time
t_correct_cent_in_sort  = cell(nFPs, nPorts);
t_correct_cent_out_sort = cell(nFPs, nPorts);
RT_correct_sort         = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        % find this condition
        ind_ij = find(FPs_correct==TargetFPs(i) & Ports_correct==Ports(j));

        t_correct_cent_in_sort{i,j}  = t_correct_cent_in(ind_ij);
        t_correct_cent_out_sort{i,j} = t_correct_cent_out(ind_ij);
        RT_correct_sort{i,j}         = RT_correct(ind_ij);

        % sort order by hold duration
        [RT_correct_sort{i,j}, ind_sort] = sort(RT_correct_sort{i,j});

        t_correct_cent_in_sort{i,j}  = t_correct_cent_in_sort{i,j}(ind_sort);
        t_correct_cent_out_sort{i,j} = t_correct_cent_out_sort{i,j}(ind_sort);
    end
end
RT_cent_in_sort  = RT_correct_sort;
RT_cent_out_sort = RT_correct_sort;

% 1.1.2 Premature trials

%% Premature responses
ind_premature        = find(strcmp(r.Behavior.Outcome, 'Premature'));
t_premature_cent_in  = t_cent_in(ind_premature);
t_premature_cent_out = t_cent_out(ind_premature);
FPs_premature        = shape_it(rb.Foreperiods(ind_premature));
Ports_premature      = shape_it(rb.PortCorrect(ind_premature));

HD_premature         = t_premature_cent_out - t_premature_cent_in;

% initialize sorted cells, sorted by FPs and Ports, ordered by hold duration
HD_premature_sort    = cell(nFPs, nPorts);
FP_premature_sort    = cell(nFPs, nPorts);
Port_premature_sort  = cell(nFPs, nPorts);
t_premature_cent_in_sort  = cell(nFPs, nPorts);
t_premature_cent_out_sort = cell(nFPs, nPorts);

for i = 1:nFPs
    for j = 1:nPorts
        % find this condition
        ind_ij = find(FPs_premature==TargetFPs(i) & Ports_premature==Ports(j));

        t_premature_cent_in_sort{i,j}  = t_premature_cent_in(ind_ij);
        t_premature_cent_out_sort{i,j} = t_premature_cent_out(ind_ij);
        FP_premature_sort{i,j}         = FPs_premature(ind_ij);
        Port_premature_sort{i,j}       = Ports_premature(ind_ij);

        HD_premature_sort{i,j}         = HD_premature(ind_ij);

        % sort order by hold duration
        [HD_premature_sort{i,j}, ind_sort] = sort(HD_premature_sort{i,j});

        t_premature_cent_in_sort{i,j}  = t_premature_cent_in_sort{i,j}(ind_sort);
        t_premature_cent_out_sort{i,j} = t_premature_cent_out_sort{i,j}(ind_sort);
        FP_premature_sort{i,j}         = FP_premature_sort{i,j}(ind_sort);
        Port_premature_sort{i,j}       = Port_premature_sort{i,j}(ind_sort);
    end
end
% gather data in sorted order
t_premature_cent_in  = cell(1, nPorts);
t_premature_cent_out = cell(1, nPorts);
HD_premature_cent_in  = cell(1, nPorts);
HD_premature_cent_out = cell(1, nPorts);
FP_premature_cent_in  = cell(1, nPorts);
FP_premature_cent_out = cell(1, nPorts);

for j = 1:nPorts
    t_premature_cent_in{j}   = cell2mat(t_premature_cent_in_sort(:,j));
    t_premature_cent_out{j}  = cell2mat(t_premature_cent_out_sort(:,j));
    HD_premature_cent_in{j}  = cell2mat(HD_premature_sort(:,j));
    HD_premature_cent_out{j} = cell2mat(HD_premature_sort(:,j));
    FP_premature_cent_in{j}  = cell2mat(FP_premature_sort(:,j));
    FP_premature_cent_out{j} = cell2mat(FP_premature_sort(:,j));

    [~, ind_sort] = sort(HD_premature_cent_in{j});

    t_premature_cent_in{j}   = t_premature_cent_in{j}(ind_sort);
    t_premature_cent_out{j}  = t_premature_cent_out{j}(ind_sort);
    HD_premature_cent_in{j}  = HD_premature_cent_in{j}(ind_sort);
    HD_premature_cent_out{j} = HD_premature_cent_out{j}(ind_sort);
    FP_premature_cent_in{j}  = FP_premature_cent_in{j}(ind_sort);
    FP_premature_cent_out{j} = FP_premature_cent_out{j}(ind_sort);
end
% 1.1.3 Late trials

%% Late response
ind_late        = find(strcmp(r.Behavior.Outcome, 'Late'));
t_late_cent_in  = t_cent_in(ind_late);
t_late_cent_out = t_cent_out(ind_late);
FPs_late        = shape_it(rb.Foreperiods(ind_late));
Ports_late      = shape_it(rb.PortCorrect(ind_late));

HD_late         = t_late_cent_out - t_late_cent_in;

% initialize sorted cells, sorted by FPs and Ports
HD_late_sort         = cell(nFPs, nPorts);
FP_late_sort         = cell(nFPs, nPorts);
Port_late_sort       = cell(nFPs, nPorts);
t_late_cent_in_sort  = cell(nFPs, nPorts);
t_late_cent_out_sort = cell(nFPs, nPorts);

for i = 1:nFPs
    for j = 1:nPorts
        % find this condition
        ind_ij = find(FPs_late==TargetFPs(i) & Ports_late==Ports(j));

        t_late_cent_in_sort{i,j}  = t_late_cent_in(ind_ij);
        t_late_cent_out_sort{i,j} = t_late_cent_out(ind_ij);
        FP_late_sort{i,j}         = FPs_late(ind_ij);
        Port_late_sort{i,j}       = Ports_late(ind_ij);

        HD_late_sort{i,j}         = HD_late(ind_ij);

        % sort by hold duration
        [HD_late_sort{i,j}, ind_sort] = sort(HD_late_sort{i,j});

        t_late_cent_in_sort{i,j}  = t_late_cent_in_sort{i,j}(ind_sort);
        t_late_cent_out_sort{i,j} = t_late_cent_out_sort{i,j}(ind_sort);
        FP_late_sort{i,j}         = FP_late_sort{i,j}(ind_sort);
        Port_late_sort{i,j}       = Port_late_sort{i,j}(ind_sort);
    end
end
% gather data in sorted order
t_late_cent_in  = cell(1, nPorts);
t_late_cent_out = cell(1, nPorts);
HD_late_cent_in  = cell(1, nPorts);
HD_late_cent_out = cell(1, nPorts);
FP_late_cent_in  = cell(1, nPorts);
FP_late_cent_out = cell(1, nPorts);

for j = 1:nPorts
    t_late_cent_in{j}  = cell2mat(t_late_cent_in_sort(:,j));
    t_late_cent_out{j} = cell2mat(t_late_cent_out_sort(:,j));
    HD_late_cent_in{j}  = cell2mat(HD_late_sort(:,j));
    HD_late_cent_out{j} = cell2mat(HD_late_sort(:,j));
    FP_late_cent_in{j}  = cell2mat(FP_late_sort(:,j));
    FP_late_cent_out{j} = cell2mat(FP_late_sort(:,j));
end

% 1.1.4 Probe trials

%% Probe trials
ind_probe        = find(strcmp(rb.Outcome, 'Probe'));
t_probe_cent_in  = t_cent_in(ind_probe);
t_probe_cent_out = t_cent_out(ind_probe);
Ports_probe      = shape_it(rb.PortCorrect(ind_probe));

HD_probe         = t_probe_cent_out - t_probe_cent_in;

% initialize sorted cells, sorted by FPs and Ports
HD_probe_sort         = cell(1, nPorts);
Port_probe_sort       = cell(1, nPorts);
t_probe_cent_in_sort  = cell(1, nPorts);
t_probe_cent_out_sort = cell(1, nPorts);

for j = 1:nPorts
    % find this condition
    ind_j = find(Ports_probe==Ports(j));

    t_probe_cent_in_sort{j}  = t_probe_cent_in(ind_j);
    t_probe_cent_out_sort{j} = t_probe_cent_out(ind_j);
    Port_probe_sort{j}       = Ports_probe(ind_j);

    HD_probe_sort{j}         = HD_probe(ind_j);

    % sort by hold duration
    [HD_probe_sort{j}, ind_sort] = sort(HD_probe_sort{j});

    t_probe_cent_in_sort{j}  = t_probe_cent_in_sort{j}(ind_sort);
    t_probe_cent_out_sort{j} = t_probe_cent_out_sort{j}(ind_sort);
    Port_probe_sort{j}       = Port_probe_sort{j}(ind_sort);
end

%% 1.2 Align to choice-in

%% Choice
tmin          = 200; % allow at least 0.2 second between a successful cent_out and poke
tmax          = 2000; % allow at most 2 second between a successful cent_out and poke

ind_choice    = find(strcmp(rb.Labels, 'PokeChoiceIn'));
t_choice      = shape_it(rb.EventTimings(rb.EventMarkers==ind_choice));
disp(['Number of Choice poke-in is ' num2str(length(t_choice))]);
% 1.2.1 Rewarded choice poke (correct trials)

FP_correct_choice   = zeros(1, length(t_choice)); % find out choice_in FP associated with each reward
Port_correct_choice = zeros(1, length(t_choice)); % find out choice_in Port associated with each reward

mt_correct          = zeros(1, length(t_choice));
for i = 1:length(t_choice)
    dt = t_choice(i) - t_correct_cent_out;
    dt = dt(dt>tmin & dt<tmax); % reward must be collected within 2 sec after a correct cent_out
    if ~isempty(dt)
        mt_correct(i) = dt(end);
        ind = find(t_correct_cent_out==t_choice(i)-dt(end));
        if ~isempty(ind)
            FP_correct_choice(i)   = FPs_correct(ind);
            Port_correct_choice(i) = Ports_correct(ind);
        else
            disp('Not found')
        end
    else
        mt_correct(i) = NaN;
    end
end

ind_valid = ~isnan(mt_correct);
mt_correct          = mt_correct(ind_valid);
t_correct_choice    = t_choice(ind_valid);
FP_correct_choice   = FP_correct_choice(ind_valid);
Port_correct_choice = Port_correct_choice(ind_valid);
% Check movement time distribution
Edges = (0:50:2000);
figure(45); clf; set(gcf, 'Units', 'Centimeters', 'Position', [5 5 6 4]);
histogram(mt_correct, Edges);
xlabel('Movement time (ms)');
ylabel('Count');

% sort choice according to FP and Port
t_correct_choice_sorted = cell(nFPs, nPorts);
mt_correct_sorted       = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        ind_ij = find(FP_correct_choice==TargetFPs(i) & Port_correct_choice==Ports(j));
        t_correct_choice_sorted{i,j} = t_correct_choice(ind_ij);
        mt_correct_sorted{i,j}       = mt_correct(ind_ij);

        % rank them
        [mt_correct_sorted{i,j}, ind_sort] = sort(mt_correct_sorted{i,j});
        t_correct_choice_sorted{i,j}       = t_correct_choice_sorted{i,j}(ind_sort);
    end
end
mt_correct = mt_correct_sorted;
t_correct_choice = t_correct_choice_sorted;

% 1.2.2 Non-rewarded choice poke (wrong or probe trials)

ind_noreward         = find(~strcmp(rb.Outcome, 'Correct'));
FP_noreward_choice   = zeros(1, length(t_choice)); % find out choice_in FP associated with each reward
Port_noreward_choice = zeros(1, length(t_choice)); % find out choice_in Port associated with each reward

t_noreward_cent_out = t_cent_out(ind_noreward);
FPs_noreward        = shape_it(rb.Foreperiods(ind_noreward));
PortChosen_noreward = shape_it(rb.PortChosen(ind_noreward));

mt_noreward         = zeros(1, length(t_choice));
for i = 1:length(t_choice)
    dt = t_choice(i) - t_noreward_cent_out;
    dt = dt(dt>tmin & dt<tmax); % reward must be collected within 2 sec after a noreward cent_out
    if ~isempty(dt)
        mt_noreward(i) = dt(end);
        ind = find(t_noreward_cent_out==t_choice(i)-dt(end));
        if ~isempty(ind)
            FP_noreward_choice(i)   = FPs_noreward(ind);
            Port_noreward_choice(i) = PortChosen_noreward(ind);
        else
            disp('Not found')
        end
    else
        mt_noreward(i) = NaN;
    end
end

ind_valid = ~isnan(mt_noreward);
mt_noreward          = mt_noreward(ind_valid);
t_noreward_choice    = t_choice(ind_valid);
FP_noreward_choice   = FP_noreward_choice(ind_valid);
Port_noreward_choice = Port_noreward_choice(ind_valid);
% Check movement time distribution
Edges = (0:50:2000);
figure(45); clf; set(gcf, 'Units', 'Centimeters', 'Position', [5 5 6 4]);
histogram(mt_noreward, Edges);
xlabel('Movement time (ms)');
ylabel('Count');

% sort by movement time
% rank them
% sort choice according to FP and Port
t_noreward_choice_sorted = cell(1, nPorts);
mt_noreward_sorted       = cell(1, nPorts);
for j = 1:nPorts
    ind_j = find(Port_noreward_choice==Ports(j));
    t_noreward_choice_sorted{j} = t_noreward_choice(ind_j);
    mt_noreward_sorted{j}       = mt_noreward(ind_j);

    % rank them
    [mt_noreward_sorted{j}, ind_sort] = sort(mt_noreward_sorted{j});
    t_noreward_choice_sorted{j}       = t_noreward_choice_sorted{j}(ind_sort);
end
mt_noreward = mt_noreward_sorted;
t_noreward_choice = t_noreward_choice_sorted;

%% 1.3 Align to trigger stimulus

%% Trigger
ind_trig = find(strcmp(rb.Labels, 'Trigger'));
t_trig   = shape_it(rb.EventTimings(rb.EventMarkers==ind_trig)); % trigger time in ms.

Types_trig = cell(1, length(t_trig));
FPs_trig   = zeros(1, length(t_trig));
Ports_trig = zeros(1, length(t_trig));
RTs_trig   = nan(1, length(t_trig)); % reaction time (used for ranking later)

for i = 1:length(t_trig)
    % find the most recent cent_in
    ind_recent_cent_in = find(t_cent_in<t_trig(i), 1, 'last');
    if ~isempty(ind_recent_cent_in) && abs(t_cent_in(ind_recent_cent_in)-t_trig(i))<2500
        % check the condition
        Types_trig{i} = rb.Outcome{ind_recent_cent_in};
        FPs_trig(i)   = rb.Foreperiods(ind_recent_cent_in);
        Ports_trig(i) = rb.PortCorrect(ind_recent_cent_in);
        ind_following_cent_out = find(t_cent_out>t_trig(i), 1, 'first');
        if ~isempty(ind_following_cent_out)
            RTs_trig(i) = t_cent_out(ind_following_cent_out) - t_trig(i);
        end
    else
        Types_trig{i} = 'NaN';
        FPs_trig(i)   = NaN;
        Ports_trig(i) = NaN;
    end
end

t_late_trig  = cell(1, nPorts);
RT_late_trig = cell(1, nPorts);
FP_late_trig = cell(1, nPorts);
for j = 1:nPorts
    % sort trigger (to plot)
    t_late_trig{j}  = t_trig(strcmp(Types_trig, 'Late') & Ports_trig==Ports(j));
    RT_late_trig{j} = RTs_trig(strcmp(Types_trig, 'Late') & Ports_trig==Ports(j));
    FP_late_trig{j} = FPs_trig(strcmp(Types_trig, 'Late') & Ports_trig==Ports(j));
    % sort according to reaction time
    [RT_late_trig{j}, ind_sort] = sort(RT_late_trig{j});
    t_late_trig{j}  = t_late_trig{j}(ind_sort);
    FP_late_trig{j}  = FP_late_trig{j}(ind_sort);
end

% trigger according to FPs
t_trig_sort  = cell(nFPs, nPorts);
RT_trig_sort = cell(nFPs, nPorts);

for i = 1:nFPs
    for j = 1:nPorts
        % sort trigger (to plot)
        t_trig_sort{i,j}  = t_trig(strcmp(Types_trig, 'Correct') & FPs_trig==TargetFPs(i) & Ports_trig==Ports(j));
        RT_trig_sort{i,j} = RTs_trig(strcmp(Types_trig, 'Correct') & FPs_trig==TargetFPs(i) & Ports_trig==Ports(j));
        % sort according to reaction time
        [RT_trig_sort{i,j}, ind_sort] = sort(RT_trig_sort{i,j});
        t_trig_sort{i,j}  = t_trig_sort{i,j}(ind_sort);
    end
end

%% 1.4 Align to init-in and init-out

ind_init_in  = find(strcmp(rb.Labels, 'PokeInitIn'));
t_init_in    = shape_it(rb.EventTimings(rb.EventMarkers==ind_init_in));

ind_init_out = find(strcmp(rb.Labels, 'PokeInitOut'));
t_init_out   = shape_it(rb.EventTimings(rb.EventMarkers==ind_init_out));

if length(t_init_in)>length(t_init_out)
    t_init_in(end) = [];
end
disp(['Number of Init poke-in is ' num2str(length(t_init_in))]);

InitDur = t_init_out - t_init_in;
[InitDur_sort, ind_sort] = sort(InitDur);
t_init_in_sort           = t_init_in(ind_sort);
t_init_out_sort          = t_init_out(ind_sort);

%% 2 Check compute range

%% Check ComputeRange
if ~isempty(ComputeRange)
    t_cent_in(t_cent_in<ComputeRange(1) | t_cent_in>ComputeRange(2))=[];
    
    for i=1:nFPs
        to_remove = find(t_correct_cent_in_sort{i}<ComputeRange(1) | t_correct_cent_in_sort{i}>ComputeRange(2) | t_correct_cent_out_sort{i}<ComputeRange(1) | t_correct_cent_out_sort{i}>ComputeRange(2));
        t_correct_cent_out_sort{i}(to_remove) = [];
        RT_cent_out_sort{i}(to_remove)        = [];
        t_correct_cent_in_sort{i}(to_remove)  = [];
        RT_cent_in_sort{i}(to_remove)         = [];
    end
    to_remove_premature = find(t_premature_cent_in<ComputeRange(1) | t_premature_cent_in>ComputeRange(2) | t_premature_cent_out<ComputeRange(1) | t_premature_cent_out>ComputeRange(2));
    t_premature_cent_in(to_remove_premature)         = [];
    HD_premature_cent_in(to_remove_premature)  = [];
    FP_premature_cent_in(to_remove_premature)       = [];
    t_premature_cent_out(to_remove_premature)        = [];
    HD_premature_cent_out(to_remove_premature) = [];
    FP_premature_cent_out(to_remove_premature)      = [];
    
    to_remove_late = find(t_late_cent_in<ComputeRange(1) | t_late_cent_in>ComputeRange(2) | t_late_cent_out<ComputeRange(1) | t_late_cent_out>ComputeRange(2));
    t_late_cent_in(to_remove_late)         = [];
    late_duration_cent_in(to_remove_late)  = [];
    FPs_late_cent_in(to_remove_late)       = [];
    t_late_cent_out(to_remove_late)        = [];
    late_duration_cent_out(to_remove_late) = [];
    FPs_late_cent_out(to_remove_late)      = [];
    
    for i=1:nFPs
        to_remove_rewardpokes =  find(t_reward_pokes{i}<ComputeRange(1) | t_reward_pokes{i}>ComputeRange(2));
        t_reward_pokes{i}(to_remove_rewardpokes) = [];
        mt_correct{i}(to_remove_rewardpokes) = [];
    end
    
    to_remove_nonrewardpokes = find(t_nonreward_pokes<ComputeRange(1) | t_nonreward_pokes>ComputeRange(2));
    t_nonreward_pokes(to_remove_nonrewardpokes) = [];
    move_time_nonreward(to_remove_nonrewardpokes) = [];
    
    for i = 1:nFPs
        to_remove = find(t_trig_sort{i}<ComputeRange(1) | t_trig_sort{i}>ComputeRange(2));
        t_trig_sort{i}(to_remove) = [];
        RT_trig_sort{i}(to_remove) = [];
    end
    
    to_remove_trigger_late = find(t_late_trig<ComputeRange(1) | t_late_trig>ComputeRange(2));
    t_late_trig(to_remove_trigger_late) = [];
    RT_late_trig(to_remove_trigger_late) = [];
    FP_late_trig(to_remove_trigger_late) = [];
end

%% 3 Summerize

%% Summarize
PSTHOut.ANM_Session = {r.BehaviorClass.Subject, r.BehaviorClass.Session};
PSTHOut.SpikeNotes  = r.Units.SpikeNotes;

%% 3.1 Cent-In

n_sort = numel(t_correct_cent_in_sort);
PSTHOut.CentIn.Labels = [repmat({'Correct'}, 1, n_sort), 'Premature', 'Premature', 'Late', 'Late', 'Probe', 'Probe', 'All'];

PSTHOut.CentIn.Time = cell(1, length(PSTHOut.CentIn.Labels));
PSTHOut.CentIn.Time(1:n_sort) = t_correct_cent_in_sort;
PSTHOut.CentIn.Time(n_sort+(1:2)) = t_premature_cent_in;
PSTHOut.CentIn.Time(n_sort+(3:4)) = t_late_cent_in;
PSTHOut.CentIn.Time(n_sort+(5:6)) = t_probe_cent_in_sort;
PSTHOut.CentIn.Time(n_sort+7) = {t_cent_in};

PSTHOut.CentIn.FP                = {repmat(TargetFPs, 1, 2), FP_premature_cent_in, FP_late_cent_in, nan, rb.Foreperiods};
PSTHOut.CentIn.Port              = {[ones(1, nFPs) 2*ones(1, nFPs)], Ports, Ports, Ports, rb.PortCorrect};
PSTHOut.CentIn.RT_Correct        = RT_cent_in_sort;
PSTHOut.CentIn.HoldDur.Premature = HD_premature_cent_in;
PSTHOut.CentIn.HoldDur.Late      = HD_late_cent_in;
PSTHOut.CentIn.HoldDur.Probe     = HD_probe_sort;

%% 3.2 Cent-Out

n_sort = numel(t_correct_cent_out_sort);
PSTHOut.CentOut.Labels = [repmat({'Correct'}, 1, n_sort), 'Premature', 'Premature', 'Late', 'Late', 'Probe', 'Probe'];

PSTHOut.CentOut.Time   = cell(1, length(PSTHOut.CentOut.Labels));
PSTHOut.CentOut.Time(1:n_sort)  = t_correct_cent_out_sort;
PSTHOut.CentOut.Time(n_sort+(1:2))  = t_premature_cent_out;
PSTHOut.CentOut.Time(n_sort+(3:4))  = t_late_cent_out;
PSTHOut.CentOut.Time(n_sort+(5:6))  = t_probe_cent_out_sort;

PSTHOut.CentOut.FP                = {repmat(TargetFPs, 1, 2), FP_premature_cent_out, FP_late_cent_out, nan};
PSTHOut.CentOut.Port              = {[ones(1, nFPs) 2*ones(1, nFPs)], Ports, Ports, Ports};
PSTHOut.CentOut.RT_Correct        = RT_cent_out_sort;
PSTHOut.CentOut.HoldDur.Premature = HD_premature_cent_out;
PSTHOut.CentOut.HoldDur.Late      = HD_late_cent_out;
PSTHOut.CentOut.HoldDur.Probe     = HD_probe_sort;

%% 3.3 Choice-In

PSTHOut.ChoiceIn.Time                    = t_choice;
PSTHOut.ChoiceIn.RewardPoke.Time         = t_correct_choice; % it is a cell now!
PSTHOut.ChoiceIn.RewardPoke.Move_Time    = mt_correct;       % it is a cell now!
PSTHOut.ChoiceIn.NonrewardPoke.Time      = t_noreward_choice;
PSTHOut.ChoiceIn.NonrewardPoke.Move_Time = mt_noreward;

%% 3.4 Trigger

n_sort = numel(t_trig_sort);
PSTHOut.Triggers.Labels          = {'TriggerTime_DifferentFPs' 'TriggerTime_Late'};
PSTHOut.Triggers.Time            = cell(1, n_sort+2);
PSTHOut.Triggers.Time(1:n_sort)  = t_trig_sort;
PSTHOut.Triggers.Time(end-1:end) = t_late_trig;

PSTHOut.Triggers.RT             = [RT_trig_sort(:)', RT_late_trig];
PSTHOut.Triggers.FP             = {repmat(TargetFPs, 1, 2), FP_late_trig};
PSTHOut.Triggers.Port           = {[ones(1, nFPs) 2*ones(1, nFPs)], Ports};

%% 3.5 InitIn & InitOut

PSTHOut.InitIn.Time     = t_init_in_sort;
PSTHOut.InitIn.InitDur  = InitDur_sort;

PSTHOut.InitOut.Time    = t_init_out_sort;
PSTHOut.InitOut.InitDur = InitDur_sort;

%% 4 Compute PSTH

%% Check how many units we need to compute
% derive PSTH from these
% go through each units if necessary
for iku = 1:length(ku_all)
    ku = ku_all(iku);
    if ku>length(r.Units.SpikeTimes)
        disp('##########################################')
        disp('########### That is all you have ##############')
        disp('##########################################')
        return
    end
    disp('##########################################')
    disp(['Computing this unit: ' num2str(ku)])
    disp('##########################################')
    PSTHOut.PSTH(iku) = Spikes.SRT.ComputePlotPSTH(r, PSTHOut, ku,...
        'CentInTimeDomain', CentInTimeDomain, ...
        'CentOutTimeDomain', CentOutTimeDomain, ...
        'ChoiceTimeDomain', ChoiceTimeDomain,...
        'TriggerTimeDomain', TriggerTimeDomain,...
        'InitInTimeDomain', InitInTimeDomain,...
        'InitOutTimeDomain', InitOutTimeDomain,...
        'ToSave', ToSave);
end

%%

if takeall
    
    r.PSTH.Events.CentIn     = PSTHOut.CentIn;
    r.PSTH.Events.CentOut    = PSTHOut.CentOut;
    r.PSTH.Events.Pokes      = PSTHOut.ChoiceIn;
    r.PSTH.Events.Triggers   = PSTHOut.Triggers;
    r.PSTH.Events.InitIn     = PSTHOut.InitIn;
    r.PSTH.Events.InitOut    = PSTHOut.InitOut;
    r.PSTH.PSTHs             = PSTHOut.PSTH;
    
    r_name = 'RTarray_'+r.BehaviorClass.Subject+'_'+r.BehaviorClass.Session+'.mat';
    save(r_name, 'r');
    psth_new_name = r.BehaviorClass.Subject+'_'+r.BehaviorClass.Session+'_PSTHs.mat';
    save(psth_new_name, 'PSTHOut');
    try
        % C:\Users\jiani\OneDrive\00_Work\03_Projects\05_Physiology\PSTHs
        thisFolder = fullfile(findonedrive, '00_Work', '03_Projects', '05_Physiology', 'Data', 'PETHs', r.BehaviorClass.Subject);
        if ~exist(thisFolder, 'dir')
            mkdir(thisFolder);
        end
        copyfile(psth_new_name, thisFolder);
    end
end

end