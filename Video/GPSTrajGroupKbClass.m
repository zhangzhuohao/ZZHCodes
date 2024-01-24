classdef GPSTrajGroupKbClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ANMInfoFile

        Sessions
        NumSessions
        ANM
        Treatment
        Dose
        Label
        Experimenter
        Task

        TrialInfo
        NumTrials
        TaskFP

        Trial
        Stage
        Performance
        PortCorrect
        PortChosen

        FP
        RT
        HD
        MT
        Cued

        TimeFromIn
        TimeFromOut
        TimeWarped

        PortVec
        PortCent

        AngleHead
        AngSpeedHead
        AngAccHead

        PosXHead
        PosYHead
        SpeedXHead
        SpeedYHead
        AccXHead
        AccYHead

        SpeedHead
        SpeedDirHead
        AccHead
        dPhiHead

        DistMatDtw = [];
        DistMatLw  = [];

%         AngleHeadTraceIn
%         AngleHeadTraceOut
%         AngleHeadTraceInTest
%         AngleHeadTraceOutTest
    end

    properties (Constant)
        Condition = ["Cue", "Uncue"];
        CueUncue    = [1 0];
        Ports = ["L", "R"];

        TimeMatIn  = -99:1:3000;
        TimeMatOut = -2099:1:1000;
        TimeMatWarp = -0.2+0.002:0.002:1.6;

%         ShuffleIters  = 200;
%         Alpha         = 0.01;
    end

    properties (Dependent)
        Ind

        TimeWarpAlign
        %% Trace matrix
%         AngleHeadMatIn
%         AngleHeadMatOut
%         AngleHeadMatWarp
%         AngSpeedHeadMatIn
%         AngSpeedHeadMatOut
%         AngSpeedHeadMatWarp
%         AngAccHeadMatIn
%         AngAccHeadMatOut
%         AngAccHeadMatWarp
% 
%         PosXHeadMatIn
%         PosXHeadMatOut
%         PosXHeadMatWarp
%         PosYHeadMatIn
%         PosYHeadMatOut
%         PosYHeadMatWarp
%         SpeedXHeadMatIn
%         SpeedXHeadMatOut
%         SpeedXHeadMatWarp
%         SpeedYHeadMatIn
%         SpeedYHeadMatOut
%         SpeedYHeadMatWarp
%         AccXHeadMatIn
%         AccXHeadMatOut
%         AccXHeadMatWarp
%         AccYHeadMatIn
%         AccYHeadMatOut
%         AccYHeadMatWarp
% 
%         SpeedHeadMatIn
%         SpeedHeadMatOut
%         SpeedHeadMatWarp
%         SpeedDirHeadMatIn
%         SpeedDirHeadMatOut
%         SpeedDirHeadMatWarp
%         AccHeadMatIn
%         AccHeadMatOut
%         AccHeadMatWarp
%         dPhiHeadMatIn
%         dPhiHeadMatOut
%         dPhiHeadMatWarp
    end

    methods
        function obj = GPSTrajGroupKbClass(TrajClassAll)
            %GPSTrajGroupKbClass Construct an instance of this class
            %% Session information
            obj.ANM             = unique(arrayfun(@(x) x.obj.ANM, TrajClassAll));
            obj.Sessions        = arrayfun(@(x) x.obj.Session, TrajClassAll);
            obj.Task            = unique(arrayfun(@(x) x.obj.Task, TrajClassAll));
            obj.NumSessions     = length(TrajClassAll);
            obj.Experimenter    = arrayfun(@(x) x.obj.Experimenter, TrajClassAll);

            obj.Treatment       = arrayfun(@(x) x.obj.Treatment, TrajClassAll);
            obj.Dose            = arrayfun(@(x) x.obj.Dose, TrajClassAll);
            obj.Label           = arrayfun(@(x) x.obj.Label, TrajClassAll);
            obj.NumTrials       = arrayfun(@(x) x.obj.NumTrials, TrajClassAll);
            obj.TaskFP          = unique(arrayfun(@(x) x.obj.SessionFP, TrajClassAll));

            %% Trial infotmation
            for i = 1:length(TrajClassAll)
                if i == 1
                    obj.TrialInfo = TrajClassAll(i).obj.TrialInfo;
                else
                    obj.TrialInfo = [obj.TrialInfo; TrajClassAll(i).obj.TrialInfo];
                end
            end

            allTrials           = arrayfun(@(x) x.obj.Trial, TrajClassAll, 'UniformOutput', false);
            obj.Trial           = cell2mat(allTrials');
            
            allStages           = arrayfun(@(x) x.obj.Stage, TrajClassAll, 'UniformOutput', false);
            obj.Stage           = cell2mat(allStages');

            allPerfs            = arrayfun(@(x) x.obj.Performance, TrajClassAll, 'UniformOutput', false);
            obj.Performance     = [];
            for i = 1:obj.NumSessions
                obj.Performance = [obj.Performance allPerfs{i}];
            end

            allPortCorrects     = arrayfun(@(x) x.obj.PortCorrect, TrajClassAll, 'UniformOutput', false);
            obj.PortCorrect     = [];
            for i = 1:obj.NumSessions
                obj.PortCorrect = [obj.PortCorrect allPortCorrects{i}];
            end
            
            allPortChosens      = arrayfun(@(x) x.obj.PortChosen, TrajClassAll, 'UniformOutput', false);
            obj.PortChosen      = [];
            for i = 1:obj.NumSessions
                obj.PortChosen  = [obj.PortChosen allPortChosens{i}];
            end
            
            

            allFPs              = arrayfun(@(x) x.obj.FP, TrajClassAll, 'UniformOutput', false);
            obj.FP              = cell2mat(allFPs');

            allHDs              = arrayfun(@(x) x.obj.HD, TrajClassAll, 'UniformOutput', false);
            obj.HD              = cell2mat(allHDs');

            allRTs              = arrayfun(@(x) x.obj.RT, TrajClassAll, 'UniformOutput', false);
            obj.RT              = cell2mat(allRTs');

            allMTs              = arrayfun(@(x) x.obj.MT, TrajClassAll, 'UniformOutput', false);
            obj.MT              = cell2mat(allMTs');

            allCued             = arrayfun(@(x) x.obj.Cued, TrajClassAll, 'UniformOutput', false);
            obj.Cued            = cell2mat(allCued');

            %% Time aligned
            allTimeIns          = arrayfun(@(x) x.obj.TimeFromIn, TrajClassAll, 'UniformOutput', false);
            obj.TimeFromIn      = [];
            for i = 1:obj.NumSessions
                obj.TimeFromIn  = [obj.TimeFromIn allTimeIns{i}];
            end

            allTimeOuts         = arrayfun(@(x) x.obj.TimeFromOut, TrajClassAll, 'UniformOutput', false);
            obj.TimeFromOut     = [];
            for i = 1:obj.NumSessions
                obj.TimeFromOut = [obj.TimeFromOut allTimeOuts{i}];
            end

            allTimeWarps        = arrayfun(@(x) x.obj.TimeWarped, TrajClassAll, 'UniformOutput', false);
            obj.TimeWarped      = [];
            for i = 1:obj.NumSessions
                obj.TimeWarped  = [obj.TimeWarped allTimeWarps{i}];
            end

            %% Port locations
            allPortVecs         = arrayfun(@(x) x.obj.PortVec, TrajClassAll, 'UniformOutput', false);
            obj.PortVec         = [];
            for i = 1:obj.NumSessions
                obj.PortVec     = [obj.PortVec allPortVecs{i}];
            end

            allPortCents        = arrayfun(@(x) x.obj.PortCent, TrajClassAll, 'UniformOutput', false);
            obj.PortCent        = [];
            for i = 1:obj.NumSessions
                obj.PortCent    = [obj.PortCent allPortCents{i}];
            end

            %% Body parts
            % angle of head
            allAngleHeads       = arrayfun(@(x) x.obj.AngleHead, TrajClassAll, 'UniformOutput', false);
            obj.AngleHead       = [];
            for i = 1:obj.NumSessions
                obj.AngleHead   = [obj.AngleHead allAngleHeads{i}];
            end

            allAngSpeedHeads    = arrayfun(@(x) x.obj.AngSpeedHead, TrajClassAll, 'UniformOutput', false);
            obj.AngSpeedHead    = [];
            for i = 1:obj.NumSessions
                obj.AngSpeedHead = [obj.AngSpeedHead allAngSpeedHeads{i}];
            end

            allAngAccHeads      = arrayfun(@(x) x.obj.AngAccHead, TrajClassAll, 'UniformOutput', false);
            obj.AngAccHead      = [];
            for i = 1:obj.NumSessions
                obj.AngAccHead  = [obj.AngAccHead allAngAccHeads{i}];
            end

            % position of head (X)
            allPosXHeads        = arrayfun(@(x) x.obj.PosXHead, TrajClassAll, 'UniformOutput', false);
            obj.PosXHead        = [];
            for i = 1:obj.NumSessions
                obj.PosXHead    = [obj.PosXHead allPosXHeads{i}];
            end

            allSpeedXHeads      = arrayfun(@(x) x.obj.SpeedXHead, TrajClassAll, 'UniformOutput', false);
            obj.SpeedXHead      = [];
            for i = 1:obj.NumSessions
                obj.SpeedXHead  = [obj.SpeedXHead allSpeedXHeads{i}];
            end

            allAccXHeads        = arrayfun(@(x) x.obj.AccXHead, TrajClassAll, 'UniformOutput', false);
            obj.AccXHead        = [];
            for i = 1:obj.NumSessions
                obj.AccXHead    = [obj.AccXHead allAccXHeads{i}];
            end

            % position of head (Y)
            allPosYHeads        = arrayfun(@(x) x.obj.PosYHead, TrajClassAll, 'UniformOutput', false);
            obj.PosYHead        = [];
            for i = 1:obj.NumSessions
                obj.PosYHead    = [obj.PosYHead allPosYHeads{i}];
            end

            allSpeedYHeads      = arrayfun(@(x) x.obj.SpeedYHead, TrajClassAll, 'UniformOutput', false);
            obj.SpeedYHead      = [];
            for i = 1:obj.NumSessions
                obj.SpeedYHead  = [obj.SpeedYHead allSpeedYHeads{i}];
            end

            allAccYHeads        = arrayfun(@(x) x.obj.AccYHead, TrajClassAll, 'UniformOutput', false);
            obj.AccYHead        = [];
            for i = 1:obj.NumSessions
                obj.AccYHead    = [obj.AccYHead allAccYHeads{i}];
            end

            % speed & acceleration
            allSpeedHeads        = arrayfun(@(x) x.obj.SpeedHead, TrajClassAll, 'UniformOutput', false);
            obj.SpeedHead        = [];
            for i = 1:obj.NumSessions
                obj.SpeedHead    = [obj.SpeedHead allSpeedHeads{i}];
            end

            allSpeedDirHeads      = arrayfun(@(x) x.obj.SpeedDirHead, TrajClassAll, 'UniformOutput', false);
            obj.SpeedDirHead      = [];
            for i = 1:obj.NumSessions
                obj.SpeedDirHead  = [obj.SpeedDirHead allSpeedDirHeads{i}];
            end

            allAccHeads        = arrayfun(@(x) x.obj.AccHead, TrajClassAll, 'UniformOutput', false);
            obj.AccHead        = [];
            for i = 1:obj.NumSessions
                obj.AccHead    = [obj.AccHead allAccHeads{i}];
            end

            alldPhiHeads        = arrayfun(@(x) x.obj.dPhiHead, TrajClassAll, 'UniformOutput', false);
            obj.dPhiHead        = [];
            for i = 1:obj.NumSessions
                obj.dPhiHead    = [obj.dPhiHead alldPhiHeads{i}];
            end

%             obj.AngleHeadTraceIn  = obj.getAngleHeadTrace("In");
%             obj.AngleHeadTraceOut = obj.getAngleHeadTrace("Out");

%             obj.AngleHeadTraceInTest  = obj.testTrace("In", obj.ShuffleIters, obj.Alpha);
%             obj.AngleHeadTraceOutTest = obj.testTrace("Out", obj.ShuffleIters, obj.Alpha);
        end

        %%
        function value = get.Ind(obj)

            stage = obj.Stage;

            ind.Cue     = obj.Cued==1 & stage;
            ind.Uncue   = obj.Cued==0 & stage;

            ind.ChooseL = obj.PortChosen=="L" & stage;
            ind.ChooseR = obj.PortChosen=="R" & stage;

            ind.TargetL = obj.PortCorrect=="L" & stage;
            ind.TargetR = obj.PortCorrect=="R" & stage;

            ind.Correct     = obj.Performance=="Correct"   & stage;
            ind.Wrong       = obj.Performance=="Wrong"     & stage;
            ind.Premature   = obj.Performance=="Premature" & stage;
            ind.Late        = obj.Performance=="Late"      & stage;

            value = ind;
        end

        function value = get.TimeWarpAlign(obj)
            
            num_max = max(cellfun(@(x) length(x), obj.TimeWarped));
            time_warp_align = linspace(0, 1, num_max);

            value = time_warp_align;
        end

        %% Aligna matrix
%         function value = get.AngleHeadMatIn(obj)
%             value = obj.alignMatrix(obj.AngleHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.AngleHeadMatOut(obj)
%             value = obj.alignMatrix(obj.AngleHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.AngleHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.AngleHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.AngSpeedHeadMatIn(obj)
%             value = obj.alignMatrix(obj.AngSpeedHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.AngSpeedHeadMatOut(obj)
%             value = obj.alignMatrix(obj.AngSpeedHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.AngSpeedHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.AngSpeedHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.AngAccHeadMatIn(obj)
%             value = obj.alignMatrix(obj.AngAccHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.AngAccHeadMatOut(obj)
%             value = obj.alignMatrix(obj.AngAccHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.AngAccHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.AngAccHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.PosXHeadMatIn(obj)
%             value = obj.alignMatrix(obj.PosXHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.PosXHeadMatOut(obj)
%             value = obj.alignMatrix(obj.PosXHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.PosXHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.PosXHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.PosYHeadMatIn(obj)
%             value = obj.alignMatrix(obj.PosYHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.PosYHeadMatOut(obj)
%             value = obj.alignMatrix(obj.PosYHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.PosYHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.PosYHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.SpeedXHeadMatIn(obj)
%             value = obj.alignMatrix(obj.SpeedXHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.SpeedXHeadMatOut(obj)
%             value = obj.alignMatrix(obj.SpeedXHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.SpeedXHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.SpeedXHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.SpeedYHeadMatIn(obj)
%             value = obj.alignMatrix(obj.SpeedYHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.SpeedYHeadMatOut(obj)
%             value = obj.alignMatrix(obj.SpeedYHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.SpeedYHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.SpeedYHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.AccXHeadMatIn(obj)
%             value = obj.alignMatrix(obj.AccXHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.AccXHeadMatOut(obj)
%             value = obj.alignMatrix(obj.AccXHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.AccXHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.AccXHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.AccYHeadMatIn(obj)
%             value = obj.alignMatrix(obj.AccYHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.AccYHeadMatOut(obj)
%             value = obj.alignMatrix(obj.AccYHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.AccYHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.AccYHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.SpeedHeadMatIn(obj)
%             value = obj.alignMatrix(obj.SpeedHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.SpeedHeadMatOut(obj)
%             value = obj.alignMatrix(obj.SpeedHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.SpeedHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.SpeedHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.SpeedDirHeadMatIn(obj)
%             value = obj.alignMatrix(obj.SpeedDirHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.SpeedDirHeadMatOut(obj)
%             value = obj.alignMatrix(obj.SpeedDirHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.SpeedDirHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.SpeedDirHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.AccHeadMatIn(obj)
%             value = obj.alignMatrix(obj.AccHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.AccHeadMatOut(obj)
%             value = obj.alignMatrix(obj.AccHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.AccHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.AccHead, obj.TimeWarped, obj.TimeMatWarp);
%         end
% 
%         function value = get.dPhiHeadMatIn(obj)
%             value = obj.alignMatrix(obj.dPhiHead, obj.TimeFromIn, obj.TimeMatIn);
%         end
%         function value = get.dPhiHeadMatOut(obj)
%             value = obj.alignMatrix(obj.dPhiHead, obj.TimeFromOut, obj.TimeMatOut);
%         end
%         function value = get.dPhiHeadMatWarp(obj)
%             value = obj.alignMatrix(obj.dPhiHead, obj.TimeWarped, obj.TimeMatWarp);
%         end

        %% Functions
        function [M, t_M] = trace2mat(~, trace, time_trace, time_matrix)
            
            num_timepoint = length(time_matrix);
            num_feature = size(trace{1}, 2);

            M = cellfun(@(trace, t) interp1(t, trace, time_matrix, "linear"), trace, time_trace, 'UniformOutput', false);
            M = M';
            M = cellfun(@(m) reshape(m, 1, num_timepoint, num_feature), M, 'UniformOutput', false);
            M = cell2mat(M);

            t_M = time_matrix;
        end

        function [T, t_T] = mat2trace(~, mat, time_matrix, time_trace)

            num_trace = size(mat, 1);
            num_timepoint = length(time_matrix);
            num_feature = size(mat, 3);

            if nargin<4
                time_trace = repmat(mat2cell(time_matrix', length(time_matrix)), 1, num_trace);
            end

            mat = mat2cell(mat, ones(num_trace, 1));
            mat = cellfun(@(m) reshape(m, num_timepoint, num_feature), mat, 'UniformOutput', false);

            T = cellfun(@(m, t) interp1(time_matrix, m, t, "linear"), mat, time_trace', 'UniformOutput', false);
            T = T';

            t_T = time_trace;
        end

        function A = trace2array(~, trace)

            A = cell2mat(trace');
        end

        function T = array2trace(~, array, time_trace)

            T = mat2cell(array, cellfun(@(x) length(x), time_trace));
            T = T';
        end

        %%
        function trace_normalized = normalizeTrace(~, trace, norm_method)
            
            if nargin<3
                norm_method = 'range';
            end
            trace_all = cell2mat(trace');
            trace_all_normalized = normalize(trace_all, 1, norm_method);

            trace_normalized = mat2cell(trace_all_normalized, cellfun(@(x) size(x, 1), trace));
            trace_normalized = trace_normalized';
        end

        function trace_gathered = gatherTrace(obj, body_parts)

            body_parts = string(body_parts);
            trace_gathered = obj.(body_parts(1));

            for i = 2:length(body_parts)
                trace_gathered = cellfun(@(x, y) [x y], trace_gathered, obj.(body_parts(i)), 'UniformOutput', false);
            end
        end

        %% Calculate distance matrix
        % dynamic time warp (DTW)
        function dist_mat = calDistDtw(obj, trace, trace_norm_method, dist_norm_method)

            if nargin > 2
                if ~strcmpi(trace_norm_method, 'none')
                    trace = obj.normalizeTrace(trace, trace_norm_method);
                end
            end

            num_trials = length(trace);
            dist_mat = zeros(num_trials);

            num_to_cal = .5*num_trials^2 - num_trials;
            num_calculated = 0;

            for m = 1:num_trials
                for n = m+1:num_trials
                    dist_mat(m, n) = dtw(trace{m}', trace{n}');
                    num_calculated = num_calculated + 1;

                    if ~mod(num_calculated, floor(num_to_cal/100))
                        fprintf('%.0f / %.0f    - %.0f%%\n', num_calculated, num_to_cal, 100*num_calculated/num_to_cal);
                    end
                end
            end
            dist_mat = dist_mat + dist_mat';

            if nargin > 3
                if ~strcmpi(dist_norm_method, 'none')
                    dist_mat = normalize(dist_mat(:), dist_norm_method);
                    dist_mat = reshape(dist_mat, num_trials, num_trials);
                end
            end
        end

        % linear warp (LW)
        function dist_mat = calDistLw(obj, trace, trace_norm_method, dist_norm_method)

            if nargin > 2
                if ~strcmpi(trace_norm_method, 'none')
                    trace = obj.normalizeTrace(trace, trace_norm_method);
                end
            end

            num_trials = length(trace);
            dist_mat = zeros(num_trials);

            num_to_cal = .5*num_trials^2 - num_trials;
            num_calculated = 0;

            for m = 1:num_trials
                for n = m+1:num_trials
                    dist_mat(m, n) = sum(sqrt(sum((trace{m}-trace{n}).^2, 2)));
                    num_calculated = num_calculated + 1;

                    if ~mod(num_calculated, floor(num_to_cal/100))
                        fprintf('%.0f / %.0f    - %.0f%%\n', num_calculated, num_to_cal, 100*num_calculated/num_to_cal);
                    end
                end
            end
            dist_mat = dist_mat + dist_mat';

            if nargin > 3
                if ~strcmpi(dist_norm_method, 'none')
                    dist_mat = normalize(dist_mat(:), dist_norm_method);
                    dist_mat = reshape(dist_mat, num_trials, num_trials);
                end
            end
        end

        %%
        function AngleHeadTrace = getAngleHeadTrace(obj, AlignTo)
            
            AlignTo = string(AlignTo);
            if ~ismember(AlignTo, ["In", "Out"])
                error("getAngleHeadTrace(obj, AlignTo): AlignTo should be one of: 'In', 'Out'.");
            end

            Perfs = ["Premature", "Correct", "Wrong", "Late"];

            for i = 1:length(obj.Condition)

                cue = obj.CueUncue(i);

                switch AlignTo
                    case {'In'}
                        time_bin_edges  = floor(obj.TimeMatIn(1)):20:1000*obj.TaskFP+300;
                        time_bin_center = time_bin_edges(1:end-1)+10;
                    case {'Out'}
                        time_bin_edges = -obj.TaskFP*1000-100:20:300;
                        time_bin_center = time_bin_edges(1:end-1)+10;
                end

                angle_head_trace.TimePoints = time_bin_center;

                for j = 1:length(obj.Ports)

                    for k = 1:length(Perfs)

                        Ind_this = find(obj.Stage==1 & obj.Cued==cue & obj.Ind.("Target"+obj.Ports(j)) & obj.Ind.(Perfs(k)));

                        if isempty(Ind_this)
                            continue
                        end

                        M = cell(1, length(time_bin_edges)-1);
                        for m = 1:length(time_bin_edges)-1

                            ang   = cellfun(@(x, y) x(y>=time_bin_edges(m) & y<time_bin_edges(m+1)), obj.AngleHead(Ind_this), obj.("TimeFrom"+AlignTo)(Ind_this), 'UniformOutput', false);
                            trial = cellfun(@(x, y) repmat(x, length(y), 1), num2cell(obj.Trial(Ind_this)), ang, 'UniformOutput', false);
                            time  = cellfun(@(y)    repmat(time_bin_center(m), length(y), 1), ang, 'UniformOutput', false);

                            port_correct = cellfun(@(x, y) repmat(x, length(y), 1), obj.PortCorrect(Ind_this), ang, 'UniformOutput', false);

                            ang   = cell2mat(ang');
                            trial = cell2mat(trial');
                            time  = cell2mat(time');

                            port_correct = cell2mat(port_correct');
                            port_correct = string(port_correct);
                            M{m} = table(ang, trial, time, port_correct);
                        end

                        angle_head_trace.(obj.Ports(j)+"_"+obj.Condition(i)+"_"+Perfs(k)) = M;

                    end

                    Ind_this = find(obj.Stage==1 & obj.Cued==cue & obj.Ind.("Target"+obj.Ports(j)));

                    if isempty(Ind_this)
                        continue
                    end

                    M = cell(1, length(time_bin_edges)-1);
                    for m = 1:length(time_bin_edges)-1

                        ang   = cellfun(@(x, y) x(y>=time_bin_edges(m) & y<time_bin_edges(m+1)), obj.AngleHead(Ind_this), obj.("TimeFrom"+AlignTo)(Ind_this), 'UniformOutput', false);
                        trial = cellfun(@(x, y) repmat(x, length(y), 1), num2cell(obj.Trial(Ind_this)), ang, 'UniformOutput', false);
                        time  = cellfun(@(y)    repmat(time_bin_center(m), length(y), 1), ang, 'UniformOutput', false);

                        port_correct = cellfun(@(x, y) repmat(x, length(y), 1), obj.PortCorrect(Ind_this), ang, 'UniformOutput', false);

                        ang   = cell2mat(ang');
                        trial = cell2mat(trial');
                        time  = cell2mat(time');

                        port_correct = cell2mat(port_correct');
                        port_correct = string(port_correct);

                        M{m} = table(ang, trial, time, port_correct);
                    end

                    angle_head_trace.(obj.Ports(j)+"_"+obj.Condition(i)) = M;
                end
            end

            AngleHeadTrace = angle_head_trace;
        end

        %% Save
        function save(obj, savepath)
            if nargin<2
                savepath = fullfile("D:\YuLab\Work\GPS\Video", obj.ANM);
            end
            save(fullfile(savepath, "GPSTrajGroupClass_"+obj.Task+"_"+upper(obj.ANM)+"_"+obj.Sessions(1)+"_"+obj.Sessions(end)), 'obj');
        end

        %% Plots
        function print(obj, Func, targetDir)
            
            if nargin==1
                Func = "HeatMap";
            end

            switch lower(Func)
                case {'heatmap'}
                    hf = obj.plotHeatMap("All");
                case {'trace'}
                    hf = obj.plotTrace();
            end

            if nargin==3
                % check if targetDir exists
                if ~contains(targetDir, '/') && ~contains(targetDir, '\')
                    % so it is a relative path
                    if ~exist(targetDir, 'dir')
                        mkdir(targetDir)
                    end
                end
            else
                targetDir = fullfile("D:\YuLab\Work\GPS\Video", obj.ANM);
            end
            savename = fullfile(targetDir, "GPSTrajectoryClass_" + obj.Task + "_" + string(Func) + "_" + upper(obj.ANM)+"_"+obj.Sessions(1)+"_"+obj.Sessions(end));
            print(hf, '-dpdf', savename, '-bestfit')
            print(hf, '-dpng', savename)
            saveas(hf, savename, 'fig')
        end

        function fig = plotHeatMap(obj, label)

            mycolormap = customcolormap_preset("red-white-blue");

            fig = figure(33); clf(33);
            set(fig, 'unit', 'centimeters', 'position', [2 2 19 25.5], 'paperpositionmode', 'auto', 'color', 'w');

            switch label
                case {"All"}
                    group_title = obj.ANM+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end);
                otherwise
                    group_title = obj.ANM+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end) + " ("+label+")";
            end
            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.05 0.9 0.9 0.08],...
                'string', group_title, 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            % Target L
            w_l = 1.5;
            ax_l = axes; colormap(mycolormap);
            set(ax_l, 'units', 'centimeters', 'position', [w_l 1.5, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHeadKb(ax_l, obj, "PortCorrect", "L", "Performance", "All", "AlignTo", "In", "Label", label);

            % Target R
            w_r = w_l + ax_l.Position(3) + .8;
            ax_r = axes; colormap(mycolormap);
            set(ax_r, 'units', 'centimeters', 'position', [w_r 1.5 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHeadKb(ax_r, obj, "PortCorrect", "R", "Performance", "All", "AlignTo", "In", "Label", label);
            set(ax_r, "ylabel", []);

            cb = colorbar(ax_r, "Units", "centimeters", "Position", [ax_r.Position(1)+ax_r.Position(3)+0.5 ax_l.Position(2) .3 ax_r.Position(4)+.2]);
            cb.Label.String = "Head angle (°)";
            cb.Label.FontSize = 9;
            cb.Label.FontWeight = "Bold";
            cb.Ticks = -90:30:90;

            cLimits = 90;
            clim(ax_l, [-1 1]*cLimits);
            clim(ax_r, [-1 1]*cLimits);

            h_fig = max([ax_l.Position(4) ax_r.Position(4)]) + 2.7;

            fig.Position(4) = h_fig;
            
        end
        
        function fig = plotTrace(obj)

            fig = figure(34); clf(34);
            set(fig, 'unit', 'centimeters', 'position', [2 2 8*obj.TaskFP+2 7.5], 'paperpositionmode', 'auto', 'color', 'w');

            group_title = obj.ANM+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end);
            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.05 0.9 0.9 0.08],...
                'string', group_title, 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            ax1 = axes();
            set(ax1, 'units', 'centimeters', 'position', [1.5 1.5 8*obj.TaskFP 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHeadKb(ax1, obj, "In");

        end

    end
end