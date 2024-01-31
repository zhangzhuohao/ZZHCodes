classdef GPSTrajGroupClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
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

        Trial
        Stage
        Performance
        PortCorrect
        PortChosen

        FP
        RT
        HD
        MT

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
% 
%         AngleHeadTraceIn
%         AngleHeadTraceOut
%         AngleHeadTraceInTest
%         AngleHeadTraceOutTest
    end

    properties (Constant)
        ForePeriods = ["Short", "Med", "Long"];
        MixedFPs    = [.5 1 1.5];
        Ports = ["L", "R"];

        TimeMatIn  = -100:10:2500;
        TimeMatOut = -1600:10:1000;
        TimeMatWarp = -0.2:0.001:1.6;

%         ShuffleIters  = 200;
%         Alpha         = 0.01;
    end

    properties (Dependent)
        Ind

        TimeWarpAlign
    end

    methods
        function obj = GPSTrajGroupClass(TrajClassAll)
            %GPSTrajGroupClass Construct an instance of this class
            
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
            
            allPortChosens     = arrayfun(@(x) x.obj.PortChosen, TrajClassAll, 'UniformOutput', false);
            obj.PortChosen     = [];
            for i = 1:obj.NumSessions
                obj.PortChosen = [obj.PortChosen allPortChosens{i}];
            end

            allFPs             = arrayfun(@(x) x.obj.FP, TrajClassAll, 'UniformOutput', false);
            obj.FP             = cell2mat(allFPs');

            allHDs             = arrayfun(@(x) x.obj.HD, TrajClassAll, 'UniformOutput', false);
            obj.HD             = cell2mat(allHDs');

            allRTs             = arrayfun(@(x) x.obj.RT, TrajClassAll, 'UniformOutput', false);
            obj.RT             = cell2mat(allRTs');

            allMTs             = arrayfun(@(x) x.obj.MT, TrajClassAll, 'UniformOutput', false);
            obj.MT             = cell2mat(allMTs');

            %% Time aligned
            allTimeIns         = arrayfun(@(x) x.obj.TimeFromIn, TrajClassAll, 'UniformOutput', false);
            obj.TimeFromIn     = [];
            for i = 1:obj.NumSessions
                obj.TimeFromIn = [obj.TimeFromIn allTimeIns{i}];
            end

            allTimeOuts        = arrayfun(@(x) x.obj.TimeFromOut, TrajClassAll, 'UniformOutput', false);
            obj.TimeFromOut    = [];
            for i = 1:obj.NumSessions
                obj.TimeFromOut = [obj.TimeFromOut allTimeOuts{i}];
            end

            allTimeWarps        = arrayfun(@(x) x.obj.TimeWarped, TrajClassAll, 'UniformOutput', false);
            obj.TimeWarped      = [];
            for i = 1:obj.NumSessions
                obj.TimeWarped  = [obj.TimeWarped allTimeWarps{i}];
            end

            %% Port locations
            allPortVecs        = arrayfun(@(x) x.obj.PortVec, TrajClassAll, 'UniformOutput', false);
            obj.PortVec        = [];
            for i = 1:obj.NumSessions
                obj.PortVec = [obj.PortVec allPortVecs{i}];
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

            ind.Short = obj.FP==0.5 & stage;
            ind.Med   = obj.FP==1.0 & stage;
            ind.Long  = obj.FP==1.5 & stage;

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
        function dist_mat = calDist(obj, trace, varargin)

            % parsing input
            P = inputParser;

            addParameter(P, 'warp_method', 'linear', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'trace_norm_method', 'none', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'dist_norm_method', 'none', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'disp', 0, @isinteger);

            parse(P, varargin{:});

            warp_method = P.Results.warp_method;
            trace_norm_method = P.Results.trace_norm_method;
            dist_norm_method = P.Results.dist_norm_method;
            disp = P.Results.disp;

            switch trace_norm_method
                case 'none'
                otherwise
                    trace = obj.normalizeTrace(trace, trace_norm_method);
            end

            num_trials = length(trace);
            dist_mat = zeros(num_trials);

            num_to_cal = .5*num_trials^2 - num_trials;
            num_calculated = 0;

            for m = 1:num_trials
                for n = m+1:num_trials
                    switch warp_method
                        case 'dtw'
                            dist_mat(m, n) = dtw(trace{m}', trace{n}');
                        case 'linear'
                            dist_mat(m, n) = sum(sqrt(sum((trace{m}-trace{n}).^2, 2)));
                    end
                    num_calculated = num_calculated + 1;

                    if disp
                        if ~mod(num_calculated, floor(num_to_cal/100))
                            fprintf('%.0f / %.0f    - %.0f%%\n', num_calculated, num_to_cal, 100*num_calculated/num_to_cal);
                        end
                    end
                end
            end
            dist_mat = dist_mat + dist_mat';

            switch dist_norm_method
                case 'none'
                otherwise
                    dist_mat = normalize(dist_mat(:), dist_norm_method);
                    dist_mat = reshape(dist_mat, num_trials, num_trials);
            end
        end

        %%
        function AngleHeadTrace = getAngleHeadTrace(obj, AlignTo)
            
            AlignTo = string(AlignTo);
            if ~ismember(AlignTo, ["In", "Out"])
                error("getAngleHeadTrace(obj, AlignTo): AlignTo should be one of: 'In', 'Out'.");
            end

            Perfs = ["Correct", "Wrong"];

            for i = 1:length(obj.ForePeriods)

                fp = obj.MixedFPs(i);

                switch AlignTo
                    case {'In'}
                        time_bin_edges  = floor(obj.TimeMatIn(1)):20:1000*fp+300;
                        time_bin_center = time_bin_edges(1:end-1)+10;
                    case {'Out'}
                        time_bin_edges = -fp*1000-100:20:200;
                        time_bin_center = time_bin_edges(1:end-1)+10;
                end

                angle_head_trace.("TimePoints_"+obj.ForePeriods(i)) = time_bin_center;

                for j = 1:length(obj.Ports)

                    for k = 1:length(Perfs)

                        Ind_this = find(ismember(obj.Performance, Perfs(k)) & obj.FP==fp & obj.Ind.("Choose"+obj.Ports(j)));
                        
                        M = cell(1, length(time_bin_edges)-1);
                        for m = 1:length(time_bin_edges)-1

                            ang   = cellfun(@(x, y) x(y>=time_bin_edges(m) & y<time_bin_edges(m+1)), obj.AngleHead(Ind_this), obj.("TimeFrom"+AlignTo)(Ind_this), 'UniformOutput', false);
                            trial = cellfun(@(x, y) repmat(x, length(y), 1), num2cell(obj.Trial(Ind_this)), ang, 'UniformOutput', false);
                            time  = cellfun(@(y)    repmat(time_bin_center(m), length(y), 1), ang, 'UniformOutput', false);

                            port_chosen  = cellfun(@(x, y) repmat(x, length(y), 1), obj.PortChosen(Ind_this), ang, 'UniformOutput', false);
                            port_correct = cellfun(@(x, y) repmat(x, length(y), 1), obj.PortCorrect(Ind_this), ang, 'UniformOutput', false);

                            ang   = cell2mat(ang');
                            trial = cell2mat(trial');
                            time  = cell2mat(time');

                            port_chosen = cell2mat(port_chosen');
                            port_chosen = string(port_chosen);
                            port_correct = cell2mat(port_correct');
                            port_correct = string(port_correct);

                            M{m} = table(ang, trial, time, port_chosen, port_correct);
                        end

                        angle_head_trace.(obj.Ports(j)+"_"+obj.ForePeriods(i)+"_"+Perfs(k)) = M;
                    end
                end
            end

            AngleHeadTrace = angle_head_trace;
        end

        function Out = testTrace(obj, align_to, shuffle_iters, alpha)
            
            switch align_to
                case {'In', 'in', 'IN'}
                    align_to = "In";
                case {'Out', 'out', 'OUT'}
                    align_to = "Out";
                otherwise
                    error("Wrong align_to input format.")
            end

            fprintf("\nTest head angle trace ...\n");
            fprintf("\n****** %s ******\n", align_to);

            trace_this = obj.("AngleHeadTrace"+string(align_to));

            for f = 1:length(obj.MixedFPs)

                fprintf("\n*** %s ***\n", obj.ForePeriods(f));

                time_points = trace_this.("TimePoints_"+obj.ForePeriods(f));

                % Performance_PortChosen
                C_L = trace_this.("L_"+obj.ForePeriods(f)+"_Correct");
                C_R = trace_this.("R_"+obj.ForePeriods(f)+"_Correct");
                W_L = trace_this.("L_"+obj.ForePeriods(f)+"_Wrong");
                W_R = trace_this.("R_"+obj.ForePeriods(f)+"_Wrong");

                % 1. to discriminant between chosen port, 
                % between correct trials (chose L) and correct trials (chose R)
                fprintf("\nTest between correct trials (chose L) and correct trials (chose R) ...\n");
                cc_table = Traj.testTrace(time_points, C_L, C_R, shuffle_iters, alpha, "Chosen");
                Out.(obj.ForePeriods(f)+"_CC") = cc_table;
           
%                 % 2. to discriminant between chosen port,
%                 % 2.1. 
%                 fprintf("\nTest between correct trials (chose L) and wrong trials (chose R) ...\n");
%                 cw_lr_table = Traj.testTrace(time_points, C_L, W_R, shuffle_iters, alpha, "Chosen");
%                 Out.(obj.ForePeriods(f)+"_CW_LR") = cw_lr_table;
% 
%                 % 2.2. between correct trials (chose R) and wrong trials (chose L)
%                 fprintf("\nTest between correct trials (chose R) and wrong trials (chose L) ...\n");
%                 cw_rl_table = Traj.testTrace(time_points, W_L, C_R, shuffle_iters, alpha, "Chosen");
%                 Out.(obj.ForePeriods(f)+"_CW_RL") = cw_rl_table;
% 
%                 % 3. to discriminant between correct port, 
%                 % 3.1. between correct trials (chose L) and wrong trials (chose L)
%                 fprintf("\nTest between correct trials (chose L) and wrong trials (chose L) ...\n");
%                 cw_ll_table = Traj.testTrace(time_points, C_L, W_L, shuffle_iters, alpha, "Correct");
%                 Out.(obj.ForePeriods(f)+"_CW_LL") = cw_ll_table;
% 
%                 % 3.2. between correct trials (chose R) and wrong trials (chose R)
%                 fprintf("\nTest between correct trials (chose R) and wrong trials (chose R) ...\n");
%                 cw_rr_table = Traj.testTrace(time_points, C_R, W_R, shuffle_iters, alpha, "Correct");
%                 Out.(obj.ForePeriods(f)+"_CW_RR") = cw_rr_table;
            end
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
            set(fig, 'unit', 'centimeters', 'position', [2 1.5 19 21.5], 'paperpositionmode', 'auto', 'color', 'w');

            switch label
                case {"All"}
                    group_title = obj.ANM+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end);
                otherwise
                    group_title = obj.ANM+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end) + " ("+label+")";
            end
            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.05 0.95 0.9 0.04],...
                'string', group_title, 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            % Chose L
            h1 = 1.5;
            ax1 = axes; colormap(mycolormap);
            set(ax1, 'units', 'centimeters', 'position', [1.5 h1, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax1, obj, "PortChosen", "L", "Performance", "Wrong", "AlignTo", "In", "Label", label);
            ax1.Title.String = [];

            h2 = h1 + ax1.Position(4) + .2;
            ax2 = axes; colormap(mycolormap);
            set(ax2, 'units', 'centimeters', 'position', [1.5 h2, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax2, obj, "PortChosen", "L", "Performance", "Correct", "AlignTo", "In", "Label", label);
            set(ax2, 'xcolor', 'none')

            ax11 = axes; colormap(mycolormap);
            set(ax11, 'units', 'centimeters', 'position', [9.5 h1, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax11, obj, "PortChosen", "L", "Performance", "Wrong", "AlignTo", "Out", "Label", label);
            set(ax11, 'ycolor', 'none');
            ax11.Position(4) = ax1.Position(4);
            ax11.Title.String = [];

            ax21 = axes; colormap(mycolormap);
            set(ax21, 'units', 'centimeters', 'position', [9.5 h2, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax21, obj, "PortChosen", "L", "Performance", "Correct", "AlignTo", "Out", "Label", label);
            set(ax21, 'ycolor', 'none');
            ax21.Position(4) = ax2.Position(4);
            set(ax21, 'xcolor', 'none')
            
            cb = colorbar(ax2, "Units", "centimeters", "Position", [ax11.Position(1)+ax11.Position(3)+0.5 ax1.Position(2) .3 ax2.Position(4)+ax1.Position(4)+.2]);
            cb.Label.String = "Head angle (°)";
            cb.Label.FontSize = 9;
            cb.Label.FontWeight = "Bold";
            cb.Ticks = -90:30:90;

            % Chose R
            h3 = h2 + ax2.Position(4) + .8;
            ax3 = axes("Parent", fig); colormap(mycolormap);
            set(ax3, 'units', 'centimeters', 'position', [1.5 h3 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax3, obj, "PortChosen", "R", "Performance", "Wrong", "AlignTo", "In", "Label", label);
            set(ax3, "xticklabel", []);
            ax3.XLabel.String = "";
            ax3.Title.String = "";

            h4 = h3 + ax3.Position(4) + .2;
            ax4 = axes("Parent", fig); colormap(mycolormap);
            set(ax4, 'units', 'centimeters', 'position', [1.5 h4 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax4, obj, "PortChosen", "R", "Performance", "Correct", "AlignTo", "In", "Label", label);
            set(ax4, "xcolor", 'none');

            ax31 = axes; colormap(mycolormap);
            set(ax31, 'units', 'centimeters', 'position', [9.5 h3 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax31, obj, "PortChosen", "R", "Performance", "Wrong", "AlignTo", "Out", "Label", label);
            set(ax31, 'ycolor', 'none');
            set(ax31, "xticklabel", []);
            ax31.XLabel.String = "";
            ax31.Title.String = "";
            ax31.Position(4) = ax3.Position(4);

            ax41 = axes("Parent", fig); colormap(mycolormap);
            set(ax41, 'units', 'centimeters', 'position', [9.5 h4 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax41, obj, "PortChosen", "R", "Performance", "Correct", "AlignTo", "Out", "Label", label);
            set(ax41, "xcolor", 'none', 'ycolor', 'none');

            % Late
            h5 = h4 + ax4.Position(4) + .8;
            ax5 = axes; colormap(mycolormap);
            set(ax5, 'units', 'centimeters', 'position', [1.5 h5 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax5, obj, "Performance", "Late", "AlignTo", "In", "Label", label);
            set(ax5, "xticklabel", []); ax5.XLabel.String = "";

            ax51 = axes; colormap(mycolormap);
            set(ax51, 'units', 'centimeters', 'position', [9.5 h5, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax51, obj, "Performance", "Late", "AlignTo", "Out", "Label", label);
            set(ax51, 'ycolor', 'none');
            set(ax51, "xticklabel", []); ax51.XLabel.String = ""; 
            ax51.Position(4) = ax5.Position(4);

            % Premature
            h6 = h5 + ax5.Position(4) + .8;
            ax6 = axes; colormap(mycolormap);
            set(ax6, 'units', 'centimeters', 'position', [1.5 h6, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax6, obj, "Performance", "Premature", "AlignTo", "In", "Label", label);
            set(ax6, "xticklabel", []); ax6.XLabel.String = "";

            ax61 = axes; colormap(mycolormap);
            set(ax61, 'units', 'centimeters', 'position', [9.5 h6, 7 4.5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax61, obj, "Performance", "Premature", "AlignTo", "Out", "Label", label);
            set(ax61, 'ycolor', 'none');
            set(ax61, "xticklabel", []); ax61.XLabel.String = ""; 
            ax61.Position(4) = ax6.Position(4);

            h_fig = h6 + ax6.Position(4) + 1.4;

%             cLimits = max([clim(ax1) clim(ax2) clim(ax3) clim(ax4)]);
            cLimits = 90;
            clim(ax1, [-1 1]*cLimits); clim(ax11, [-1 1]*cLimits);
            clim(ax2, [-1 1]*cLimits); clim(ax21, [-1 1]*cLimits);
            clim(ax3, [-1 1]*cLimits); clim(ax31, [-1 1]*cLimits);
            clim(ax4, [-1 1]*cLimits); clim(ax41, [-1 1]*cLimits);
            clim(ax5, [-1 1]*cLimits); clim(ax51, [-1 1]*cLimits);
            clim(ax6, [-1 1]*cLimits); clim(ax61, [-1 1]*cLimits);

            fig.Position(4) = h_fig;
            
        end
        
        function fig = plotTrace(obj)

            fig = figure(34); clf(34);
            set(fig, 'unit', 'centimeters', 'position', [2 2 24 16.3], 'paperpositionmode', 'auto', 'color', 'w');

            group_title = obj.ANM+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end);
            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.2 0.95 0.6 0.04],...
                'string', group_title, 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            h1 = 1.5;
            ax1 = axes();
            set(ax1, 'units', 'centimeters', 'position', [1.5 1.5 19/2 4], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHead(ax1, obj, "Long", "In");

            h2 = h1 + ax1.Position(4) + .8;
            ax2 = axes();
            set(ax2, 'units', 'centimeters', 'position', [1.5 h2 14/2 4], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHead(ax2, obj, "Med", "In");
            set(ax2, 'xticklabel', []); ax2.XLabel.String = "";

            h3 = h2 + ax2.Position(4) + .8;
            ax3 = axes();
            set(ax3, 'units', 'centimeters', 'position', [1.5 h3 9/2 4], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHead(ax3, obj, "Short", "In");
            set(ax3, 'xticklabel', []); ax3.XLabel.String = "";

            w1 = 1.5 + ax1.Position(3) + 2;
            ax4 = axes();
            set(ax4, 'units', 'centimeters', 'position', [w1 h1 19/2 4], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHead(ax4, obj, "Long", "Out");
            set(ax4, 'yaxislocation', 'right');

            w2 = w1 + (ax1.Position(3) - ax2.Position(3));
            ax5 = axes();
            set(ax5, 'units', 'centimeters', 'position', [w2 h2 14/2 4], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHead(ax5, obj, "Med", "Out");
            set(ax5, 'yaxislocation', 'right', 'xticklabel', []); ax5.XLabel.String = "";

            w3 = w1 + (ax1.Position(3) - ax3.Position(3));
            ax6 = axes();
            set(ax6, 'units', 'centimeters', 'position', [w3 h3 9/2 4], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHead(ax6, obj, "Short", "Out");
            set(ax6, 'yaxislocation', 'right', 'xticklabel', []); ax6.XLabel.String = "";


        end

    end
end