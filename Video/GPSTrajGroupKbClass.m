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
        Trial
        Stage
        Performance
        PortCorrect
        PortChosen

        TaskFP
        FP
        RT
        HD
        MT
        Cued

        TimeFromIn
        TimeFromOut

        PortVec
        AngleHead

        AngleHeadMatIn
        AngleHeadMatOut

        AngleHeadTraceIn
        AngleHeadTraceOut
        AngleHeadTraceInTest
        AngleHeadTraceOutTest
    end

    properties (Constant)
        Condition = ["Cue", "Uncue"];
        CueUncue    = [1 0];
        Ports = ["L", "R"];
        TimePointsIn  = -99:1:2500;
        TimePointsOut = -1599:1:1000;
        ShuffleIters  = 200;
        Alpha         = 0.01;
    end

    properties (Dependent)
        Ind
    end

    methods
        function obj = GPSTrajGroupKbClass(TrajClassAll, AnmInfoFile)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.ANMInfoFile     = AnmInfoFile;
            obj.ANM             = unique(arrayfun(@(x) x.obj.ANM, TrajClassAll));
            obj.Sessions        = arrayfun(@(x) x.obj.Session, TrajClassAll);
            obj.Task            = unique(arrayfun(@(x) x.obj.Task, TrajClassAll));
            obj.NumSessions     = length(TrajClassAll);
            obj.Experimenter    = arrayfun(@(x) x.obj.Experimenter, TrajClassAll);

            obj.Treatment       = arrayfun(@(x) x.obj.Treatment, TrajClassAll);
            obj.Dose            = arrayfun(@(x) x.obj.Dose, TrajClassAll);
            obj.Label           = arrayfun(@(x) x.obj.Label, TrajClassAll);
            obj.NumTrials       = arrayfun(@(x) x.obj.NumTrials, TrajClassAll);

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
            
            obj.TaskFP         = unique(arrayfun(@(x) x.obj.SessionFP, TrajClassAll));

            allFPs             = arrayfun(@(x) x.obj.FP, TrajClassAll, 'UniformOutput', false);
            obj.FP             = cell2mat(allFPs');

            allHDs             = arrayfun(@(x) x.obj.HD, TrajClassAll, 'UniformOutput', false);
            obj.HD             = cell2mat(allHDs');

            allRTs             = arrayfun(@(x) x.obj.RT, TrajClassAll, 'UniformOutput', false);
            obj.RT             = cell2mat(allRTs');

            allMTs             = arrayfun(@(x) x.obj.MT, TrajClassAll, 'UniformOutput', false);
            obj.MT             = cell2mat(allMTs');

            allCued             = arrayfun(@(x) x.obj.Cued, TrajClassAll, 'UniformOutput', false);
            obj.Cued             = cell2mat(allCued');

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

            allPortVecs        = arrayfun(@(x) x.obj.PortVec, TrajClassAll, 'UniformOutput', false);
            obj.PortVec        = [];
            for i = 1:obj.NumSessions
                obj.PortVec = [obj.PortVec allPortVecs{i}];
            end

            allAngleHeads        = arrayfun(@(x) x.obj.AngleHead, TrajClassAll, 'UniformOutput', false);
            obj.AngleHead        = [];
            for i = 1:obj.NumSessions
                obj.AngleHead = [obj.AngleHead allAngleHeads{i}];
            end

            allAngleHeadMatIns   = arrayfun(@(x) x.obj.AngleHeadMatIn, TrajClassAll, 'UniformOutput', false);
            obj.AngleHeadMatIn   = cell2mat(allAngleHeadMatIns);
            
            allAngleHeadMatOuts  = arrayfun(@(x) x.obj.AngleHeadMatOut, TrajClassAll, 'UniformOutput', false);
            obj.AngleHeadMatOut  = cell2mat(allAngleHeadMatOuts);

            obj.AngleHeadTraceIn  = obj.getAngleHeadTrace("In");
            obj.AngleHeadTraceOut = obj.getAngleHeadTrace("Out");

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
                        time_bin_edges  = floor(obj.TimePointsIn(1)):20:1000*obj.TaskFP+300;
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