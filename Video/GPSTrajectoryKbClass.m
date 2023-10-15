classdef GPSTrajectoryKbClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        DLCTrackingOutFile
        ANMInfoFile

        Session
        ANM
        Treatment
        Dose
        Label
        Experimenter
        Task

        DLCTracking
        BehTable

        AngleHeadTraceIn
        AngleHeadTraceOut

        AngleHeadTraceInTest
        AngleHeadTraceOutTest
    end

    properties (Constant)

        Condition = ["Cue", "Uncue"];
        CueUncue = [1 0];
        Ports = ["L", "R"];
        TimePointsIn  = -99:1:2500;
        TimePointsOut = -1599:1:1000;

    end

    properties (Dependent)

        TrialInfo
        NumTrials
        Trial
        Stage
        Performance
        PortCorrect
        PortChosen
        SessionFP
        FP
        RT
        HD
        MT
        Cued

        Ind
        TimeFromIn
        TimeFromOut

        PortVec

        AngleHead
        AngleHeadMatIn
        AngleHeadMatOut

    end

    methods
        function obj = GPSTrajectoryKbClass(DLCTrackingOutFile, AnmInfoFile)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            obj.DLCTrackingOutFile = DLCTrackingOutFile;
            file_info = split(string(DLCTrackingOutFile), filesep);
            obj.Session = file_info(end-3);
            obj.ANM     = file_info(end-5);

            obj.ANMInfoFile = AnmInfoFile;

            ANMInfoTable = readtable(obj.ANMInfoFile, "Sheet", obj.ANM, "TextType", "string");
            SessionInfo  = ANMInfoTable(ANMInfoTable.Session==str2double(obj.Session), :);
            obj.Treatment    = SessionInfo.Treatment;
            obj.Dose         = SessionInfo.Dose;
            obj.Label        = SessionInfo.Label;
            obj.Experimenter = SessionInfo.Experimenter;
            obj.Task         = SessionInfo.Task;

            BehTableFile = dir(fullfile(SessionInfo.SessionFolder, "*SessionTable*"));
            obj.BehTable = readtable(fullfile(SessionInfo.SessionFolder, BehTableFile.name));

            load(DLCTrackingOutFile, "DLCTrackingOut");
            obj.DLCTracking = DLCTrackingOut;

            obj = obj.removeOddTrials;

            obj.AngleHeadTraceIn  = obj.getAngleHeadTrace("In");
            obj.AngleHeadTraceOut = obj.getAngleHeadTrace("Out");

        end

        %%
        function obj = removeOddTrials(obj)

            angle_in = cellfun(@(a, t) a(find(t==0, 1)), obj.AngleHead, obj.TimeFromIn);
            ind_odd  = find(angle_in > mean(angle_in)+5*std(angle_in) | angle_in < mean(angle_in)-5*std(angle_in));

            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];

            zero_point = cellfun(@(t) sum(t==0), obj.TimeFromIn);
            ind_odd = find(zero_point~=1);
            
            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];
        end

        %%
        function value = get.TrialInfo(obj)

            trial_info = table(repmat(obj.Session, obj.NumTrials, 1), repmat(obj.Treatment, obj.NumTrials, 1), repmat(obj.Dose, obj.NumTrials, 1), repmat(obj.Label, obj.NumTrials, 1), ...
                obj.Trial', obj.Stage', obj.Cued', obj.Performance', obj.PortCorrect', obj.PortChosen', obj.RT', obj.MT', obj.HD', ...
                'VariableNames', ["Session", "Treatment", "Dose", "Label", "Trial", "Stage", "Cued", "Performance", "PortCorrect", "PortChosen", "RT", "MT", "HD"]);

            value = trial_info;

        end

        function value = get.NumTrials(obj)

            num_trials = length(obj.Trial);

            value = num_trials;
        end

        function value = get.Trial(obj)

            trial = obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:);

            value = trial;
        end

        function value = get.Stage(obj)

            stage = obj.BehTable.Stage(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            stage = stage';

            value = stage;
        end

        function value = get.Performance(obj)

            perf = string(obj.DLCTracking.PoseTracking(1).Performance);
            perf(ismember(perf, ["LateMiss", "LateWrong", "LateCorrect"])) = "Late";

            value = perf;
        end

        function value = get.PortCorrect(obj)

            port_correct = obj.Ports(obj.BehTable.PortCorrect(obj.Trial));

            value = port_correct;
        end

        function value = get.PortChosen(obj)

            port_chosen  = strings(1, length(obj.Trial));
            chosen_ind   = ~isnan(obj.BehTable.PortChosen(obj.Trial));
            chosen_trial = obj.Trial(chosen_ind);
            port_chosen(chosen_ind) = obj.Ports(obj.BehTable.PortChosen(chosen_trial));

            value = port_chosen;
        end

        function value = get.FP(obj)

            fp = obj.BehTable.FP(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            fp = fp';

            value = fp;
        end

        function value = get.SessionFP(obj)

            value = unique(obj.FP);
        end

        function value = get.RT(obj)

            rt = obj.BehTable.RT(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            rt = rt';

            value = rt;
        end

        function value = get.HD(obj)

            hd = obj.BehTable.HoldDuration(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            hd = hd';

            value = hd;
        end

        function value = get.MT(obj)

            mt = obj.BehTable.MovementTime(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            mt = mt';

            value = mt;
        end

        function value = get.Cued(obj)

            cued = obj.BehTable.Cued(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            cued = cued';

            value = cued;
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
        function value = get.TimeFromIn(obj)

            time_from_in = cellfun(@(x) x(:, 4), obj.DLCTracking.PoseTracking(1).PosData, 'UniformOutput', false);

            value = time_from_in;
        end

        function value = get.TimeFromOut(obj)

            in2out = num2cell(1000 * (obj.BehTable.CentOutTime - obj.BehTable.CentInTime));
            in2out = in2out(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            time_from_out = cellfun(@(x, y) x - y, obj.TimeFromIn, in2out', 'UniformOutput', false);

            value = time_from_out;
        end

        function value = get.PortVec(obj)

            port_vec = arrayfun(@(x) x.R - x.L, obj.DLCTracking.PortLoc, 'UniformOutput', false);

            value = port_vec;

        end

        %%
        function value = get.AngleHead(obj)

            indL = find(strcmp("ear_base_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_base_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData;
            posR = obj.DLCTracking.PoseTracking(indR).PosData;

            head_vec = cellfun(@(x, y) y(:, 1:2) - x(:, 1:2), posL, posR, 'UniformOutput', false);

            head_angle = cellfun(@(x, y) calAngle(x, y), head_vec, obj.PortVec, 'UniformOutput', false);

            angle_sign = cellfun(@(x, y) 2*(x(:, 2)>=y(2))-1, head_vec, obj.PortVec, 'UniformOutput', false);

            head_angle = cellfun(@(x, y) x.*y, head_angle, angle_sign, 'UniformOutput', false);

            head_angle = cellfun(@(x) smoothdata(x, "gaussian", 5), head_angle, 'UniformOutput', false);

            value = head_angle;
        end

        function value = get.AngleHeadMatIn(obj)

            angle_head_mat_in = cellfun(@(a, t) interp1(t, a, obj.TimePointsIn, "linear"), obj.AngleHead, obj.TimeFromIn, 'UniformOutput', false);
            angle_head_mat_in = angle_head_mat_in';
            angle_head_mat_in = cell2mat(angle_head_mat_in);

            value = angle_head_mat_in;
        end

        function value = get.AngleHeadMatOut(obj)

            angle_head_mat_out = cellfun(@(a, t) interp1(t, a, obj.TimePointsOut, "linear"), obj.AngleHead, obj.TimeFromOut, 'UniformOutput', false);
            angle_head_mat_out = angle_head_mat_out';
            angle_head_mat_out = cell2mat(angle_head_mat_out);

            value = angle_head_mat_out;
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
                        time_bin_edges  = floor(obj.TimePointsIn(1)):20:1000*obj.SessionFP+300;
                        time_bin_center = time_bin_edges(1:end-1)+10;
                    case {'Out'}
                        time_bin_edges = -obj.SessionFP*1000-100:20:300;
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
        function save(obj, targetDir)

            [savepath, ~] = fileparts(obj.DLCTrackingOutFile);
            save(fullfile(savepath, "GPSTrajectoryClass_"+obj.Task+"_"+upper(obj.ANM)+"_"+obj.Session), 'obj');

            if nargin==2
                save(fullfile(targetDir, "GPSTrajectoryClass_"+obj.Task+"_"+upper(obj.ANM)+"_"+obj.Session), 'obj');
            end
        end

        %% Plots
        function print(obj, Func, targetDir)

            if nargin==1
                Func = "HeatMap";
            end

            switch lower(Func)
                case {'heatmap'}
                    hf = obj.plotHeatMap();
                case {'trace'}
                    hf = obj.plotTrace();
            end

            [savepath, ~] = fileparts(obj.DLCTrackingOutFile);
            savename = fullfile(savepath, "GPSTrajectoryClass_" + obj.Task + "_" + string(Func) + "_" + upper(obj.ANM) + "_" + obj.Session);
            print(hf, '-dpdf', savename, '-bestfit')
            print(hf, '-dpng', savename)
            saveas(hf, savename, 'fig')

            if nargin==3
                % check if targetDir exists
                if ~contains(targetDir, '/') && ~contains(targetDir, '\')
                    % so it is a relative path
                    if ~exist(targetDir, 'dir')
                        mkdir(targetDir)
                    end
                end
                savename = fullfile(targetDir, "GPSTrajectoryClass_" + obj.Task + "_" + string(Func) + "_" + upper(obj.ANM) + "_" + obj.Session);
                print(hf, '-dpdf', savename, '-bestfit')
                print(hf, '-dpng', savename)
                saveas(hf, savename, 'fig')
            end

        end

        function fig = plotHeatMap(obj)

            mycolormap = customcolormap_preset("red-white-blue");

            fig = figure(33); clf(33);
            set(fig, 'unit', 'centimeters', 'position', [2 2 19 25], 'paperpositionmode', 'auto', 'color', 'w');

            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.05 0.9 0.9 0.08],...
                'string', obj.ANM+" / "+obj.Session+" / "+obj.Task+" / "+char(obj.Treatment), 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            % Target L
            w_l = 1.5;
            ax_l = axes; colormap(mycolormap);
            set(ax_l, 'units', 'centimeters', 'position', [w_l 1.5, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHeadKb(ax_l, obj, "PortCorrect", "L", "Performance", "All", "AlignTo", "In");

            % Target R
            w_r = w_l + ax_l.Position(3) + .8;
            ax_r = axes; colormap(mycolormap);
            set(ax_r, 'units', 'centimeters', 'position', [w_r 1.5 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHeadKb(ax_r, obj, "PortCorrect", "R", "Performance", "All", "AlignTo", "In");
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
            set(fig, 'unit', 'centimeters', 'position', [2 2 8*obj.SessionFP+2 7.5], 'paperpositionmode', 'auto', 'color', 'w');

            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.05 0.9 0.9 0.08],...
                'string', obj.ANM+" / "+obj.Session+" / "+obj.Task+" / "+char(obj.Treatment), 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            ax1 = axes();
            set(ax1, 'units', 'centimeters', 'position', [1.5 1.5 8*obj.SessionFP 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.traceAngleHeadKb(ax1, obj, "In");

        end

    end
end