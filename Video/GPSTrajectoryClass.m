classdef GPSTrajectoryClass
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
        ForePeriods = ["Short", "Med", "Long"];
        MixedFPs    = [.5 1 1.5];
        Ports = ["L", "R"];

        TimeMatIn  = -100:10:2500;
        TimeMatOut = -1600:10:1000;
    end

    properties (Dependent)
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

        Ind
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

        AngleHeadMatIn
        AngleHeadMatOut
        AngSpeedHeadMatIn
        AngSpeedHeadMatOut
        AngAccHeadMatIn
        AngAccHeadMatOut

        PosXHeadMatIn
        PosXHeadMatOut
        PosYHeadMatIn
        PosYHeadMatOut
        SpeedXHeadMatIn
        SpeedXHeadMatOut
        SpeedYHeadMatIn
        SpeedYHeadMatOut
        AccXHeadMatIn
        AccXHeadMatOut
        AccYHeadMatIn
        AccYHeadMatOut

        SpeedHeadMatIn
        SpeedHeadMatOut
        SpeedDirHeadMatIn
        SpeedDirHeadMatOut
        AccHeadMatIn
        AccHeadMatOut
        dPhiHeadMatIn
        dPhiHeadMatOut
    end

    methods
        function obj = GPSTrajectoryClass(DLCTrackingOutFile, AnmInfoFile)
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

%             obj.AngleHeadTraceIn  = obj.getAngleHeadTrace("In");
%             obj.AngleHeadTraceOut = obj.getAngleHeadTrace("Out");
            
        end

        %%
        function obj = removeOddTrials(obj)
            
            angle_in = cellfun(@(a, t) a(t==0), obj.AngleHead, obj.TimeFromIn);
            ind_odd  = find(angle_in > mean(angle_in)+5*std(angle_in) | angle_in < mean(angle_in)-5*std(angle_in));

            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];

            time_error = cellfun(@(t) sum(diff(t)<=0), obj.TimeFromIn);
            ind_odd = find(time_error>0);
            
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
                obj.Trial', obj.Stage', obj.Performance', obj.PortCorrect', obj.PortChosen', obj.FP', obj.RT', obj.MT', obj.HD', ...
                'VariableNames', ["Session", "Treatment", "Dose", "Label", "Trial", "Stage", "Performance", "PortCorrect", "PortChosen", "FP", "RT", "MT", "HD"]);

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

        function value = get.TimeWarped(obj)
            
            in2out = num2cell(1000 * (obj.BehTable.CentOutTime - obj.BehTable.CentInTime));
            in2out = in2out(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            time_warped = cellfun(@(x, y) x ./ y, obj.TimeFromIn, in2out', 'UniformOutput', false);

            value = time_warped;
        end

        function value = get.PortVec(obj)
            
            port_vec = arrayfun(@(x) x.R - x.L, obj.DLCTracking.PortLoc, 'UniformOutput', false);
            value = port_vec;
        end

        function value = get.PortCent(obj)

            port_cent = arrayfun(@(x) (x.R + x.L) / 2, obj.DLCTracking.PortLoc, 'UniformOutput', false);
            value = port_cent;
        end

        %% Angle of head
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
            value = obj.alignMatrix(obj.AngleHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.AngleHeadMatOut(obj)
            value = obj.alignMatrix(obj.AngleHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Angle speed of head
        function value = get.AngSpeedHead(obj)
            
            d_a = cellfun(@(a) a([2:end end]) - a([1 1:end-1]), obj.AngleHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            ang_speed_head = cellfun(@(da, dt) da ./ dt, d_a, d_t, 'UniformOutput', false);
            ang_speed_head = cellfun(@(x) smoothdata(x, "gaussian", 5), ang_speed_head, 'UniformOutput', false);

            value = ang_speed_head;
        end

        function value = get.AngSpeedHeadMatIn(obj)
            value = obj.alignMatrix(obj.AngSpeedHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.AngSpeedHeadMatOut(obj)
            value = obj.alignMatrix(obj.AngSpeedHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Angle acceleration of head
        function value = get.AngAccHead(obj)
            
            d_a = cellfun(@(a) a([2:end end]) - a([1 1:end-1]), obj.AngleHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            ang_speed_head = cellfun(@(da, dt) da ./ dt, d_a, d_t, 'UniformOutput', false);

            d_d_v = cellfun(@(v) v([2:end end]) - v([1 1:end-1]), ang_speed_head, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            ang_acc_head = cellfun(@(ddv, dt) ddv ./ dt, d_d_v, d_t, 'UniformOutput', false);
            ang_acc_head = cellfun(@(x) smoothdata(x, "gaussian", 5), ang_acc_head, 'UniformOutput', false);

            value = ang_acc_head;
        end

        function value = get.AngAccHeadMatIn(obj)
            value = obj.alignMatrix(obj.AngAccHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.AngAccHeadMatOut(obj)
            value = obj.alignMatrix(obj.AngAccHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Position of head (X)
        function value = get.PosXHead(obj)

            indL = find(strcmp("ear_base_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_base_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData;
            posR = obj.DLCTracking.PoseTracking(indR).PosData;

            head_pos_x = cellfun(@(x, y, z) (y(:, 1) + x(:, 1))/2 - z(1), posL, posR, obj.PortCent, 'UniformOutput', false);
            head_pos_x = cellfun(@(x) smoothdata(x, "gaussian", 5), head_pos_x, 'UniformOutput', false);

            value = head_pos_x;
        end

        function value = get.PosXHeadMatIn(obj)
            value = obj.alignMatrix(obj.PosXHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.PosXHeadMatOut(obj)
            value = obj.alignMatrix(obj.PosXHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Position of head (Y)
        function value = get.PosYHead(obj)

            indL = find(strcmp("ear_base_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_base_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData;
            posR = obj.DLCTracking.PoseTracking(indR).PosData;

            head_pos_y = cellfun(@(x, y, z) (y(:, 2) + x(:, 2))/2 - z(2), posL, posR, obj.PortCent, 'UniformOutput', false);
            head_pos_y = cellfun(@(x) smoothdata(x, "gaussian", 5), head_pos_y, 'UniformOutput', false);

            value = head_pos_y;
        end

        function value = get.PosYHeadMatIn(obj)
            value = obj.alignMatrix(obj.PosYHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.PosYHeadMatOut(obj)
            value = obj.alignMatrix(obj.PosYHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Speed (X) of head
        function value = get.SpeedXHead(obj)
            
            d_x = cellfun(@(x) x([2:end end]) - x([1 1:end-1]), obj.PosXHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            speed_x_head = cellfun(@(dx, dt) dx ./ dt, d_x, d_t, 'UniformOutput', false);
            speed_x_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_x_head, 'UniformOutput', false);

            value = speed_x_head;
        end

        function value = get.SpeedXHeadMatIn(obj)
            value = obj.alignMatrix(obj.SpeedXHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.SpeedXHeadMatOut(obj)
            value = obj.alignMatrix(obj.SpeedXHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Speed (Y) of head
        function value = get.SpeedYHead(obj)
            
            d_y = cellfun(@(y) y([2:end end]) - y([1 1:end-1]), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            speed_y_head = cellfun(@(dx, dt) dx ./ dt, d_y, d_t, 'UniformOutput', false);
            speed_y_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_y_head, 'UniformOutput', false);

            value = speed_y_head;
        end

        function value = get.SpeedYHeadMatIn(obj)
            value = obj.alignMatrix(obj.SpeedYHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.SpeedYHeadMatOut(obj)
            value = obj.alignMatrix(obj.SpeedYHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Acceleration (X) of head
        function value = get.AccXHead(obj)
            
            d_x = cellfun(@(x) x([2:end end]) - x([1 1:end-1]), obj.PosXHead, 'UniformOutput', false);
            dd_x = cellfun(@(dx) dx([2:end end]) - dx([1 1:end-1]), d_x, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            acc_x_head = cellfun(@(ddx, dt) ddx ./ (dt.^2), dd_x, d_t, 'UniformOutput', false);
            acc_x_head = cellfun(@(x) smoothdata(x, "gaussian", 5), acc_x_head, 'UniformOutput', false);

            value = acc_x_head;
        end

        function value = get.AccXHeadMatIn(obj)
            value = obj.alignMatrix(obj.AccXHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.AccXHeadMatOut(obj)
            value = obj.alignMatrix(obj.AccXHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Acceleration (Y) of head
        function value = get.AccYHead(obj)
            
            d_y = cellfun(@(y) y([2:end end]) - y([1 1:end-1]), obj.PosYHead, 'UniformOutput', false);
            dd_y = cellfun(@(dy) dy([2:end end]) - dy([1 1:end-1]), d_y, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            acc_y_head = cellfun(@(ddy, dt) ddy ./ (dt.^2), dd_y, d_t, 'UniformOutput', false);
            acc_y_head = cellfun(@(x) smoothdata(x, "gaussian", 5), acc_y_head, 'UniformOutput', false);

            value = acc_y_head;
        end

        function value = get.AccYHeadMatIn(obj)
            value = obj.alignMatrix(obj.AccYHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.AccYHeadMatOut(obj)
            value = obj.alignMatrix(obj.AccYHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Speed of head
        function value = get.SpeedHead(obj)
            
            d_x = cellfun(@(x) x([2:end end]) - x([1 1:end-1]), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) y([2:end end]) - y([1 1:end-1]), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            speed_head = cellfun(@(dx, dy, dt) vecnorm([dx, dy], 2, 2) ./ dt, d_x, d_y, d_t, 'UniformOutput', false);
            speed_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_head, 'UniformOutput', false);

            value = speed_head;
        end

        function value = get.SpeedHeadMatIn(obj)
            value = obj.alignMatrix(obj.SpeedHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.SpeedHeadMatOut(obj)
            value = obj.alignMatrix(obj.SpeedHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Speed direction of head
        function value = get.SpeedDirHead(obj)
            
            d_x = cellfun(@(x) x([2:end end]) - x([1 1:end-1]), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) y([2:end end]) - y([1 1:end-1]), obj.PosYHead, 'UniformOutput', false);

            speed_vec = cellfun(@(dx, dy) [dx, dy], d_x, d_y, 'UniformOutput', false);
            port_vec_vert = cellfun(@(pv) [1 -1] .* [pv(2) pv(1)], obj.PortVec, 'UniformOutput', false);

            speed_dir_head = cellfun(@(x, y) calAngle(x, y), speed_vec, port_vec_vert, 'UniformOutput', false);

            angle_sign = cellfun(@(x, y) 2*(x(:, 1)>=y(1))-1, speed_vec, port_vec_vert, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x, y) x.*y, speed_dir_head, angle_sign, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_dir_head, 'UniformOutput', false);

            value = speed_dir_head;
        end

        function value = get.SpeedDirHeadMatIn(obj)
            value = obj.alignMatrix(obj.SpeedDirHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.SpeedDirHeadMatOut(obj)
            value = obj.alignMatrix(obj.SpeedDirHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Acceleration of head
        function value = get.AccHead(obj)
            
            d_x = cellfun(@(x) x([2:end end]) - x([1 1:end-1]), obj.PosXHead, 'UniformOutput', false);
            d_d_x = cellfun(@(dx) dx([2:end end]) - dx([1 1:end-1]), d_x, 'UniformOutput', false);
            d_y = cellfun(@(y) y([2:end end]) - y([1 1:end-1]), obj.PosYHead, 'UniformOutput', false);
            d_d_y = cellfun(@(dy) dy([2:end end]) - dy([1 1:end-1]), d_y, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            acc_head = cellfun(@(ddx, ddy, dt) vecnorm([ddx, ddy], 2, 2) ./ (dt.^2), d_d_x, d_d_y, d_t, 'UniformOutput', false);
            acc_head = cellfun(@(x) smoothdata(x, "gaussian", 5), acc_head, 'UniformOutput', false);

            value = acc_head;
        end

        function value = get.AccHeadMatIn(obj)
            value = obj.alignMatrix(obj.AccHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.AccHeadMatOut(obj)
            value = obj.alignMatrix(obj.AccHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %% Change of speed direction (Phi) of head
        function value = get.dPhiHead(obj)
            
            d_x = cellfun(@(x) x([2:end end]) - x([1 1:end-1]), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) y([2:end end]) - y([1 1:end-1]), obj.PosYHead, 'UniformOutput', false);

            speed_vec = cellfun(@(dx, dy) [dx, dy], d_x, d_y, 'UniformOutput', false);
            port_vec_vert = cellfun(@(pv) [1 -1] .* [pv(2) pv(1)], obj.PortVec, 'UniformOutput', false);

            speed_dir_head = cellfun(@(x, y) calAngle(x, y), speed_vec, port_vec_vert, 'UniformOutput', false);

            angle_sign = cellfun(@(x, y) 2*(x(:, 1)>=y(1))-1, speed_vec, port_vec_vert, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x, y) x.*y, speed_dir_head, angle_sign, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x) unwrap(x, 180), speed_dir_head, 'UniformOutput', false);

            d_dir = cellfun(@(dir) dir([2:end end]) - dir([1 1:end-1]), speed_dir_head, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);
            dphi_head = cellfun(@(ddir, dt) ddir ./ (dt.^2), d_dir, d_t, 'UniformOutput', false);
            dphi_head = cellfun(@(x) smoothdata(x, "gaussian", 5), dphi_head, 'UniformOutput', false);

            value = dphi_head;
        end

        function value = get.dPhiHeadMatIn(obj)
            value = obj.alignMatrix(obj.dPhiHead, obj.TimeFromIn, obj.TimeMatIn);
        end

        function value = get.dPhiHeadMatOut(obj)
            value = obj.alignMatrix(obj.dPhiHead, obj.TimeFromOut, obj.TimeMatOut);
        end

        %%
        function M = alignMatrix(~, trace, time_points, time_matrix)
            
            M = cellfun(@(x, t) interp1(t, x, time_matrix, "linear"), trace, time_points, 'UniformOutput', false);
            M = M';
            M = cell2mat(M);
        end

        function trace_normalized = normalizeTrace(~, trace, norm_method)
            
            if nargin<3
                norm_method = 'range';
            end
            trace_all = cell2mat(trace');
            trace_all_normalized = normalize(trace_all, norm_method);

            trace_normalized = mat2cell(trace_all_normalized, cellfun(@(x) length(x), trace));
            trace_normalized = trace_normalized';
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
                        time_bin_edges = -fp*1000-100:20:300;
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

        function Out = testTrace(obj, align_to, shuffle_iters)
            
            switch align_to
                case {'In', 'in', 'IN'}
                    align_to = "In";
                    TimeBins = [0 250; 250 500; 500 750; 750 1000; 1000 1250; 1250 1500];
                case {'Out', 'out', 'OUT'}
                    align_to = "Out";
                    TimeBins = fliplr(-1*[0 250; 250 500; 500 750; 750 1000; 1000 1250; 1250 1500]);
                otherwise
                    error("Wrong align_to input format.")
            end

            for f = 1:length(obj.MixedFPs)

                % 1. to discriminant between chosen port, 
                % between correct trials (chose L) and correct trials (chose R)
                indL = obj.Ind.ChooseL & obj.Ind.(obj.ForePeriods(f)) & obj.Ind.Correct;
                indR = obj.Ind.ChooseR & obj.Ind.(obj.ForePeriods(f)) & obj.Ind.Correct;
                TimeL = obj.("TimeFrom"+align_to)(indL);
                TimeR = obj.("TimeFrom"+align_to)(indR);
                AngleL = obj.AngleHead(indL);
                AngleR = obj.AngleHead(indR);

                disp(obj.ForePeriods(f));
                switch obj.ForePeriods(f)
                    case {'Short'}
                        time_bins = TimeBins(1:2, :);
                    case {'Med'}
                        time_bins = TimeBins(1:4, :);
                    case {'Long'}
                        time_bins = TimeBins;
                end

                YL = repmat("L", length(TimeL), 1);
                YR = repmat("R", length(TimeR), 1);
                Y  = [YL; YR];

                Y_shuffle = strings(length(Y), shuffle_iters);
                for s = 1:shuffle_iters
                    Y_shuffle(:, s) = Y(randperm(length(Y)));
                end

                auc = zeros(size(time_bins, 1), 1);
                pval = zeros(size(time_bins, 1), 1);
                for t = 1:size(time_bins, 1)
                    XL = cellfun(@(ang, time) mean(ang(time>=time_bins(t, 1) & time<time_bins(t, 2)), 'omitnan'), AngleL, TimeL);
                    XR = cellfun(@(ang, time) mean(ang(time>=time_bins(t, 1) & time<time_bins(t, 2)), 'omitnan'), AngleR, TimeR);
                    X  = [XL'; XR'];

                    roc_this = rocmetrics(Y, X, "L");
                    auc(t) = roc_this.AUC;
                    auc_this = 0.5 + abs(roc_this.AUC - 0.5);

                    auc_shuffle = zeros(1, shuffle_iters);
                    for s = 1:shuffle_iters
                        roc_shuffle = rocmetrics(Y_shuffle(:, s), X, "L");
                        auc_shuffle(s) = 0.5 + abs(roc_shuffle.AUC - 0.5);
                    end
                    pval_this = sum(auc_shuffle >= auc_this) / shuffle_iters;
                    pval(t) = pval_this;
                    disp(t);
                end

                Out.(obj.ForePeriods(f)) = table(time_bins(:, 1), time_bins(:, 2), auc, pval, ...
                    'VariableNames', {'TimeBeg', 'TimeEnd', 'AUC', 'pVal'});
                if align_to=="Out"
                    Out.(obj.ForePeriods(f)) = flip(Out.(obj.ForePeriods(f)));
                end
            end
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
            set(fig, 'unit', 'centimeters', 'position', [2 2 19 21.2], 'paperpositionmode', 'auto', 'color', 'w');

            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.05 0.95 0.9 0.04],...
                'string', obj.ANM+" / "+obj.Session+" / "+obj.Task+" / "+char(obj.Treatment), 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

            % Chose L
            h1 = 1.5;
            ax1 = axes; colormap(mycolormap);
            set(ax1, 'units', 'centimeters', 'position', [1.5 h1, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax1, obj, "PortChosen", "L", "Performance", "Wrong", "AlignTo", "In");
            ax1.Title.String = [];

            h2 = h1 + ax1.Position(4) + .2;
            ax2 = axes; colormap(mycolormap);
            set(ax2, 'units', 'centimeters', 'position', [1.5 h2, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax2, obj, "PortChosen", "L", "Performance", "Correct", "AlignTo", "In");
            set(ax2, 'xcolor', 'none')

            ax11 = axes; colormap(mycolormap);
            set(ax11, 'units', 'centimeters', 'position', [9.5 h1, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax11, obj, "PortChosen", "L", "Performance", "Wrong", "AlignTo", "Out");
            set(ax11, 'ycolor', 'none');
            ax11.Position(4) = ax1.Position(4);
            ax11.Title.String = [];

            ax21 = axes; colormap(mycolormap);
            set(ax21, 'units', 'centimeters', 'position', [9.5 h2, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax21, obj, "PortChosen", "L", "Performance", "Correct", "AlignTo", "Out");
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
            set(ax3, 'units', 'centimeters', 'position', [1.5 h3 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax3, obj, "PortChosen", "R", "Performance", "Wrong", "AlignTo", "In");
            set(ax3, "xticklabel", []);
            ax3.XLabel.String = "";
            ax3.Title.String = "";

            h4 = h3 + ax3.Position(4) + .2;
            ax4 = axes("Parent", fig); colormap(mycolormap);
            set(ax4, 'units', 'centimeters', 'position', [1.5 h4 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax4, obj, "PortChosen", "R", "Performance", "Correct", "AlignTo", "In");
            set(ax4, "xcolor", 'none');

            ax31 = axes; colormap(mycolormap);
            set(ax31, 'units', 'centimeters', 'position', [9.5 h3 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax31, obj, "PortChosen", "R", "Performance", "Wrong", "AlignTo", "Out");
            set(ax31, 'ycolor', 'none');
            set(ax31, "xticklabel", []);
            ax31.XLabel.String = "";
            ax31.Title.String = "";
            ax31.Position(4) = ax3.Position(4);

            ax41 = axes("Parent", fig); colormap(mycolormap);
            set(ax41, 'units', 'centimeters', 'position', [9.5 h4 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax41, obj, "PortChosen", "R", "Performance", "Correct", "AlignTo", "Out");
            set(ax41, "xcolor", 'none', 'ycolor', 'none');

            % Late
            h5 = h4 + ax4.Position(4) + .8;
            ax5 = axes; colormap(mycolormap);
            set(ax5, 'units', 'centimeters', 'position', [1.5 h5 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax5, obj, "Performance", "Late", "AlignTo", "In");
            set(ax5, "xticklabel", []); ax5.XLabel.String = "";

            ax51 = axes; colormap(mycolormap);
            set(ax51, 'units', 'centimeters', 'position', [9.5 h5, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax51, obj, "Performance", "Late", "AlignTo", "Out");
            set(ax51, 'ycolor', 'none');
            set(ax51, "xticklabel", []); ax51.XLabel.String = ""; 
            ax51.Position(4) = ax5.Position(4);

            % Premature
            h6 = h5 + ax5.Position(4) + .8;
            ax6 = axes; colormap(mycolormap);
            set(ax6, 'units', 'centimeters', 'position', [1.5 h6, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax6, obj, "Performance", "Premature", "AlignTo", "In");
            set(ax6, "xticklabel", []); ax6.XLabel.String = "";

            ax61 = axes; colormap(mycolormap);
            set(ax61, 'units', 'centimeters', 'position', [9.5 h6, 7 5], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
            Plot.heatmapAngleHead(ax61, obj, "Performance", "Premature", "AlignTo", "Out");
            set(ax61, 'ycolor', 'none');
            set(ax61, "xticklabel", []); ax61.XLabel.String = ""; 
            ax61.Position(4) = ax6.Position(4);

            h_fig = h6 + ax6.Position(4) + 1.3;

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
            set(fig, 'unit', 'centimeters', 'position', [2 2 24 16], 'paperpositionmode', 'auto', 'color', 'w');

            uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.2 0.95 0.6 0.04],...
                'string', obj.ANM+" / "+obj.Session+" / "+obj.Task+" / "+char(obj.Treatment), 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

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