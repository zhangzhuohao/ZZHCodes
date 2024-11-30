classdef GPSTrajSessionClass < GPSTrajClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        SessionFolder

        DLCTrackingOutFile
        BehClassFile

        Session
        Subject
        Task
        Protocol
        Treatment
        Dose
        Label
        Experimenter
        TargetFP
        
        DLCTracking
        BehTable

        TraceMatrix
        TraceInterp
        TraceMedian

        DistMatLw
        DistMatDtw
    end

    properties (Dependent)
        SaveName

        NumTrials
        Trial
        TrialInfo

        TimeFromIn
        TimeFromOut
        TimeFromCue
        TimeWarpHD

        PortLeft
        PortRight
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
    end

    methods
        function obj = GPSTrajSessionClass(DLCTrackingOutFile, BehClassFile)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.SessionFolder = fileparts(DLCTrackingOutFile);

            obj.DLCTrackingOutFile = DLCTrackingOutFile;
            obj.BehClassFile = BehClassFile;
            BehClass = load(BehClassFile, "obj");

            obj.Session      = BehClass.obj.Session;
            obj.Subject      = BehClass.obj.Subject;

            obj.Task         = BehClass.obj.Task;
            obj.Protocol     = BehClass.obj.Protocol;
            obj.TargetFP     = BehClass.obj.TargetFP;

            obj.Treatment    = BehClass.obj.Treatment;
            obj.Dose         = BehClass.obj.Dose;
            obj.Label        = BehClass.obj.Label;
            obj.Experimenter = BehClass.obj.Experimenter;

            obj.BehTable     = BehClass.obj.BehavTable;

            load(DLCTrackingOutFile, "DLCTrackingOut");
            obj.DLCTracking = DLCTrackingOut;

            fprintf("\n%s\n", obj.Session);

            obj.removeOddTrials();

            obj.TraceMatrix.In  = obj.get_matrix(obj.Features, obj.TimeFromIn, obj.TimeMatIn);   % time aligned to cent-poke-in time
            obj.TraceMatrix.Out = obj.get_matrix(obj.Features, obj.TimeFromOut, obj.TimeMatOut); % time aligned to cent-poke-out time
%             obj.TraceMatrix.Cue = obj.get_matrix(obj.Features, obj.TimeFromCue, obj.TimeMatCue); % time aligned to trigger-cue time
            obj.TraceMatrix.HD  = obj.get_matrix(obj.Features, obj.TimeWarpHD, obj.TimeMatHD); % time linear warped between cent-poke-in and cent-poke-out
        end

        %%
        function removeOddTrials(obj)

            % check time mapping error
            time_error = cellfun(@(t) sum(diff(t)<=0), obj.TimeFromIn);
            ind_odd = find(time_error>0);
            if ~isempty(ind_odd)
                fprintf("Remove %d trials for time mapping error\n", length(ind_odd));
            end
            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];
            
            % check odd cent-in pose
            angle_in = cellfun(@(a, t) a(t==0), obj.AngleHead, obj.TimeFromIn);
            ind_odd  = find(angle_in > mean(angle_in)+5*std(angle_in) | angle_in < mean(angle_in)-5*std(angle_in));
            if ~isempty(ind_odd)
                fprintf("Remove %d trials for wrong cent-in position\n", length(ind_odd));
            end
            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];

            % check odd cent-out location
            loc_out_pre = cellfun(@(a, t) a(find(t<=0, 5, 'last')), obj.PosXHead, obj.TimeFromOut, 'UniformOutput', false);
            loc_out_pre = cell2mat(loc_out_pre);
            loc_out_pre_m = median(loc_out_pre, 1, "omitnan");
            ind_odd = find(any(abs(loc_out_pre - loc_out_pre_m)>=50, 2));
%             ind_odd = find(loc_out_post>median(loc_out_pre)+50 | loc_out_post<median(loc_out_pre)-50);
            if ~isempty(ind_odd)
                fprintf("Remove %d trials for wrong cent-out location\n", length(ind_odd));
            end
            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];
        end

        %% 
        function save_name = get.SaveName(obj)
            save_name = sprintf("GPSTrajSessionClass_%s_%s_%s", obj.Protocol, upper(obj.Subject), obj.Session);
        end

        %%
        function trial_info = get.TrialInfo(obj)
            trial_info = obj.BehTable(obj.Trial, :);
        end

        function num_trials = get.NumTrials(obj)
            num_trials = length(obj.Trial);
        end

        function trial = get.Trial(obj)
            trial = obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:)';
        end

        %%
        function time_from_in = get.TimeFromIn(obj)
            time_from_in = cellfun(@(x) x(:, 4)', obj.DLCTracking.PoseTracking(1).PosData, 'UniformOutput', false)';
        end

        function time_from_out = get.TimeFromOut(obj)
            in2out = num2cell(1000 * (obj.TrialInfo.CentOutTime - obj.TrialInfo.CentInTime));
            time_from_out = cellfun(@(x, y) x - y, obj.TimeFromIn, in2out, 'UniformOutput', false);
        end

        function time_from_cue = get.TimeFromCue(obj)
            in2cue = zeros(obj.NumTrials, 1);
            
            ind_non_pre = obj.TrialInfo.Outcome~="Premature";
            in2cue(ind_non_pre) = 1000 * (obj.TrialInfo.TriggerCueTime(ind_non_pre) - obj.TrialInfo.CentInTime(ind_non_pre));
            ind_no_cue = obj.TrialInfo.Outcome=="Premature" | obj.TrialInfo.Cued==0;
            in2cue(ind_no_cue) = 1000 * obj.TrialInfo.FP(ind_no_cue);

            in2cue = num2cell(in2cue);
            time_from_cue = cellfun(@(x, y) x - y, obj.TimeFromIn, in2cue, 'UniformOutput', false);
        end

        function time_warped = get.TimeWarpHD(obj)
            in2out = num2cell(1000 * (obj.TrialInfo.CentOutTime - obj.TrialInfo.CentInTime));
            time_warped = cellfun(@(x, y) x ./ y, obj.TimeFromIn, in2out, 'UniformOutput', false);
        end

        function port_left = get.PortLeft(obj)
            port_left = arrayfun(@(x) (x.L - x.R) / 2, obj.DLCTracking.PortLoc', 'UniformOutput', false);
        end

        function port_right = get.PortRight(obj)
            port_right = arrayfun(@(x) (x.R - x.L) / 2, obj.DLCTracking.PortLoc', 'UniformOutput', false);
        end

        function port_vec = get.PortVec(obj)
            port_vec = arrayfun(@(x) x.R - x.L, obj.DLCTracking.PortLoc', 'UniformOutput', false);
        end

        function port_cent = get.PortCent(obj)
            port_cent = arrayfun(@(x) (x.R + x.L) / 2, obj.DLCTracking.PortLoc', 'UniformOutput', false);
        end

        %% Angle of head
        function head_angle = get.AngleHead(obj)

            indL = find(strcmp("ear_base_left" , obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_base_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData';
            posR = obj.DLCTracking.PoseTracking(indR).PosData';

            head_vec = cellfun(@(x, y) y(:, 1:2) - x(:, 1:2), posL, posR, 'UniformOutput', false);
            head_angle = cellfun(@(x, y) calAngle(x, y), head_vec, obj.PortVec, 'UniformOutput', false);

            angle_sign = cellfun(@(x, y) 2*(x(:, 2)>=y(2))-1, head_vec, obj.PortVec, 'UniformOutput', false);
            head_angle = cellfun(@(x, y) x.*y, head_angle, angle_sign, 'UniformOutput', false);
            head_angle = cellfun(@(x) smoothdata(x', "gaussian", 5), head_angle, 'UniformOutput', false);
        end

        %% Angle speed of head
        function ang_speed_head = get.AngSpeedHead(obj)
            
            d_a = cellfun(@(a) a([2:end end]) - a([1 1:end-1]), obj.AngleHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            ang_speed_head = cellfun(@(da, dt) da ./ dt, d_a, d_t, 'UniformOutput', false);
            ang_speed_head = cellfun(@(x) smoothdata(x, "gaussian", 5), ang_speed_head, 'UniformOutput', false);

        end

        %% Angle acceleration of head
        function ang_acc_head = get.AngAccHead(obj)
            
            d_a = cellfun(@(a) a([2:end end]) - a([1 1:end-1]), obj.AngleHead, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            ang_speed_head = cellfun(@(da, dt) da ./ dt, d_a, d_t, 'UniformOutput', false);

            d_d_v = cellfun(@(v) v([2:end end]) - v([1 1:end-1]), ang_speed_head, 'UniformOutput', false);
            d_t = cellfun(@(t) t([2:end end]) - t([1 1:end-1]), obj.TimeFromIn, 'UniformOutput', false);

            ang_acc_head = cellfun(@(ddv, dt) ddv ./ dt, d_d_v, d_t, 'UniformOutput', false);
            ang_acc_head = cellfun(@(x) smoothdata(x, "gaussian", 5), ang_acc_head, 'UniformOutput', false);

        end

        %% Position of head (X)
        function head_pos_x = get.PosXHead(obj)

            indL = find(strcmp("ear_base_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_base_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData';
            posR = obj.DLCTracking.PoseTracking(indR).PosData';

            head_pos_x = cellfun(@(x, y, z) (y(:, 1) + x(:, 1))/2 - z(1), posL, posR, obj.PortCent, 'UniformOutput', false);
            head_pos_x = cellfun(@(x) smoothdata(x', "gaussian", 5), head_pos_x, 'UniformOutput', false);
        end

        %% Position of head (Y)
        function head_pos_y = get.PosYHead(obj)

            indL = find(strcmp("ear_base_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_base_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData';
            posR = obj.DLCTracking.PoseTracking(indR).PosData';

            head_pos_y = cellfun(@(x, y, z) (y(:, 2) + x(:, 2))/2 - z(2), posL, posR, obj.PortCent, 'UniformOutput', false);
            head_pos_y = cellfun(@(x) smoothdata(x', "gaussian", 5), head_pos_y, 'UniformOutput', false);
        end

        %% Speed (X) of head
        function speed_x_head = get.SpeedXHead(obj)
            
            d_x = cellfun(@(x) obj.cal_diff(x), obj.PosXHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t), obj.TimeFromIn, 'UniformOutput', false);

            speed_x_head = cellfun(@(dx, dt) dx ./ dt, d_x, d_t, 'UniformOutput', false);
            speed_x_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_x_head, 'UniformOutput', false);
        end

        %% Speed (Y) of head
        function speed_y_head = get.SpeedYHead(obj)
            
            d_y = cellfun(@(y) obj.cal_diff(y), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t), obj.TimeFromIn, 'UniformOutput', false);

            speed_y_head = cellfun(@(dx, dt) dx ./ dt, d_y, d_t, 'UniformOutput', false);
            speed_y_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_y_head, 'UniformOutput', false);
        end

        %% Acceleration (X) of head
        function acc_x_head = get.AccXHead(obj)
            
            d_d_x = cellfun(@(x) obj.cal_diff(x', 2), obj.PosXHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);

            acc_x_head = cellfun(@(ddx, dt) ddx ./ (dt.^2), d_d_x, d_t, 'UniformOutput', false);
            acc_x_head = cellfun(@(x) smoothdata(x, "gaussian", 5), acc_x_head, 'UniformOutput', false);
        end

        %% Acceleration (Y) of head
        function acc_y_head = get.AccYHead(obj)
            
            d_d_y = cellfun(@(y) obj.cal_diff(y', 2), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);

            acc_y_head = cellfun(@(ddy, dt) ddy ./ (dt.^2), d_d_y, d_t, 'UniformOutput', false);
            acc_y_head = cellfun(@(x) smoothdata(x, "gaussian", 5), acc_y_head, 'UniformOutput', false);
        end

        %% Speed of head
        function speed_head = get.SpeedHead(obj)
            
            d_x = cellfun(@(x) obj.cal_diff(x'), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) obj.cal_diff(y'), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);

            speed_head = cellfun(@(dx, dy, dt) vecnorm([dx, dy], 2, 2) ./ dt, d_x, d_y, d_t, 'UniformOutput', false);
            speed_head = cellfun(@(x) smoothdata(x', "gaussian", 5), speed_head, 'UniformOutput', false);
        end

        %% Speed direction of head
        function speed_dir_head = get.SpeedDirHead(obj)
            
            d_x = cellfun(@(x) obj.cal_diff(x'), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) obj.cal_diff(y'), obj.PosYHead, 'UniformOutput', false);

            speed_vec = cellfun(@(dx, dy) [dx, dy], d_x, d_y, 'UniformOutput', false);
            port_vec_vert = cellfun(@(pv) [1 -1] .* [pv(2) pv(1)], obj.PortVec, 'UniformOutput', false);

            speed_dir_head = cellfun(@(x, y) calAngle(x, y), speed_vec, port_vec_vert, 'UniformOutput', false);

            angle_sign = cellfun(@(x, y) 2*(x(:, 1)>=y(1))-1, speed_vec, port_vec_vert, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x, y) x.*y, speed_dir_head, angle_sign, 'UniformOutput', false);
%             speed_dir_head = cellfun(@(x) unwrap(x, 180), speed_dir_head, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x) smoothdata(x', "gaussian", 5), speed_dir_head, 'UniformOutput', false);
        end

        %% Acceleration of head
        function acc_head = get.AccHead(obj)
            
            d_d_x = cellfun(@(x) obj.cal_diff(x', 2), obj.PosXHead, 'UniformOutput', false);
            d_d_y = cellfun(@(y) obj.cal_diff(y', 2), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);

            acc_head = cellfun(@(ddx, ddy, dt) vecnorm([ddx, ddy], 2, 2) ./ (dt.^2), d_d_x, d_d_y, d_t, 'UniformOutput', false);
            acc_head = cellfun(@(x) smoothdata(x', "gaussian", 5), acc_head, 'UniformOutput', false);
        
        end

        %% Change of speed direction (dPhi) of head
        function dphi_head = get.dPhiHead(obj)
            
            d_x = cellfun(@(x) obj.cal_diff(x'), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) obj.cal_diff(y'), obj.PosYHead, 'UniformOutput', false);

            speed_vec = cellfun(@(dx, dy) [dx, dy], d_x, d_y, 'UniformOutput', false);
            port_vec_vert = cellfun(@(pv) [1 -1] .* [pv(2) pv(1)], obj.PortVec, 'UniformOutput', false);

            speed_dir_head = cellfun(@(x, y) calAngle(x, y), speed_vec, port_vec_vert, 'UniformOutput', false);

            angle_sign = cellfun(@(x, y) 2*(x(:, 1)>=y(1))-1, speed_vec, port_vec_vert, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x, y) x.*y, speed_dir_head, angle_sign, 'UniformOutput', false);
            speed_dir_head = cellfun(@(x) unwrap(x, 180), speed_dir_head, 'UniformOutput', false);

            d_dir = cellfun(@(dir) obj.cal_diff(dir), speed_dir_head, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);
            dphi_head = cellfun(@(ddir, dt) ddir ./ (dt.^2), d_dir, d_t, 'UniformOutput', false);
            dphi_head = cellfun(@(x) smoothdata(x', "gaussian", 5), dphi_head, 'UniformOutput', false);

        end

        %% Save
        function save(obj, copy_dir)
            save_path = fullfile(obj.SessionFolder, obj.SaveName);
            save(save_path, 'obj');

            if nargin==2
                copy_path = fullfile(copy_dir, obj.SaveName);
                copyfile(save_path+".mat", copy_path+".mat");
            end
        end % save

        %% Plots
        function print(obj, copy_dir)
            save_path = fullfile(obj.SessionFolder, obj.SaveName);
            fig = obj.plot();
            exportgraphics(fig, save_path+".png", 'Resolution', 600);
            exportgraphics(fig, save_path+".pdf", 'ContentType', 'vector');
            saveas(fig, save_path, 'fig');

            if nargin==2
                if ~isfolder(copy_dir)
                    mkdir(copy_dir);
                end
                copy_path = fullfile(copy_dir, obj.SaveName);
                copyfile(save_path+".pdf", copy_path+".pdf");
                copyfile(save_path+".png", copy_path+".png");
                copyfile(save_path+".fig", copy_path+".fig");
            end
        end % print

        function fig = plot(obj)
            try
                set_matlab_default
            catch
                disp('You do not have "set_matlab_default"' )
            end

            fig = feval(obj.Task+".plotTrajSession", obj);
        end % plot

    end
end