classdef GPSTrajSessionClass < GPSTrajClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        SessionFolder

        DLCTrackingOutFile
        BehClassFile

        ANMInfoFile

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

        MatIn
        MatOut
        MatWarp

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

            obj = obj.removeOddTrials;

            obj.MatIn   = obj.get_matrix(obj.Features, obj.TimeFromIn, obj.TimeMatIn);
            obj.MatOut  = obj.get_matrix(obj.Features, obj.TimeFromOut, obj.TimeMatOut);
            obj.MatWarp = obj.get_matrix(obj.Features, obj.TimeWarped, obj.TimeMatWarp);
        end

        %%
        function obj = removeOddTrials(obj)

            time_error = cellfun(@(t) sum(diff(t)<=0), obj.TimeFromIn);
            ind_odd = find(time_error>0);
            
            for i = 1:length(obj.DLCTracking.PoseTracking)
                obj.DLCTracking.PoseTracking(i).PosData(ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).BpodEventIndex(:, ind_odd) = [];
                obj.DLCTracking.PoseTracking(i).Performance(ind_odd) = [];
            end

            obj.DLCTracking.PortLoc(ind_odd) = [];
            
            angle_in = cellfun(@(a, t) a(t==0), obj.AngleHead, obj.TimeFromIn);
            ind_odd  = find(angle_in > mean(angle_in)+5*std(angle_in) | angle_in < mean(angle_in)-5*std(angle_in));

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
            in2out = num2cell(1000 * (obj.BehTable.CentOutTime - obj.BehTable.CentInTime));
            in2out = in2out(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            time_from_out = cellfun(@(x, y) x - y, obj.TimeFromIn, in2out, 'UniformOutput', false);
        end

        function time_warped = get.TimeWarped(obj)
            in2out = num2cell(1000 * (obj.BehTable.CentOutTime - obj.BehTable.CentInTime));
            in2out = in2out(obj.DLCTracking.PoseTracking(1).BpodEventIndex(1,:));
            time_warped = cellfun(@(x, y) x ./ y, obj.TimeFromIn, in2out, 'UniformOutput', false);
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
            
            d_x = cellfun(@(x) obj.cal_diff(x'), obj.PosXHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);

            speed_x_head = cellfun(@(dx, dt) dx ./ dt, d_x, d_t, 'UniformOutput', false);
            speed_x_head = cellfun(@(x) smoothdata(x, "gaussian", 5), speed_x_head, 'UniformOutput', false);
        end

        %% Speed (Y) of head
        function speed_y_head = get.SpeedYHead(obj)
            
            d_y = cellfun(@(y) obj.cal_diff(y'), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromIn, 'UniformOutput', false);

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

        %% Align features to matrix
        function feature_mat = get_matrix(obj, features, time_trace, time_matrix)
            feature_mat = struct();
            for i = 1:length(features)
                feature_mat.(features(i)) = obj.trace2mat(obj.(features(i)), time_trace, time_matrix);
            end
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
            savename = fullfile(savepath, "GPSTrajectoryClass_" + obj.Task + "_" + string(Func) + "_" + upper(obj.Subject) + "_" + obj.Session);
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
                savename = fullfile(targetDir, "GPSTrajectoryClass_" + obj.Task + "_" + string(Func) + "_" + upper(obj.Subject) + "_" + obj.Session);
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
                'string', obj.Subject+" / "+obj.Session+" / "+obj.Task+" / "+char(obj.Treatment), 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

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
                'string', obj.Subject+" / "+obj.Session+" / "+obj.Task+" / "+char(obj.Treatment), 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

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