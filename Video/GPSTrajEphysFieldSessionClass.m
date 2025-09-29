classdef GPSTrajEphysFieldSessionClass < GPSTrajClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        SessionFolder

        DLCTracking
        r

        Experimenter
        Session
        Subject
        Task
        Protocol

        TargetFP
        BehTable

        TraceMatrix
    end

    properties (Dependent)
        SaveName

        NumTrials
        Trial
        TrialInfo

        TimeFramesE
        TimeFramesB

        BodyParts
        RefLocations

        EarLX
        EarLY

        EarRX
        EarRY

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
        function obj = GPSTrajFieldSessionClass(DLCTrackingOut, r)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.DLCTracking = DLCTrackingOut;
            BehClass = r.BehaviorClass;

            obj.Experimenter = BehClass.Experimenter;
            obj.Session      = BehClass.Session;
            obj.Subject      = BehClass.Subject;
            obj.Task         = BehClass.Task;
            obj.Protocol     = BehClass.Protocol;

            obj.TargetFP     = BehClass.TargetFP;
            obj.BehTable     = BehClass.BehavTable;

            obj.DLCTracking = DLCTrackingOut;

            fprintf("\n%s\n", obj.Session);

            obj.removeOddTrials();
        end

        %%
        function removeOddTrials(obj)

            % check time mapping error
            time_error = cellfun(@(t) sum(diff(t)<=0), obj.TimeFromCentIn);
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
        end

        %% 
        function save_name = get.SaveName(obj)
            save_name = sprintf("TrajField_%s_%s", obj.Subject, obj.Session);
        end

        %%
        function trial_info = get.TrialInfo(obj)
            trial_info = obj.BehTable(obj.Trial, :);
        end

        function num_trials = get.NumTrials(obj)
            num_trials = length(obj.Trial);
        end

        function trial = get.Trial(obj)
            trial = obj.DLCTracking.trial_id;
        end

        %%
        function t_frame_e = get.TimeFramesE(obj)
            t_frame_e = obj.DLCTracking.t_frame_e;
        end

        function t_frame_b = get.TimeFramesB(obj)
            t_frame_b = obj.DLCTracking.t_frame_b;
        end

        %%
        function body_parts = get.BodyParts(obj)
            body_parts = fields(obj.DLCTracking.BodyPart);
        end

        function reference_loc = get.RefLocations(obj)
            reference_loc = fields(obj.DLCTracking.RefLoc);
        end

        %% Position of head (X)
        function head_pos_x = get.PosXHead(obj)
            ear_l = obj.DLCTracking.BodyPart.EarLTop;
            ear_r = obj.DLCTracking.BodyPart.EarLTop;

            indL = find(strcmp("ear_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData';
            posR = obj.DLCTracking.PoseTracking(indR).PosData';

            head_pos_x = cellfun(@(x, y) (y(:, 1) + x(:, 1))/2, posL, posR, 'UniformOutput', false);
            head_pos_x = cellfun(@(x) smoothdata(x', "gaussian", 3), head_pos_x, 'UniformOutput', false);
        end

        %% Position of head (Y)
        function head_pos_y = get.PosYHead(obj)
            indL = find(strcmp("ear_left", obj.DLCTracking.BodyParts));
            indR = find(strcmp("ear_right", obj.DLCTracking.BodyParts));

            posL = obj.DLCTracking.PoseTracking(indL).PosData';
            posR = obj.DLCTracking.PoseTracking(indR).PosData';

            head_pos_y = cellfun(@(x, y) (y(:, 2) + x(:, 2))/2, posL, posR, 'UniformOutput', false);
            head_pos_y = cellfun(@(x) smoothdata(x', "gaussian", 3), head_pos_y, 'UniformOutput', false);
        end

        %% Speed (X) of head
        function speed_x_head = get.SpeedXHead(obj)
            
            d_x = cellfun(@(x) obj.cal_diff(x), obj.PosXHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t), obj.TimeFromCentIn, 'UniformOutput', false);

            speed_x_head = cellfun(@(dx, dt) dx ./ dt, d_x, d_t, 'UniformOutput', false);
            speed_x_head = cellfun(@(x) smoothdata(x, "gaussian", 3), speed_x_head, 'UniformOutput', false);
        end

        %% Speed (Y) of head
        function speed_y_head = get.SpeedYHead(obj)
            
            d_y = cellfun(@(y) obj.cal_diff(y), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t), obj.TimeFromCentIn, 'UniformOutput', false);

            speed_y_head = cellfun(@(dx, dt) dx ./ dt, d_y, d_t, 'UniformOutput', false);
            speed_y_head = cellfun(@(x) smoothdata(x, "gaussian", 3), speed_y_head, 'UniformOutput', false);
        end

        %% Acceleration (X) of head
        function acc_x_head = get.AccXHead(obj)
            
            d_d_x = cellfun(@(x) obj.cal_diff(x', 2), obj.PosXHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromCentIn, 'UniformOutput', false);

            acc_x_head = cellfun(@(ddx, dt) ddx ./ (dt.^2), d_d_x, d_t, 'UniformOutput', false);
            acc_x_head = cellfun(@(x) smoothdata(x, "gaussian", 3), acc_x_head, 'UniformOutput', false);
        end

        %% Acceleration (Y) of head
        function acc_y_head = get.AccYHead(obj)
            
            d_d_y = cellfun(@(y) obj.cal_diff(y', 2), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromCentIn, 'UniformOutput', false);

            acc_y_head = cellfun(@(ddy, dt) ddy ./ (dt.^2), d_d_y, d_t, 'UniformOutput', false);
            acc_y_head = cellfun(@(x) smoothdata(x, "gaussian", 3), acc_y_head, 'UniformOutput', false);
        end

        %% Speed of head
        function speed_head = get.SpeedHead(obj)
            
            d_x = cellfun(@(x) obj.cal_diff(x'), obj.PosXHead, 'UniformOutput', false);
            d_y = cellfun(@(y) obj.cal_diff(y'), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromCentIn, 'UniformOutput', false);

            speed_head = cellfun(@(dx, dy, dt) vecnorm([dx, dy], 2, 2) ./ dt, d_x, d_y, d_t, 'UniformOutput', false);
            speed_head = cellfun(@(x) smoothdata(x', "gaussian", 3), speed_head, 'UniformOutput', false);
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
            speed_dir_head = cellfun(@(x) smoothdata(x', "gaussian", 3), speed_dir_head, 'UniformOutput', false);
        end

        %% Acceleration of head
        function acc_head = get.AccHead(obj)
            
            d_d_x = cellfun(@(x) obj.cal_diff(x', 2), obj.PosXHead, 'UniformOutput', false);
            d_d_y = cellfun(@(y) obj.cal_diff(y', 2), obj.PosYHead, 'UniformOutput', false);
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromCentIn, 'UniformOutput', false);

            acc_head = cellfun(@(ddx, ddy, dt) vecnorm([ddx, ddy], 2, 2) ./ (dt.^2), d_d_x, d_d_y, d_t, 'UniformOutput', false);
            acc_head = cellfun(@(x) smoothdata(x', "gaussian", 3), acc_head, 'UniformOutput', false);
        
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
            d_t = cellfun(@(t) obj.cal_diff(t'), obj.TimeFromCentIn, 'UniformOutput', false);
            dphi_head = cellfun(@(ddir, dt) ddir ./ (dt.^2), d_dir, d_t, 'UniformOutput', false);
            dphi_head = cellfun(@(x) smoothdata(x', "gaussian", 3), dphi_head, 'UniformOutput', false);

        end

        %% Save
        function save(obj, dir_name)
            if nargin <= 1
                dir_name = './';
            end
            TrajField = obj;
            save(fullfile(dir_name, obj.SaveName), 'TrajField');
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

            fig = figure(11); clf(fig);
            set(fig, 'Visible', 'on', 'Units', 'centimeters', 'Position', [5 5 12.5 5.5], 'Color', 'w', 'toolbar', 'none');

%             fig_title = sprintf("%s / %s / %s / %s", obj.Subject, obj.Session, obj.Protocol, obj.Label);
%             set_fig_title(fig, fig_title);
% 
%             ax = axes(fig, 'Units', 'centimeters', 'Position', [.5 .5 10 4], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', ...
%                  'XLim', [-150 350], 'YLim', [-25 175], 'NextPlot', 'add', 'YDir', 'reverse');
%             line(ax, [-60 320], [0 0], 'LineWidth', 1, 'Color', 'k');
%             line(ax, [-60 320], [1 1] * obj.CorridorWidth, 'LineWidth', 1, 'Color', 'k');
%             line(ax, [1 1] * obj.CorridorLength, [-10 160], 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':');
%             line(ax, [0 0], [-10 160], 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':');
%             text(ax, 0, -20, 'Out', 'HorizontalAlignment', 'center');
%             text(ax, obj.CorridorLength, -20, 'In', 'HorizontalAlignment', 'center');
%             for i = 1:obj.NumTrials
%                 patch(ax, 'XData', [obj.PosHeadT{i}(1,:) nan], 'YData', [obj.PosHeadT{i}(2,:) nan], 'CData', [obj.TimeFromZoneIn{i} nan], 'LineWidth', 1, 'EdgeAlpha', .3, 'EdgeColor', 'flat');
%             end
%             colormap(ax, 'jet');
%             %             ax.CLim = [-.2 1.2];
%             ax.CLim = [-200 800];
% 
%             cb = colorbar(ax, 'Units', 'centimeters', 'Position', [10.8 .6 .3 3.8]);
%             cb.Label.String = 'Time from cross-in (ms)';
% %             cb.Ticks = [0 1];
% %             cb.TickLabels = {'In', 'Out'};
% %             cb.Label.String = 'norm. time';

        end % plot

    end
end