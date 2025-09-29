classdef GPSTrajClass < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        PerformanceType = ["Correct", "Premature", "Late", "Wrong", "Probe"];
        PerformanceCode = [1, -1, -2, 0, -3];
        Definition = [
            "ShuttleTime  (ST): Time of moving from port_init to port_cent";
            "HoldDuration (HD): Time from port_cent poke out to port_choice poke in (for hold paradigm only)";
            "ReactionTime (RT): For hold paradigm, it is time from tone to port_cent poke-out; for free paradigm, it is the time from tone to port_choice poke-in";
            "MovementTime (MT): For hold paradigm, it is the time from port_cent poke-out to port_choice poke-in; for free paradigm, leave it empy; for autoshaping paradigm, it is the time from port_cent poke-in to port_choice poke-in";
            "ChoiceTime   (CT):  Time from tone to port_choice poke in (for hold paradigm only)";
            ];
        Ports      = ["L", "R"];
        LeftRight  = [ 1 ,  2 ]; % Correct port, 1 for left; 2 for right
        CueUncue   = [ 1 ,  0 ];
        Guidance   = [ 1 ,  0,  -1]; % Guided, Empty, Filled
        BandWidth  = 0.05;

        % Time points for interpolation
        TimeMatIn   = -100:20:3000;
        TimeMatOut  = -2100:20:1000;
        TimeMatCue  = -2000:20:1100;
        TimeMatHD   = -0.25:0.01:1.25;

        % Movement features extracted
        Features = [
            "AngleHead";
            "AngSpeedHead";
            "AngAccHead";
            "PosXHead";
            "PosYHead";
            "SpeedXHead";
            "SpeedYHead";
            "AccXHead";
            "AccYHead";
            "SpeedHead";
            "SpeedDirHead";
            "AccHead";
            "dPhiHead"
            ];

        % Features for distance calculating
        FeatureDist = [
            "AngleHead";
            "PosXHead";
            "PosYHead";
            ];
    end

    methods
        function obj = GPSTrajClass(); end

        %% Date pre-processing
        function data_sorted = sort_data(~, data_origin, sort_refs, sort_codes)
            b_class = GPSBehClass();
            data_sorted = b_class.sort_data(data_origin, sort_refs, sort_codes);
        end % sort_data

        function table_sorted = sort_table(~, table_raw, sort_vars, sort_codes)
            b_class = GPSBehClass();
            table_sorted = b_class.sort_table(table_raw, sort_vars, sort_codes);
        end % sort_table

        function data_spliced = splice_data(~, data_cell)
            b_class = GPSBehClass();
            data_spliced = b_class.splice_data(data_cell);
        end % sort_table

        %% Calculate diff
        function d_x = cal_diff(obj, x, n)
            if nargin==2
                n = 1;
            end
            d_x = x([2:end end]) - x([1 1:end-1]);
            d_x(2:end-1) = d_x(2:end-1) ./ 2;
            if n > 1
                d_x = obj.cal_diff(d_x, n-1);
            end
        end % cal_diff

        %% Form transformation 
        function [M, t_M] = trace2mat(~, trace, time_trace, time_matrix)
            
            num_timepoint = length(time_matrix);

            M = cellfun(@(trace, t) interp1(t, trace', time_matrix, "linear"), trace, time_trace, 'UniformOutput', false);
            M = cellfun(@(m) reshape(m, 1, num_timepoint, []), M, 'UniformOutput', false);
            M = cell2mat(M);

            t_M = time_matrix;
        end % trace2mat

        function [T, t_T] = mat2trace(~, T, time_matrix, time_trace)

            [num_trace, num_timepoint, num_feature] = size(T, [1 2 3]);

            if nargin<4
                time_trace = repmat(mat2cell(time_matrix, 1, length(time_matrix)), num_trace, 1);
            end

            T = mat2cell(T, ones(num_trace, 1));
            T = cellfun(@(m) reshape(m, num_timepoint, num_feature), T, 'UniformOutput', false);
            T = cellfun(@(m, t) interp1(time_matrix', m, t', "linear"), T, time_trace, 'UniformOutput', false);
            T = cellfun(@(t) t', T, 'UniformOutput', false);

            t_T = time_trace;
        end % mat2trace

        function A = trace2array(~, trace)
            A = cellfun(@(x) x', trace, 'UniformOutput', false);
            A = cell2mat(A);
        end % trace2array

        function T = array2trace(~, array, time_trace)
            T = mat2cell(array, cellfun(@(x) length(x), time_trace));
            T = cellfun(@(x) x', T, 'UniformOutput', false);
        end % array2trace

        %% Trace processiong
        function trace_normalized = normalize_trace(~, trace, norm_method)
            
            if nargin<3
                norm_method = 'range';
            end
            trace_all = cell2mat(trace');
            trace_all_normalized = normalize(trace_all, 2, norm_method);

            trace_normalized = mat2cell(trace_all_normalized, size(trace_all_normalized, 1), cellfun(@(x) size(x, 2), trace));
            trace_normalized = trace_normalized';
        end % normalize_trace

        function trace_gathered = gather_trace(obj, body_parts)

            body_parts = string(body_parts);
            trace_gathered = obj.(body_parts(1));

            for i = 2:length(body_parts)
                trace_gathered = cellfun(@(x, y) [x; y], trace_gathered, obj.(body_parts(i)), 'UniformOutput', false);
            end
        end % gather_trace

        function [trace_interp, t_interp] = interp_trace(obj, trace, time_trace, time_matrix)
            mat = obj.trace2mat(trace, time_trace, time_matrix);
            trace_interp = obj.mat2trace(mat, time_matrix);
            t_interp = time_matrix;
        end % interpolate_trace

        function [trace_trim, time_trim] = trim_trace(~, trace, time_trace, varargin)

            % input
            P = inputParser;

            addRequired(P, 'trace', @iscell);
            addRequired(P, 'time_trace', @(x) iscell(x) && length(x)==length(trace));
            addOptional(P, 'bound_l', [], @isnumeric);
            addOptional(P, 'bound_u', [], @isnumeric);

            parse(P, trace, time_trace, varargin{:});

            trace = P.Results.trace;
            time_trace = P.Results.time_trace;
            bound_l = P.Results.bound_l;
            bound_u = P.Results.bound_u;

            %
            if isempty(bound_l)
                trace_trim = trace;
                time_trim  = time_trace;
            elseif length(bound_l)==1
                trace_trim = cellfun(@(x, t) x(:, t>=bound_l), trace, time_trace, 'UniformOutput', false);
                time_trim  = cellfun(@(t) t(:, t>=bound_l), time_trace, 'UniformOutput', false);
            elseif length(bound_l)==length(trace)
                trace_trim = cellfun(@(x, t, b) x(:, t>=b), trace, time_trace, num2cell(bound_l), 'UniformOutput', false);
                time_trim  = cellfun(@(t, b) t(:, t>=b), time_trace, num2cell(bound_l), 'UniformOutput', false);
            else
                error('Length of trace and bound_l should be matched');
            end

            %
            if isempty(bound_u)

            elseif length(bound_u)==1
                trace_trim = cellfun(@(x, t) x(:, t<=bound_u), trace_trim, time_trim, 'UniformOutput', false);
                time_trim  = cellfun(@(t) t(:, t<=bound_u), time_trim, 'UniformOutput', false);
            elseif length(bound_u)==length(trace)
                trace_trim = cellfun(@(x, t, b) x(:, t<=b), trace_trim, time_trim, num2cell(bound_u), 'UniformOutput', false);
                time_trim  = cellfun(@(t, b) t(:, t<=b), time_trim, num2cell(bound_u), 'UniformOutput', false);
            else
                error('Length of trace and bound_l should be matched');
            end
        end % trim_trace

        %% Map trace (x,y) location
        function [trace_mapped, map_sz, x_c, y_c] = map_trace(~, x, y, x_bins, y_bins, varargin)
            % input
            P = inputParser;

            addOptional(P, 'smooth', true, @(x) isnumeric(x) || islogical(x));
            addOptional(P, 'sigma', 0.5, @(x) isnumeric(x) && x>0);

            parse(P, varargin{:});

            smooth = P.Results.smooth; % whether to smooth the trajectory
            sigma  = P.Results.sigma;  % the sigma value of gaussian filter kernel

            % get trial number, the length of x-pos and y-pos should be equal
            n_trial = length(x);
            if length(y) ~= n_trial
                error("length of x and y do not match");
            end

            % initialize the map
            map_sz_x = length(x_bins) - 1;
            map_sz_y = length(y_bins) - 1;
            map_sz = [map_sz_y map_sz_x];

            trace_mapped = zeros(map_sz_y, map_sz_x, n_trial);

            % get bin centers
            x_c = (x_bins(1:end-1) + x_bins(2:end)) ./ 2;
            y_c = (y_bins(1:end-1) + y_bins(2:end)) ./ 2;

            % map trajectory for each trial
            for i = 1:n_trial
                % assign traj positions to discrete bins
                x_i = discretize(x{i}, x_bins);
                y_i = discretize(y{i}, y_bins);
                if length(x_i)~=length(y_i)
                    error("length of x and y in trial %d do not match", i);
                end

                % add visit count to the map
                for j = 1:length(x_i)
                    if ~isnan(y_i(j)) && ~isnan(x_i(j))
                        trace_mapped(y_i(j), x_i(j), i) = trace_mapped(y_i(j), x_i(j), i) + 1;
                    end
                end
            end

            % smooth the trajectory
            if smooth
                for i = 1:n_trial
                    % gaussian filter
                    trace_mapped(:, :, i) = imgaussfilt(trace_mapped(:, :, i), sigma);
                end
            end
        end % map_trace

        %% Align features to matrix
        function feature_mat = get_matrix(obj, features, time_trace, time_matrix, ind)
            if nargin < 5
                ind = 1:length(time_trace);
            end
            
            feature_mat = struct();
            for i = 1:length(features)
                feature_mat.(features(i)) = obj.trace2mat(obj.(features(i))(ind), time_trace(ind), time_matrix);
            end
        end % get_matrix

        %% Calculate distance matrix
        function dist_mat_lw = get_dist_lw(obj, features)
            dist_mat_lw = struct();
            trace_all = obj.gather_trace(features);
            periods = ["AP", "FP", "HD", "MT", "CT"];
            info = obj.TrialInfo;

            for i = 1:length(periods)
                period_i = periods(i);
                switch period_i
                    case "AP"
                        info_i = info;
                        trace = cellfun(@(x, t) x(:, t<0), trace_all, obj.TimeFromIn, 'UniformOutput', false);
                    case "FP"
                        ind = find(info.Outcome=="Correct");
                        info_i = info(ind, :);
                        trace = cellfun(@(x, t, fp) x(:, t>=0 & t<=fp), trace_all(ind), obj.TimeFromIn(ind), num2cell(info.FP(ind)*1000), 'UniformOutput', false);
                    case "HD"
                        info_i = info;
                        trace = cellfun(@(x, t) x(:, t>=0 & t<=1), trace_all, obj.TimeWarpHD, 'UniformOutput', false);
                    case "MT"
                        ind = find(info.Outcome=="Correct");
                        info_i = info(ind, :);
                        trace = cellfun(@(x, t, mt) x(:, t>=0 & t<=mt), trace_all(ind), obj.TimeFromOut(ind), num2cell(info.MT(ind)*1000), 'UniformOutput', false);
                    case "CT"
                        ind = find(info.Outcome=="Correct");
                        info_i = info(ind, :);
                        trace = cellfun(@(x, t, ct) x(:, t>=0 & t<=ct), trace_all(ind), obj.TimeFromCue(ind), num2cell(info.CT(ind)*1000), 'UniformOutput', false);
                end
                dist_mat_lw.(period_i) = struct();
                dist_mat_lw.(period_i).info = info_i;
                dist_mat_lw.(period_i).dist = obj.cal_dist(trace, 'warp_method', 'linear');
            end
        end % get_dist_lw

        %% Calculate trace median with 95% ci
        function trace_median = cal_trace_median(~, trace_matrix, time_matrix, t_bound, alpha)

            if nargin < 4 || isempty(t_bound)
                t_bound = [min(time_matrix) max(time_matrix)];
            end
            if nargin < 5
                alpha = 0.05;
            end

            if length(time_matrix) ~= size(trace_matrix)
                error('Length of time and trace dont match');
            end

            ind = time_matrix>=t_bound(1) & time_matrix<=t_bound(2);

            trace_median.time = time_matrix(ind);
            
            trace_matrix = trace_matrix(:, ind);
            trace_median.trace = median(trace_matrix, 1, 'omitnan');
            trace_median.trace = smoothdata(trace_median.trace, "gaussian", 5);

            num_valid = sum(~isnan(trace_matrix), 1);
            ind_valid = find(num_valid>5);
            trace_median.ci = nan(2, length(trace_median));
            trace_median.ci(:, ind_valid) = bootci(1000, {@(x) median(x, 'omitnan'), trace_matrix(:, ind_valid)}, 'type', 'cper', 'alpha', alpha);
            trace_median.ci = smoothdata(trace_median.ci, 2, "gaussian", 5);
        end % cal_trace_median

        %% Calculate distance matrix
        function dist_mat = cal_dist(obj, trace, varargin)
            % parsing input
            P = inputParser;

            addParameter(P, 'warp_method', 'linear', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'trace_norm_method', 'zscore', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'dist_norm_method', 'range', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'disp_iter', false, @islogical);

            parse(P, varargin{:});

            warp_method = P.Results.warp_method;
            trace_norm_method = P.Results.trace_norm_method;
            dist_norm_method = P.Results.dist_norm_method;
            disp_iter = P.Results.disp_iter;

            % warp trace
            if strcmp(warp_method, 'linear')
                t_origin = cellfun(@(x) linspace(0, 1, size(x, 2)), trace, 'UniformOutput', false);
                trace_len = cellfun(@(x) size(x, 2), trace);
                t_warped = linspace(0, 1, max(trace_len));
                trace = obj.interp_trace(trace, t_origin, t_warped);
            end

            % normalize trace
            switch trace_norm_method
                case 'none'
                otherwise
                    trace = obj.normalize_trace(trace, trace_norm_method);
            end

            % calculate distance
            num_trials = length(trace);
            dist_mat = zeros(num_trials);

            num_to_cal = .5*num_trials^2 - num_trials;
            num_calculated = 0;

            for m = 1:num_trials
                for n = m+1:num_trials
                    if any(isnan(trace{m}), "all") || any(isnan(trace{n}), "all")
                        dist_mat(m, n) = nan;
                    else
                        num_points = max([size(trace{m}, 2), size(trace{n}, 2)]);
                        switch warp_method
                            case 'dtw'
                                dist_mat(m, n) = dtw(trace{m}, trace{n}) / num_points;
                            case 'linear'
                                dist_mat(m, n) = sum(sqrt(sum((trace{m}-trace{n}).^2, 1))) / num_points;
                            case 'correlation'
                                corr_coef = corrcoef(trace{m}, trace{n});
                                dist_mat(m, n) = corr_coef(1,2);
                        end
                        num_calculated = num_calculated + 1;

                        if disp_iter
                            if ~mod(num_calculated, floor(num_to_cal/10))
                                fprintf('%.0f / %.0f    - %.0f%%\n', num_calculated, num_to_cal, 100*num_calculated/num_to_cal);
                            end
                        end
                    end
                end
            end
            dist_mat = dist_mat + dist_mat';
            if warp_method=="correlation"
                for m = 1:num_trials
                    dist_mat(m, m) = 1;
                end
            end

            % normalize distance matrix
            switch dist_norm_method
                case 'none'
                otherwise
                    dist_mat = normalize(dist_mat(:), dist_norm_method);
                    dist_mat = reshape(dist_mat, num_trials, num_trials);
            end
        end % cal_dist

        %% plot function
        function ax = plot_trace(~, ax, time_trace, trace, varargin)
            % check aquired input
            if length(time_trace)~=length(trace)
                error("Length of time_trace and trace should be the same.");
            end
            n_trace = length(trace);

            % parsing input
            P = inputParser;

            addParameter(P, 'color'  , [.6 .6 .6], @(x) (isnumeric(x) && size(x, 2)==3));
            addParameter(P, 'lw'     , 1.2       , @isnumeric);
            addParameter(P, 'ls'     , "-"       , @(x) all(ismember(x, ["-", "--", "-.", ":"])));
            addParameter(P, 'alpha'  , 0.3       , @isnumeric);
            addParameter(P, 'shuffle', true      , @islogical);

            parse(P, varargin{:});

            color = P.Results.color;
            lw = P.Results.lw;
            ls = P.Results.ls;
            alpha = P.Results.alpha;
            shuffle = P.Results.shuffle;

            % check optional input
            if size(color, 1)==1
                color = repmat(color, n_trace, 1);
            elseif length(color)~=n_trace
                error("Unmatched length between color and trace");
            end
            if length(lw)==1
                lw = repmat(lw, n_trace, 1);
            elseif length(lw)~=n_trace
                error("Unmatched length between lw and trace");
            end
            if length(ls)==1
                ls = repmat(ls, n_trace, 1);
            elseif length(ls)~=n_trace
                error("Unmatched length between ls and trace");
            end
            if length(alpha)==1
                alpha = repmat(alpha, n_trace, 1);
            elseif length(alpha)~=n_trace
                error("Unmatched length between alpha and trace");
            end

            % plot trace
            set(ax, 'NextPlot', 'add');
            for i = 1:n_trace
                plot(ax, time_trace{i}, trace{i}, 'Color', [color(i, :) alpha(i)], 'LineStyle', ls(i), 'LineWidth', lw(i));
            end
            if shuffle
                trace_plot = ax.Children;
                id_shuffle = randperm(length(trace_plot));
                ax.Children = trace_plot(id_shuffle);
            end
        end % plot_trace

        function ax = plot_trace_multi(obj, ax, time_cell, trace_cell, varargin)
            % check aquired input
            size_error = @(x, y) ndims(x)~=ndims(y) || any(size(x)~=size(y));
            if size_error(time_cell, trace_cell)
                error("Size of time_cell and trace_cell should be the same.");
            elseif ~ismatrix(trace_cell)
                error("Please dont put in more than 2 conditions.")
            end
            [n_1, n_2] = size(trace_cell, [1 2]);

            % default param value
            default_color   = repmat({[.6 .6 .6]}, [n_1, n_2]);
            default_lw      = repmat({1.2}, [n_1, n_2]);
            default_ls      = repmat({'-'}, [n_1, n_2]);
            default_alpha   = repmat({0.3}, [n_1, n_2]);

            % parsing input
            P = inputParser;

            addParameter(P, 'color'  , default_color);
            addParameter(P, 'lw'     , default_lw);
            addParameter(P, 'ls'     , default_ls);
            addParameter(P, 'alpha'  , default_alpha);
            addParameter(P, 'shuffle', true, @islogical);

            parse(P, varargin{:});

            color = P.Results.color;
            lw = P.Results.lw;
            ls = P.Results.ls;
            alpha = P.Results.alpha;
            shuffle = P.Results.shuffle;

            % check optional input
            if ~iscell(color)
                color = repmat({color}, [n_1, n_2]);
            elseif length(color)==1
                color = repmat(color, [n_1, n_2]);
            elseif size_error(color, trace_cell)
                error("Unmatched length between color and trace");
            end
            if ~iscell(lw)
                lw = repmat({lw}, [n_1, n_2]);
            elseif length(lw)==1
                lw = repmat(lw, [n_1, n_2]);
            elseif size_error(lw, trace_cell)
                error("Unmatched length between lw and trace");
            end
            if ~iscell(ls)
                ls = repmat({ls}, [n_1, n_2]);
            elseif length(ls)==1
                ls = repmat(ls, [n_1, n_2]);
            elseif size_error(ls, trace_cell)
                error("Unmatched length between ls and trace");
            end
            if ~iscell(alpha)
                alpha = repmat({alpha}, [n_1, n_2]);
            elseif length(alpha)==1
                alpha = repmat(alpha, [n_1, n_2]);
            elseif size_error(alpha, trace_cell)
                error("Unmatched length between alpha and trace");
            end

            % plot trace
            set(ax, 'NextPlot', 'add');
            for i = 1:n_1
                for j = 1:n_2
                    obj.plot_trace(ax, time_cell{i,j}, trace_cell{i,j}, 'color', color{i,j}, 'lw', lw{i,j}, 'ls', ls{i,j}, 'alpha', alpha{i,j}, 'shuffle', false);
                end
            end
            if shuffle
                ax.Children = ax.Children(randperm(length(ax.Children)));
            end
        end % plot_trace_multi

    end
end