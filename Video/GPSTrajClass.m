classdef GPSTrajClass < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        PerformanceType = ["Correct", "Premature", "Late", "Wrong"];
        PerformanceCode = [1, -1, -2, 0];
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
        TimeMatIn   = -100:10:3000;
        TimeMatOut  = -2100:10:1000;
        TimeMatWarp = -0.25:0.005:1.25;

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

        %% Functions
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
        end

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

        function trace_interp = interp_trace(obj, trace, time_trace, time_matrix)
            mat = obj.trace2mat(trace, time_trace, time_matrix);
            trace_interp = obj.mat2trace(mat, time_matrix);
        end % interpolate_trace

        %% Calculate distance matrix
        function dist_mat = cal_dist(obj, trace, varargin)
            % parsing input
            P = inputParser;

            addParameter(P, 'warp_method', 'linear', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'trace_norm_method', 'range', @(x) (ischar(x) || isstring(x)));
            addParameter(P, 'dist_norm_method', 'range', @(x) (ischar(x) || isstring(x)));

            parse(P, varargin{:});

            warp_method = P.Results.warp_method;
            trace_norm_method = P.Results.trace_norm_method;
            dist_norm_method = P.Results.dist_norm_method;

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
                    num_points = max([length(trace{m}), length(trace{n})]);
                    switch warp_method
                        case 'dtw'
                            dist_mat(m, n) = dtw(trace{m}', trace{n}') / num_points;
                        case 'linear'
                            dist_mat(m, n) = sum(sqrt(sum((trace{m}-trace{n}).^2, 2))) / num_points;
                    end
                    num_calculated = num_calculated + 1;

                    if ~mod(num_calculated, floor(num_to_cal/100))
                        fprintf('%.0f / %.0f    - %.0f%%\n', num_calculated, num_to_cal, 100*num_calculated/num_to_cal);
                    end
                end
            end
            dist_mat = dist_mat + dist_mat';

            % normalize distance matrix
            switch dist_norm_method
                case 'none'
                otherwise
                    dist_mat = normalize(dist_mat(:), dist_norm_method);
                    dist_mat = reshape(dist_mat, num_trials, num_trials);
            end
        end % calDist

    end
end