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
        TimeMatIn  = -100:10:3000;
        TimeMatOut = -2100:10:1000;
        TimeMatWarp = -0.25:0.005:1.25;
    end

    methods
        function obj = GPSTrajClass(); end

        %% Functions
        %% Form transformation 
        function [M, t_M] = trace2mat(~, trace, time_trace, time_matrix)
            
            num_timepoint = length(time_matrix);
            num_feature = size(trace{1}, 2);

            M = cellfun(@(trace, t) interp1(t, trace, time_matrix, "linear"), trace, time_trace, 'UniformOutput', false);
            M = M';
            M = cellfun(@(m) reshape(m, 1, num_timepoint, num_feature), M, 'UniformOutput', false);
            M = cell2mat(M);

            t_M = time_matrix;
        end % trace2mat

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
        end % mat2trace

        function A = trace2array(~, trace)
            A = cell2mat(trace');
        end % trace2array

        function T = array2trace(~, array, time_trace)
            T = mat2cell(array, cellfun(@(x) length(x), time_trace));
            T = T';
        end % array2trace

        %% Trace processiong
        function trace_normalized = normalize_trace(~, trace, norm_method)
            
            if nargin<3
                norm_method = 'range';
            end
            trace_all = cell2mat(trace');
            trace_all_normalized = normalize(trace_all, 1, norm_method);

            trace_normalized = mat2cell(trace_all_normalized, cellfun(@(x) size(x, 1), trace));
            trace_normalized = trace_normalized';
        end % normalizeTrace

        function trace_gathered = gather_trace(obj, body_parts)

            body_parts = string(body_parts);
            trace_gathered = obj.(body_parts(1));

            for i = 2:length(body_parts)
                trace_gathered = cellfun(@(x, y) [x y], trace_gathered, obj.(body_parts(i)), 'UniformOutput', false);
            end
        end % gatherTrace

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
                    switch warp_method
                        case 'dtw'
                            dist_mat(m, n) = dtw(trace{m}', trace{n}');
                        case 'linear'
                            dist_mat(m, n) = sum(sqrt(sum((trace{m}-trace{n}).^2, 2)));
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