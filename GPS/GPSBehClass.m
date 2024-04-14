classdef GPSBehClass < handle
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
        PhaseCount = 50; % Trial number to determine the early or late training phase
        Ports      = ["L", "R"];
        LeftRight  = [ 1 ,  2 ]; % Correct port, 1 for left; 2 for right
        CueUncue   = [ 1 ,  0 ];
        Guidance   = [ 1 ,  0,  -1]; % Guided, Empty, Filled
        BandWidth  = .05;
    end

    properties (Dependent)
        Bins
    end

    methods
        %% Initiation
        function obj = GPSBehClass(); end
        
        %% Bin-edges for ksdensity functions
        function bins = get.Bins(~)
            bins.width = .002;

            bins.LogST = 0:bins.width:3;
            bins.HD    = 0:bins.width:3;
            bins.RT    = 0:bins.width:2;
            bins.MT    = 0:bins.width:2;
            bins.CT    = 0:bins.width:2.5;

            bins.widthInter     = 0.02;
            bins.Interruption   = 0:bins.widthInter:3;
        end

        %% get PDFs and CDFs
        function data_pdf = get_kde(obj, data, bin_edges, func, var_name, cal_ci)

            if cal_ci
                fprintf("... Calculate 95CI for %s of %s ... \n", func, var_name);
            end
            data_pdf = cellfun(@(x) obj.cal_kde(x, bin_edges, func, cal_ci), data, 'UniformOutput', false);

        end % get_kde

        function data_kde = cal_kde(obj, data_this, bin_edges, func, cal_ci)
            kde = @(x) ksdensity(x, bin_edges, 'Function', func, 'BandWidth', obj.BandWidth);

            data_kde.x = bin_edges;
            data_kde.ci = [];
            if length(data_this) > 5
                data_kde.f = kde(data_this);
                if cal_ci
                    data_kde.ci = bootci(1000, {kde, data_this}, 'type', 'cper', 'alpha', .05);
                end
            else
                data_kde.f = zeros(1, length(bin_edges));
            end
        end % calKDE

        %% Data pre-processing
        function data_sorted = sort_data(obj, data_origin, sort_refs, sort_codes)

            % Check inputs
            if ~isvector(data_origin)
                error("data_origin should be a vector");
            end
            data_origin = data_origin(:);

            if isempty(sort_refs)
                data_sorted = {data_origin};
                return
            end

            n_refs = cellfun(@length, sort_refs);
            n_data = length(data_origin);
            if ~all(n_refs==n_data)
                error("Length of reference data should be equal to the data to be sorted");
            end

            % Get output size
            sz = cellfun(@length, sort_codes);
            if length(sz)==1
                data_sorted = cell([sz, 1]);
            else
                data_sorted = cell(sz);
            end

            code_this = sort_codes{1};
            ref_this = sort_refs{1};
            if length(sz)>1 % recursion, if more than one reference variable remains
                for i = 1:length(code_this)
                    ind_this = ref_this==code_this(i);
                    data_now = data_origin(ind_this);
                    refs_now = cellfun(@(x) x(ind_this), sort_refs(2:end), 'UniformOutput', false);

                    data_sorted_i = obj.sort_data(data_now, refs_now, sort_codes(2:end));
                    data_sorted(i,:) = data_sorted_i(:);
                end
            else
                for i = 1:length(code_this)
                    ind_this = ref_this==code_this(i);
                    data_sorted{i} = data_origin(ind_this);
                end
            end
        end % sort_data

        function data_splited = split_data(obj, data_origin, n)
            data_splited = cell(n, 1);

            data_this = data_origin(obj.Stage==1);
            num_data = length(data_this);
            ind = 1:num_data;
%             [~, ~, indrmv] = rmoutliers_custome(data_origin);

            for i = 1:n
                ind_range = [i-1 i] * num_data/n;
                ind_this = ind>ind_range(1) & ind<=ind_range(2);
                data_splited{i} = data_this(ind_this);
                data_splited{i}(isnan(data_splited{i})) = [];
            end

        end % split_data

        function data_spliced = splice_data(obj, data_cell)
            sz = size(data_cell{1});
            if length(sz)==1
                data_spliced = cell([sz, 1]);
            else
                data_spliced = cell(sz);
            end

            if length(sz)==2 && sz(2)==1
                for i = 1:sz
                    data_this = cellfun(@(x) x{i}, data_cell, 'UniformOutput', false);
                    data_spliced{i} = vertcat(data_this{:});
                end
            else  % recursion
                sz_now = sz(2:end);
                if length(sz_now)==1
                    sz_now = [sz_now, 1];
                end
                for i = 1:sz(1)
                    data_now = cellfun(@(x) squeeze(x(i,:)), data_cell, 'UniformOutput', false);
                    data_now = cellfun(@(x) reshape(x, sz_now), data_now, 'UniformOutput', false);
                    data_splice_i = obj.splice_data(data_now);
                    data_spliced(i,:) = data_splice_i(:);
                end
            end
        end % splice_data

        %% Get statistics
        function stat_table = get_stat(obj, data_origin, sort_refs, sort_codes, sort_vars, cal_ci)

            data_sorted = obj.sort_data(data_origin, sort_refs, sort_codes);

            data_origin(isnan(data_origin)) = [];
            data_2575 = prctile(data_origin, [25, 75]);
            interq = data_2575(2) - data_2575(1);
            c = 5; % threshold for removing outliers
            bound_l = data_2575(1)-interq*c;
            bound_u = data_2575(2)+interq*c;
            
            stat_table_cell = cellfun(@(x) obj.cal_stat(x(x>=bound_l & x<=bound_u), cal_ci), data_sorted, 'UniformOutput', false);
            stat_table = vertcat(stat_table_cell{:});

            stat_table = obj.add_sort_info(stat_table, sort_codes, sort_vars);
        end % get_stat

        function stat_this = cal_stat(~, data_this, cal_ci)

            dur_this = calDur(data_this, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);

            N = sum(~isnan(data_this));

            Mean = dur_this.mean;
            STD = dur_this.std;
            SEM = dur_this.sem;

            Median = dur_this.median;
            Median_kde = dur_this.median_ksdensity;
            Q1 = dur_this.q1;
            Q3 = dur_this.q3;
            IQR = Q3 - Q1;

            if ~cal_ci
                stat_this = table(N, Mean, STD, SEM, Median, Median_kde, Q1, Q3, IQR);
            else
                if ~isnan(Median)
                    med_ci = bootci(1000, {@(x) median(x, 'omitnan'), data_this}, 'Type', 'cper');
                    iqr_ci = bootci(1000, {@iqr, data_this}, 'Type', 'cper');
                else
                    med_ci = [nan nan];
                    iqr_ci = [nan nan];
                end
                Median_ci_l = med_ci(1);
                Median_ci_u = med_ci(2);
                IQR_ci_l = iqr_ci(1);
                IQR_ci_u = iqr_ci(2);
                stat_this = table(N, Mean, STD, SEM, Median, Median_ci_l, Median_ci_u, Median_kde, Q1, Q3, IQR, IQR_ci_l, IQR_ci_u);
            end
            
        end % cal_stat

        %% Get performance ratio
        function perf_table = get_perfs(obj, data_outcome, sort_refs, sort_codes, sort_vars)

            data_sorted = obj.sort_data(data_outcome, sort_refs, sort_codes);

            perf_table_cell = cellfun(@(x) obj.cal_ratio(x, obj.PerformanceType), data_sorted, 'UniformOutput', false);
            perf_table = vertcat(perf_table_cell{:});

            perf_table = obj.add_sort_info(perf_table, sort_codes, sort_vars);
        end % get_perfs

        function ratio_table = cal_ratio(~, data_origin, labels)

            n_labels = length(labels);

            counts = zeros(n_labels, 1);
            for i = 1:n_labels
                counts(i) = sum(data_origin==labels(i));
            end
            counts_total = sum(counts);

            ratios = counts ./ counts_total;
            ratio_table = table('Size', [1, 1+n_labels], 'VariableTypes', repmat("double", 1, 1+n_labels), ...
                'VariableNames', ["N", labels]);
            ratio_table.N = counts_total;
            for i = 1:n_labels
                ratio_table.(labels(i)) = ratios(i);
            end
        end % cal_ratio

        %% Get tracks
        function perf_track = get_perf_track(obj, data_outcome, data_index, sort_refs, sort_codes)
            outcome_sorted = obj.sort_data(data_outcome, sort_refs, sort_codes);
            index_sorted = obj.sort_data(data_index, sort_refs, sort_codes);

            win_sz = cellfun(@(x) floor(length(x)/5), outcome_sorted, 'UniformOutput', false);
            step_sz = cellfun(@(x) max(1, floor(x/5)), win_sz, 'UniformOutput', false);

            outcome_track = cellfun(@(x, w_z, s_z) obj.cal_ratio_track(x, obj.PerformanceType, w_z, s_z), ...
                outcome_sorted, win_sz, step_sz, 'UniformOutput', false);
            index_track = cellfun(@(x, w_z, s_z) obj.cal_data_track(x, w_z, s_z), ...
                index_sorted, win_sz, step_sz, 'UniformOutput', false);

            perf_track = cellfun(@(o, i) addvars(o, i, 'NewVariableNames', "Pos", 'Before', "N"), ...
                outcome_track, index_track, 'UniformOutput', false);
        end % get_perf_track

        function ratio_track_table = cal_ratio_track(obj, data_origin, labels, win_sz, step_sz)
            n_labels = length(labels);
            num_trials = length(data_origin);
            length_track = floor((num_trials-win_sz+1)/step_sz);

            ratio_track_table = table('Size', [length_track, 1+n_labels], ...
                'VariableTypes', repmat("double", 1, 1+n_labels), ...
                'VariableNames', ["N", labels]);

            count_start = 1;
            i = 1;
            while count_start+win_sz-1 < num_trials
                win_this = count_start:(count_start+win_sz-1);
                ratio_track_table(i,:) = obj.cal_ratio(data_origin(win_this), labels);
                count_start = count_start + step_sz;
                i = i+1;
            end
        end % cal_ratio_track

        function data_track = cal_data_track(~, data_origin, win_sz, step_sz)
            num_trials = length(data_origin);
            length_track = floor((num_trials-win_sz+1)/step_sz);
            data_track = zeros(length_track, 1);

            count_start = 1;
            i = 1;
            while count_start+win_sz-1 < num_trials
                win_this = count_start:(count_start+win_sz-1);
                data_track(i) = sum(data_origin(win_this)) / win_sz;
                count_start = count_start + step_sz;
                i = i+1;
            end
        end % cal_data_track

        %% Add information to table
        function table_new = add_sort_info(~, table_raw, sort_codes, sort_vars)

            if isempty(sort_codes)
                table_new = table_raw;
                return;
            end

            sz = cellfun(@length, sort_codes);
            var_types = cellfun(@class, sort_codes, 'UniformOutput', false);
            info = table('Size', [prod(sz), length(sz)], ...
                'VariableTypes', string(var_types), ...
                'VariableNames', sort_vars);

            for i = 1:length(sz)
                info_this = repmat(sort_codes{i}, prod(sz(1:i-1)), prod(sz(i+1:end)));
                info.(sort_vars(i)) = info_this(:);
            end

            table_new = horzcat(info, table_raw);
        end % add_sort_info

        %% Plot functions
        function plot_trials(~, ax, data, varargin)

            defaultColor = [0 0 0];
            defaultSize = 18;
            defaultNumShow = 50;
            defaultFP = zeros(1, length(data));

            p = inputParser;
            addRequired(p,'ax', @ishandle);
            addRequired(p,'data', @iscell);
            addOptional(p,'n_show', defaultNumShow, @isnumeric);
            addParameter(p,'FP', defaultFP, @(x) (isvector(x) && length(x)==length(data)));
            addParameter(p,'FaceColor', defaultColor);
            addParameter(p,'MarkerColor', defaultColor);
            addParameter(p,'MarkerSize', defaultSize);

            parse(p, ax, data, varargin{:});

            ax = p.Results.ax;
            data = p.Results.data;
            n_show = p.Results.n_show;
            FP = p.Results.FP;
            FaceColor = p.Results.FaceColor;
            MarkerColor = p.Results.MarkerColor;
            MarkerSize = p.Results.MarkerSize;

            n_show = min([n_show min(cellfun(@(x) length(x), data))]);

            for i=1:length(data)
                if FP(i)~=0
                    fill(ax, [0 FP(i) FP(i) 0], [0 0 1 1]+i-1, 'r', ...
                        'FaceColor', FaceColor, 'FaceAlpha', .4*FP(i)/max(FP), 'EdgeColor', 'none');
                end
                id_show = randperm(length(data{i}), n_show);
                scatter(ax, data{i}(id_show), linspace(.1, .9, n_show)+i-1, MarkerSize, "filled", ...
                    "MarkerEdgeColor", "none", "MarkerFaceColor", MarkerColor, "MarkerFaceAlpha", .75*sqrt(FP(i)/max(FP)));
            end
            set(ax, 'XColor', 'none', 'YColor', 'none');
        end % plot_trials

        %
        function plot_distr(~, ax, data, varargin)
            defaultColor = [0 0 0];
            defaultLw = repmat(1.2, 1, length(data));
            defaultLs = repmat("-", 1, length(data));
            defaultFP = zeros(1, length(data));

            p = inputParser;
            addRequired(p,'ax', @ishandle);
            addRequired(p,'data', @iscell);
            addParameter(p,'FP', defaultFP, @(x) (isvector(x) && length(x)==length(data)));
            addParameter(p,'Color', defaultColor);
            addParameter(p,'LineWidth', defaultLw, @isvector);
            addParameter(p,'LineStyle', defaultLs, @isvector);

            parse(p, ax, data, varargin{:});

            ax = p.Results.ax;
            data = p.Results.data;
            FP = p.Results.FP;
            Color = p.Results.Color;
            LineWidth = p.Results.LineWidth;
            LineStyle = p.Results.LineStyle;

            if length(LineWidth)==1
                LineWidth = repmat(LineWidth, 1, length(data));
            elseif length(LineWidth)~=length(data)
                error("Unmatched length between lw and data");
            end
            if length(LineStyle)==1
                LineStyle = repmat(LineStyle, 1, length(data));
            elseif length(LineStyle)~=length(data)
                error("Unmatched length between ls and data");
            end

            for i=1:length(data)
                if FP(i)~=0
                    xline(ax, FP(i), 'LineStyle', LineStyle(i), 'Color', [.5 .5 .5]);
                end
                if ~isempty(data{i}.ci)
                    fill(ax, [data{i}.x flip(data{i}.x)], [data{i}.ci(1,:) flip(data{i}.ci(2,:))], 'r', ...
                        'FaceColor', Color, 'FaceAlpha', .4, 'EdgeColor', 'none');
                end
                plot(ax, data{i}.x, data{i}.f, 'Color', Color, 'LineWidth', LineWidth(i), 'LineStyle', LineStyle(i));
            end
        end % plot_distr

        %
        function plot_violin_compare(~, ax, datain_1, datain_2, band_width, varargin)
            defaultCate = string(1:length(datain_1));
            defaultColor = repmat({[0 0 0]}, 1, 2);
            defaultScale = 1;

            p = inputParser;
            addRequired(p,'ax', @ishandle);
            addRequired(p,'datain_1', @iscell);
            addRequired(p,'datain_2', @iscell);
            addRequired(p,'band_width', @isnumeric);
            addParameter(p,'cate_name', defaultCate);
            addParameter(p,'Color', defaultColor);
            addParameter(p,'Scale', defaultScale);
            addParameter(p,'TextMedian', false);

            parse(p, ax, datain_1, datain_2, band_width, varargin{:});

            ax = p.Results.ax;
            datain_1 = p.Results.datain_1;
            datain_2 = p.Results.datain_2;
            band_width = p.Results.band_width;
            cate_name = p.Results.cate_name;
            Color = p.Results.Color;
            Scale = p.Results.Scale;
            TextMedian = p.Results.TextMedian;

            n_data_1 = max(cellfun(@length, datain_1));
            n_data_2 = max(cellfun(@length, datain_2));
            n_data = max([n_data_1 n_data_2]);

            if length(n_data_1)~=length(n_data_2)
                error("input data dont match");
            else
                n_cate = length(datain_1);
            end

            %
            data_1 = nan(n_data, n_cate);
            data_2 = nan(n_data, n_cate);

            for i = 1:n_cate
                data_1(1:length(datain_1{i}), i) = Scale*datain_1{i};
                data_2(1:length(datain_2{i}), i) = Scale*datain_2{i};
            end

            data_all = [data_1 data_2];
            [~, ~, ind_rmv] = rmoutliers_custome(data_all(:));
            data_all(ind_rmv) = nan;
            data_1 = data_all(:, 1:n_cate);
            data_2 = data_all(:, n_cate+1:end);

            data_med_ci_1 = bootci(1000, @(x) median(x, 'omitnan'), data_1);
            data_med_ci_2 = bootci(1000, @(x) median(x, 'omitnan'), data_2);

            data_med_1 = median(data_1, 'omitnan');
            data_med_2 = median(data_2, 'omitnan');

            %
            axes(ax);
            v = violinplot({data_2, data_1}, cate_name, ...
                'ViolinColor', Color([2 1]), 'ViolinAlpha', {.1, .1}, ...
                'BandWidth', Scale*band_width, 'MarkerSize', 2, ...
                'ShowMedian', false, 'ShowWhiskers', false, 'ShowBox', false);
            xlim(ax, [.5 n_cate+.5]);

            for i = 1:n_cate
                v(i).ScatterPlot.MarkerFaceAlpha = .2;
                v(i).ScatterPlot2.MarkerFaceAlpha = .2;

                scatter(ax, i-.05, data_med_1(i), 16, Color{1}, "filled", ...
                    'LineWidth', 1, 'MarkerEdgeColor', 'flat', ...
                    'MarkerFaceColor', 'none');
                scatter(ax, i+.05, data_med_2(i), 16, Color{2}, "filled", ...
                    'LineWidth', 1, 'MarkerEdgeColor', 'flat', ...
                    'MarkerFaceColor', 'none');
                plot(ax, [i-.05 i-.05], data_med_ci_1(:,i), 'Color', Color{1}, 'LineWidth', 1);
                plot(ax, [i+.05 i+.05], data_med_ci_2(:,i), 'Color', Color{2}, 'LineWidth', 1);

                if TextMedian
                    switch Scale
                        case 1
                            median_format = '%.2f';
                        case 1000
                            median_format = '%.0f';
                    end
                    text(ax, i-.05, ax.YLim(1) + 0.6*range(ax.YLim), sprintf(median_format, data_med_1(i)), 'FontSize', 5, 'HorizontalAlignment', 'right');
                    text(ax, i+.05, ax.YLim(1) + 0.6*range(ax.YLim), sprintf(median_format, data_med_2(i)), 'FontSize', 5, 'HorizontalAlignment', 'left');
                end
                set(ax, 'Box', 'off');
            end
        end % plot_violin_compare

        %% For easy plot
        function [ax, pos, ax_sz] = assign_axes(~, fig, pos, ax_sz, n_row, n_col)
            ax = cell(n_row, n_col);
            [dist_w, dist_h] = get_plot_dist(ax, pos, ax_sz);

            for i = 1:n_row
                for j = 1:n_col
                    x_now = pos(1) + dist_w * (j-1);
                    y_now = pos(2) + dist_h * (n_row-i);

                    ax{i,j} = axes(fig, "Units", "centimeters", "Position", [x_now y_now ax_sz], ...
                        'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');
                end
            end
        end % assign_axes

        function add_serial(~, fig, ax_pos, serial_text, dist)
            pos = [ax_pos(1)+dist(1), ax_pos(2)+ax_pos(4)+dist(2), 1 1];
            ax = axes(fig, "Units", "centimeters", "Position", pos, ...
                'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out', ...
                'Color', 'none', 'XColor', 'none', 'YColor', 'none',...
                'XLim', [0 1], 'YLim', [0 1]); %
            text(ax, 0, 0, serial_text, 'FontSize', 9, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            drawnow;
        end % add_serial

        function [dist_w, dist_h] = get_plot_dist(~, ax_cell, pos, ax_sz)
            [n_h, n_w] = size(ax_cell);
            if n_w>1
                sep_w = ( pos(3) - n_w*ax_sz(1) ) / (n_w-1);
            else
                sep_w = 0;
            end
            if n_h > 1
                sep_h = ( pos(4) - n_h*ax_sz(2) ) / (n_h-1);
            else
                sep_h = 0;
            end

            dist_w = sep_w + ax_sz(1);
            dist_h = sep_h + ax_sz(2);
        end % get_plot_dist
    end
end