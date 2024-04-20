classdef GPSPlot < handle

    properties
    end

    methods
        function obj = GPSPlot(); end

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
            defaultColor =  repmat({[0 0 0]}, 1, length(data));
            defaultLw = repmat(1.2, 1, length(data));
            defaultLs = repmat("-", 1, length(data));
            defaultFP = zeros(1, length(data));

            p = inputParser;
            addRequired(p,'ax', @ishandle);
            addRequired(p,'data', @iscell);
            addParameter(p,'FP', defaultFP, @(x) (isvector(x) && length(x)==length(data)));
            addParameter(p,'Color', defaultColor, @iscell);
            addParameter(p,'LineWidth', defaultLw, @isvector);
            addParameter(p,'LineStyle', defaultLs, @isvector);
            addParameter(p,'Reversed', 0);

            parse(p, ax, data, varargin{:});

            ax = p.Results.ax;
            data = p.Results.data;
            FP = p.Results.FP;
            Color = p.Results.Color;
            LineWidth = p.Results.LineWidth;
            LineStyle = p.Results.LineStyle;
            Reversed = p.Results.Reversed;

            if length(Color)==1
                Color = repmat(Color, 1, length(data));
            elseif length(Color)~=length(data)
                error("Unmatched length between color and data");
            end
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
                    if Reversed
                        yline(ax, FP(i), 'LineStyle', LineStyle(i), 'Color', [.5 .5 .5]);
                    else
                        xline(ax, FP(i), 'LineStyle', LineStyle(i), 'Color', [.5 .5 .5]);
                    end
                end
                if ~isempty(data{i}.ci)
                    if Reversed
                        fill(ax, [data{i}.ci(1,:) flip(data{i}.ci(2,:))], [data{i}.x flip(data{i}.x)], 'r', ...
                            'FaceColor', Color{i}, 'FaceAlpha', .4, 'EdgeColor', 'none');
                    else
                        fill(ax, [data{i}.x flip(data{i}.x)], [data{i}.ci(1,:) flip(data{i}.ci(2,:))], 'r', ...
                            'FaceColor', Color{i}, 'FaceAlpha', .4, 'EdgeColor', 'none');
                    end
                end
                if Reversed
                    plot(ax, data{i}.f, data{i}.x, 'Color', Color{i}, 'LineWidth', LineWidth(i), 'LineStyle', LineStyle(i));
                else
                    plot(ax, data{i}.x, data{i}.f, 'Color', Color{i}, 'LineWidth', LineWidth(i), 'LineStyle', LineStyle(i));
                end
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
        function [ax, pos, ax_sz] = assign_ax_to_fig(obj, fig, n_row, n_col, pos, ax_sz)
            ax = cell(n_row, n_col);
            [dist_w, dist_h] = obj.get_plot_dist(ax, pos, ax_sz);

            for i = 1:n_row
                for j = 1:n_col
                    x_now = pos(1) + dist_w * (j-1);
                    y_now = pos(2) + dist_h * (n_row-i);

                    ax{i,j} = axes(fig, "Units", "centimeters", "Position", [x_now y_now ax_sz], ...
                        'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');
                end
            end
        end % assign_ax_to_fig

        function data_assigned = assign_data_to_ax(~, ax, data_cell)
            sz_ax = size(ax);
            sz_data = size(data_cell);
            if length(sz_ax)>2
                error("Dimensions of axes is more than 2");
            end
            if length(sz_data)>2
                error("Dimensions of data cell is more than 2");
            end

            data_assigned = cell(sz_ax);
            if sz_ax(1)==1
                if numel(data_cell)==numel(ax)
                    for i = 1:sz_ax(2)
                        data_assigned{i} = data_cell(i);
                    end
                else
                    if sz_data(2)==sz_ax(2)
                        for i = 1:sz_ax(2)
                            data_assigned{i} = data_cell(:,i);
                        end
                    else
                        error("Size of axes (%d, %d) is not equal to size of data cell (%d, %d)", sz_ax(1), sz_ax(2), sz_data(1), sz_data(2));
                    end
                end
            elseif sz_ax(2)==1
                if numel(data_cell)==numel(ax)
                    for i = 1:sz_ax(1)
                        data_assigned{i} = data_cell(i);
                    end
                else
                    if sz_data(1)==sz_ax(1)
                        for i = 1:sz_ax(1)
                            data_assigned{i} = data_cell(i,:);
                        end
                    else
                        error("Size of axes (%d, %d) is not equal to size of data cell (%d, %d)", sz_ax(1), sz_ax(2), sz_data(1), sz_data(2));
                    end
                end
            else
                if numel(data_cell)~=numel(ax)
                    error("Size of axes (%d, %d) is not equal to size of data cell (%d, %d)", sz_ax(1), sz_ax(2), sz_data(1), sz_data(2));
                else
                    if sz_ax(1)~=sz_data(1) && sz_ax(1)~=sz_data(2)
                        data_cell = data_cell';
                    end
                end
                for i = 1:sz_ax(1)
                    for j = 1:sz_sz(2)
                        data_assigned{i, j} = data_cell(i, j);
                    end
                end
            end
        end % assign_data_to_ax

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
