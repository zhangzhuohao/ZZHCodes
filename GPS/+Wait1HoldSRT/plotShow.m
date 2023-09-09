function fig = plotShow(obj)

%%
opts.color = GPSColor(); % Color class for GPS
opts.mk = ["o", "x"]; % Scatter marker for correct and others
opts.ls = [":", "-.", "-"]; % Line style for [ShortFP, MedFP, LongFP]
opts.lw = [0.5, 1, 1.5]; % Line width for [ShortFP, MedFP, LongFP]

opts.session_date = char(obj.Sessions);
opts.session_date = opts.session_date(:, 5:8);

opts.plotsize = [8    4;
    8.5  3.5;
    8    2;
    8.5  5;
    5    5;
    4    4];

opts.sep_col = 1.2;
opts.sep_row = 1.5;

%%
fig = figure(24); clf(24);
set(gcf, 'unit', 'centimeters', 'position', [2 2 23 13], 'paperpositionmode', 'auto', 'color', 'w');

uicontrol('Style', 'text', 'parent', 24, 'units', 'normalized', 'position', [0.3 0.95 0.4 0.04],...
    'string', obj.Subject+" / "+obj.Task, 'fontsize', 11, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set axes and plot

% Performance track
ha1 = axes;
set(ha1, 'units', 'centimeters', 'position', [1.5 8, 7 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_performance_progress(ha1, obj, opts);

% legend
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [8.5 9, 1.5 3], ...
    'nextplot', 'add', 'fontsize', 7, 'tickdir', 'out');
line(ha11, [0 .6], [0 0], 'Color', opts.color.Correct, 'LineWidth', 1);
text(ha11, .7, 0, 'Correct', 'FontSize', 8, 'Color', opts.color.Correct);
line(ha11, [0 .6], [1 1], 'Color', opts.color.Premature, 'LineWidth', 1);
text(ha11, .7, 1, 'Premature', 'FontSize', 8, 'Color', opts.color.Premature);
line(ha11, [0 .6], [2 2], 'Color', opts.color.Wrong, 'LineWidth', 1);
text(ha11, .7, 2, 'Wrong', 'FontSize', 8, 'Color', opts.color.Wrong);
line(ha11, [0 .6], [3 3], 'Color', opts.color.Late, 'LineWidth', 1);
text(ha11, .7, 3, 'Late', 'FontSize', 8, 'Color', opts.color.Late);

set(ha11, 'xlim', [0 2], 'ylim', [0 7], 'xcolor', 'none', 'ycolor', 'none', 'ydir', 'reverse', 'color', 'none');

ha2 = axes;
set(ha2, 'units', 'centimeters', 'position', [1.5 1.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_reaction_time_violin(ha2, obj, opts, 'early');

ha21 = axes;
set(ha21, 'units', 'centimeters', 'position', [6.5 1.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_reaction_time_violin(ha21, obj, opts, 'late');
set(ha21, 'yticklabel', [], 'ylabel', []);

ha3 = axes;
set(ha3, 'units', 'centimeters', 'position', [13 2, 9 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_hold_duration_pdf_early_late(ha3, obj, "l", opts)

ha31 = axes;
set(ha31, 'units', 'centimeters', 'position', [13 6.5, 9 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_hold_duration_pdf_early_late(ha31, obj, "r", opts)
set(ha31, 'xticklabel', [], 'xlabel', []);

ha3.YLim(2) = max([ha3.YLim(2) ha31.YLim(2)]);
ha31.YLim(2) = ha3.YLim(2);

ha32 = axes;
set(ha32, 'units', 'centimeters', 'position', [19 10, 1.5 1.2], ...
    'nextplot', 'add', 'fontsize', 9, 'color', 'none', 'xcolor', 'none', 'ycolor', 'none', 'xlim', [0 3], 'ylim', [0 3], 'ticklength', [0.01 0.1], 'tickdir', 'out');
line(ha32, [1 2], [2 2], 'Color', opts.color.PhaseEarly, 'LineWidth', 1.5, 'LineStyle', '-');
line(ha32, [1 2], [1 1], 'Color', opts.color.PhaseLate, 'LineWidth', 1.5, 'LineStyle', '-');
text(ha32, 2.2, 2, "Early", 'FontSize', 9);
text(ha32, 2.2, 1, "Late", 'FontSize', 9);

%% ha1. Make a diagram of the setup
    function plot_performance_progress(ax, obj, opts)

        session_id = 1:obj.NumSessions;

        ind_this = find(obj.Performance.Foreperiod==0 & obj.Performance.TargetPort=="Both");

        plot(ax, session_id, obj.Performance.PrematureRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Premature, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w', 'MarkerSize', 6);
        plot(ax, session_id, obj.Performance.LateRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Late, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Late, 'markeredgecolor', 'w', 'MarkerSize', 6);
        plot(ax, session_id, obj.Performance.WrongRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Wrong, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Wrong, 'markeredgecolor', 'w', 'MarkerSize', 6);
        plot(ax, session_id, obj.Performance.CorrectRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Correct, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Correct, 'markeredgecolor', 'w', 'MarkerSize', 6);

        ax.XLabel.FontWeight = "Bold";
        ax.YLabel.FontWeight = "Bold";
        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Performance (%)';
        set(ax, 'xlim', [0.5 obj.NumSessions+0.5], 'xtick', session_id, 'tickdir', 'out');
    end

%% Reaction time violin
    function plot_reaction_time_violin(ax, obj, opts, phase)

        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--')

        num_violin = length(obj.Ports)*length(obj.MixedFP);
        thisRT = nan(obj.PhaseCount, num_violin);

        for fp_this = 1:length(obj.MixedFP)
            for p_this = 1:length(obj.Ports)
                switch phase
                    case {'early'}
                        ax.Title.String = "Early";
                        ax.Title.Color  = opts.color.PhaseEarly;
                        thisRT(:, 2*(fp_this-1)+p_this) = obj.RTSortedAll{fp_this, p_this}(1:obj.PhaseCount);
                    case {'late'}
                        ax.Title.String = "Late";
                        ax.Title.Color  = opts.color.PhaseLate;
                        thisRT(:, 2*(fp_this-1)+p_this) = obj.RTSortedAll{fp_this, p_this}(end-obj.PhaseCount+1:end);
                end
            end
        end
% 
%         violinplot({thisRT(:, 2:2:end), thisRT(:, 1:2:end)}, obj.Sessions, ...
%             'ViolinColor', {repmat(opts.color.PortR, num_violin/2, 1), repmat(opts.color.PortL, num_violin/2, 1)}, ...
%             'ScatterSize', 8, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false, 'Bandwidth', 0.1, 'ShowData', false);

        medians = median(thisRT);
        iqr_pos = prctile(thisRT, 75)-medians;
        iqr_neg = medians-prctile(thisRT, 25);
        
       
        errorbar(ax, (1:3)+0.1, medians(2:2:end), iqr_neg(2:2:end), iqr_pos(2:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
            'Color', opts.color.PortR, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.PortR, 'MarkerSize', 5);
        errorbar(ax, (1:3)-0.1, medians(1:2:end), iqr_neg(1:2:end), iqr_pos(1:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
            'Color', opts.color.PortL, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.PortL, 'MarkerSize', 5);

        ax.Title.FontWeight = "Bold";
        ax.XLabel.FontWeight = "Bold";
        ax.YLabel.FontWeight = "Bold";
        ax.XLabel.String = 'Foreperiod (s)';
        ax.YLabel.String = 'Reaction time (s)';
        set(ax, 'xlim', [.5 length(obj.MixedFP)+.5], 'ylim', [0 max(medians+iqr_pos)+0.02], 'xtick', 1:length(obj.MixedFP), ...
            'xticklabel', obj.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
    end

    function plot_hold_duration_pdf_early_late(ax, obj, port, opts)

        p_this = obj.Ports==upper(string(port));

        if all(all(cellfun(@(x) length(x)>=100, obj.HDSortedNone)))
            HDs = obj.HDSortedNone;
        else
            HDs = cellfun(@(x, y) [x; y], obj.HDSortedNone, obj.HDSortedSaline, 'UniformOutput', false);
        end

        if ~any(any(cellfun(@(x) length(x)>=100, HDs)))
            return
        end

        for fp_this = 1:length(obj.MixedFP)

            xline(ax, obj.MixedFP(fp_this), 'color', [.7 .7 .7], 'linewidth', opts.lw(fp_this), 'LineStyle', '-');

            hd_this = HDs{fp_this, p_this};

            if length(hd_this) < 2*obj.PhaseCount
                continue;
            end

            hd_early_pdf = ksdensity(hd_this(1:obj.PhaseCount), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_late_pdf  = ksdensity(hd_this(end-obj.PhaseCount+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');

            plot(ax, obj.Bins.HoldDuration, hd_early_pdf, ...
                'color', opts.color.PhaseEarly, 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
            plot(ax, obj.Bins.HoldDuration, hd_late_pdf, ...
                'color', opts.color.PhaseLate, 'linewidth', opts.lw(fp_this), 'LineStyle', '-');

        end

        switch upper(string(port))
            case {"L"}
                ax.YLabel.String = "Prob. density Left (1/s)";
                ax.YLabel.Color  = opts.color.PortL;
            case {"R"}
                ax.YLabel.String = "Prob. density Right (1/s)";
                ax.YLabel.Color  = opts.color.PortR;
        end

        ax.XLabel.FontWeight = "Bold";
        ax.YLabel.FontWeight = "Bold";
        ax.XLabel.String = "Hold duration (s)";
        set(ax, 'xlim', [0 2], 'ylimmode', 'auto', 'tickdir', 'out');
    end

end
