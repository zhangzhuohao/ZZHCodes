function fig = plotSession(obj)

%% set the options
opts.color = GPSColor(); % Color class for GPS
opts.mk = ["o", "x"]; % Scatter marker for correct and others
opts.ls = ["-", ":"]; % Line style for [Cue, Uncue]
opts.lw = [0.5, 1, 1.5]; % Line width for [ShortFP, MedFP, LongFP]
opts.ticks = obj.TrialStartTime + cellfun(@(x)x(1), obj.CentPokeInTime); % Ticks for each trial

%%
fig = figure(22); clf(22);
set(gcf, 'unit', 'centimeters', 'position', [2 2 24.5 21.2], 'paperpositionmode', 'auto', 'color', 'w');

uicontrol('Style', 'text', 'parent', 22, 'units', 'normalized', 'position', [0.05 0.95 0.9 0.04],...
    'string', [obj.Subject ' / ' obj.Session ' / ' obj.Task ' / ' char(obj.Treatment)], 'fontsize', 10, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set the axes and plot
plot_size1 = [8 4];
plot_size2 = [2 4];
plot_size3 = [3.5 4];

% Maze diagram
ha1 = axes;
set(ha1, 'units', 'centimeters', 'position', [1.2 16, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_diagram(ha1, opts);

% Cumulative plot
ha2 = axes;
set(ha2, 'units', 'centimeters', 'position', [1.5 16-plot_size1(2)-1, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_cumulative(ha2, opts);

% Performance track
ha3 = axes;
set(ha3, 'units', 'centimeters', 'position', [1.5 16-2*plot_size1(2)-2, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_performance(ha3, obj, opts);

% Shuttle time
ha4 = axes;
set(ha4, 'units', 'centimeters', 'position', [1.5+plot_size1(1)+1 16, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_shuttle_time(ha4, obj, opts);

% PDF of shuttle time
ha5 = axes;
set(ha5, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8 16, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_shuttle_time_pdf(ha5, obj);

% Hold duration
ha6 = axes;
set(ha6, 'units', 'centimeters', 'position', [1.5+plot_size1(1)+1 16-plot_size1(2)-1, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration(ha6, obj, opts);

% PDF for hold duration
ha7 = axes;
set(ha7, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8 16-plot_size1(2)-1, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_pdf_l(ha7, obj, opts);

% CDF for hold duration
ha71 = axes;
set(ha71, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8+plot_size2(1)+0.8 16-plot_size1(2)-1, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_pdf_r(ha71, obj, opts);

% Reaction time
ha8 = axes;
set(ha8, 'units', 'centimeters', 'position', [1.5+plot_size1(1)+1 16-2*plot_size1(2)-2, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time(ha8, obj, opts);

% PDF for reaction time to PortL
ha9 = axes;
set(ha9, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8 16-2*plot_size1(2)-2, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_pdf_l(ha9, obj, opts);

% PDF for reaction time to PortR
ha91 = axes;
set(ha91, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8+plot_size2(1)+0.8 16-2*plot_size1(2)-2, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_pdf_r(ha91, obj, opts);

ha9.XLim(2) = max([ha9.XLim(2) ha91.XLim(2)]);
ha91.XLim   = ha9.XLim;

% Movement time
ha10 = axes;
set(ha10, 'units', 'centimeters', 'position', [1.5+plot_size1(1)+1 16-3*plot_size1(2)-3, plot_size1], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time(ha10, obj, opts);

% PDF for movement time to PortL
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8 16-3*plot_size1(2)-3, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time_pdf_l(ha11, obj, opts);

% PDF for movement time to PortR
ha111 = axes;
set(ha111, 'units', 'centimeters', 'position', [1.5+2*plot_size1(1)+1.8+plot_size2(1)+0.8 16-3*plot_size1(2)-3, plot_size2], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time_pdf_r(ha111, obj, opts);

ha11.XLim(2) = max([ha11.XLim(2) ha111.XLim(2)]);
ha111.XLim   = ha11.XLim;

% Collected performance
ha12 = axes;
set(ha12, 'units', 'centimeters', 'position', [1.5, 16-3*plot_size1(2)-3, plot_size3], 'nextplot', 'add', 'fontsize', 8);
plot_performance_collected(ha12, obj, opts);

% Legend for collected performance
ha121 = axes;
set(ha121, 'units', 'centimeters', 'position', [3.2, 16-3*plot_size1(2)-3+2*plot_size3(2)/3, plot_size3/3], 'nextplot', 'add', 'fontsize', 8);
legend_performance(ha121, opts);

% Violin plot for collected reaction time
ha13 = axes;
set(ha13, 'units', 'centimeters', 'position', [1.5+plot_size3(1)+1, 16-3*plot_size1(2)-3, plot_size3], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_collected(ha13, obj, opts);

%% ha1. Make a diagram of the setup
    function plot_diagram(ax, opts)

        line(ax, [1 1], [0 6], 'color', 'k', 'linewidth', 2);
        line(ax, [1 3], [0 0], 'color', 'k', 'linewidth', 2);
        line(ax, [3 5], [0 2], 'color', 'k', 'linewidth', 2);
        line(ax, [5 8], [2 2], 'color', 'k', 'linewidth', 2);
        line(ax, [1 3], [6 6], 'color', 'k', 'linewidth', 2);
        line(ax, [3 5], [6 4], 'color', 'k', 'linewidth', 2);
        line(ax, [5 8], [4 4], 'color', 'k', 'linewidth', 2);

        line(ax, [8 9], [4 4.8], 'color', 'k', 'linewidth', 2);
        line(ax, [8 9], [2 1.2], 'color', 'k', 'linewidth', 2);
        line(ax, [9 11], [4.8 4.8], 'color', 'k', 'linewidth', 2);
        line(ax, [9 11], [1.2 1.2], 'color', 'k', 'linewidth', 2);
        line(ax, [11 11], [1.2 4.8], 'color', 'k', 'linewidth', 2);

        viscircles(ax, [9.6, 3], 0.3, 'color', 'k', 'LineWidth', 1);
        text(ax, 8.6, 3, 'Init', 'FontWeight','bold', 'HorizontalAlignment', 'center');
        viscircles(ax, [3.5, 3], 0.3,  'color', 'k', 'LineWidth', 1);
        text(ax, 3.9, 3, 'Cent', 'FontWeight', 'bold');

        viscircles(ax, [3, 1.5], 0.3, 'color', opts.color.PortL);
        text(ax, 1.2, 1.5, 'Left', 'FontWeight','bold', 'HorizontalAlignment', 'left', 'Color', opts.color.PortL);
        viscircles(ax, [3, 4.5], 0.3,  'color', opts.color.PortR);
        text(ax, 1.2, 4.5, 'Right', 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Color', opts.color.PortR);

        set(ax, 'XColor', 'none', 'YColor', 'none', 'ylim', [0 6], 'xlim', [0 12]);
    end

%% ha2. Cumulative plot
    function plot_cumulative(ax, opts)

        stairs(opts.ticks, 1:length(opts.ticks), 'color', 'k', 'LineWidth', 1.5);

        ax.XLabel.String = 'Time in session (s)';
        ax.YLabel.String = 'Trial (count)';
        set(ax, 'XLim', [0 opts.ticks(end)+5], 'YLim', [0 length(opts.ticks)+5]);
    end

%% ha3. Plot performance over time
    function plot_performance(ax, obj, opts)

        plot(ax, obj.PerformanceTrackCue.WinPos, obj.PerformanceTrackCue.CorrectRatio,   'o', 'linestyle', '-', 'color', opts.color.Correct, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Correct,   'markeredgecolor', 'w');
        plot(ax, obj.PerformanceTrackCue.WinPos, obj.PerformanceTrackCue.WrongRatio,     'o', 'linestyle', '-', 'color', opts.color.Wrong, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Wrong,     'markeredgecolor', 'w');
        plot(ax, obj.PerformanceTrackCue.WinPos, obj.PerformanceTrackCue.PrematureRatio, 'o', 'linestyle', '-', 'color', opts.color.Premature, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
        plot(ax, obj.PerformanceTrackCue.WinPos, obj.PerformanceTrackCue.LateRatio,      'o', 'linestyle', '-', 'color', opts.color.Late, ...
            'markersize', 5, 'linewidth', 1.2, 'markerfacecolor', opts.color.Late,      'markeredgecolor', 'w');
        
        plot(ax, obj.PerformanceTrackUncue.WinPos, obj.PerformanceTrackUncue.CorrectRatio,   'o', 'linestyle', ':', 'color', opts.color.Correct, ...
            'markersize', 4, 'linewidth', 1.2, 'markeredgecolor', opts.color.Correct,   'markerfacecolor', 'w');
        plot(ax, obj.PerformanceTrackUncue.WinPos, obj.PerformanceTrackUncue.WrongRatio,     'o', 'linestyle', ':', 'color', opts.color.Wrong, ...
            'markersize', 4, 'linewidth', 1.2, 'markeredgecolor', opts.color.Wrong,     'markerfacecolor', 'w');
        plot(ax, obj.PerformanceTrackUncue.WinPos, obj.PerformanceTrackUncue.PrematureRatio, 'o', 'linestyle', ':', 'color', opts.color.Premature, ...
            'markersize', 4, 'linewidth', 1.2, 'markeredgecolor', opts.color.Premature, 'markerfacecolor', 'w');
        plot(ax, obj.PerformanceTrackUncue.WinPos, obj.PerformanceTrackUncue.LateRatio,      'o', 'linestyle', ':', 'color', opts.color.Late, ...
            'markersize', 4, 'linewidth', 1.2, 'markeredgecolor', opts.color.Late,      'markerfacecolor', 'w');

        ax.XLabel.String = 'Time in session (s)';
        ax.YLabel.String = 'Performance (%)';
        set(ax, 'XLim', [0 opts.ticks(end)+5], 'YLim', [0 100]);
    end

%% ha4. Plot shuttle time
    function plot_shuttle_time(ax, obj, opts)

        line(ax, [opts.ticks opts.ticks]', [0 2/20], 'color', 'K')

        ST_log = log10(obj.ShuttleTime);
        scatter(ax, opts.ticks, ST_log, ...
            18, 'k', opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        ax.XLabel.String = 'Time in session (s)';
        ax.YLabel.String = 'Shuttle time (s)';
        set(ax, 'xlim', [0 opts.ticks(end)+5], 'ylim', [0 2], 'ticklength', [0.01 0.1], 'ytick', 1:3, 'yticklabel', {'10^1', '10^2', '10^3'});
    end

%% ha5. Plot shuttle time PDF on the right
    function plot_shuttle_time_pdf(ax, obj) 
        
        binEdges = obj.Bins.ShuttleTimeLog;
        ST_log = log10(obj.ShuttleTime);
        PDF_PortL = ksdensity(ST_log, binEdges);

        plot(ax, PDF_PortL, binEdges, 'color', 'k', 'linewidth', 1.5);

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'ylim', [0 3], 'xlimmode', 'auto', 'yticklabel', [], 'ticklength', [0.02 0.1], 'ytick', 1:3)
    end

%% ha6. Plot hold duration and FP
    function plot_hold_duration(ax, obj, opts)

        yline(ax, obj.MixedFP+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        yline(ax, obj.MixedFP+0.8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');
        line(ax, [opts.ticks opts.ticks]', [0 2.5/20], 'color', 'K')

        % FPs
        scatter(ax, opts.ticks, obj.FP, 'LineWidth', 1, 'Marker', '_', 'MarkerEdgeColor', [.7 .7 .7], 'SizeData', 4);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.ticks(obj.Ind.correctL & obj.Ind.cue), ...
            obj.HoldDuration(obj.Ind.correctL & obj.Ind.cue), ...
            24, opts.color.PortL, opts.mk(1),'Markerfacealpha', 0.8, 'linewidth', 1.5, 'MarkerFaceColor', opts.color.PortL);

        % port 2 correct
        scatter(ax, opts.ticks(obj.Ind.correctR & obj.Ind.cue), ...
            obj.HoldDuration(obj.Ind.correctR & obj.Ind.cue), ...
            24, opts.color.PortR, opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5, 'MarkerFaceColor', opts.color.PortR);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.ticks(obj.Ind.correctL & obj.Ind.uncue), ...
            obj.HoldDuration(obj.Ind.correctL & obj.Ind.uncue), ...
            24, opts.color.PortL, opts.mk(1), 'MarkerFaceAlpha', 0.8, 'linewidth', 1.5);

        % port 2 correct
        scatter(ax, opts.ticks(obj.Ind.correctR & obj.Ind.uncue), ...
            obj.HoldDuration(obj.Ind.correctR & obj.Ind.uncue), ...
            24, opts.color.PortR, opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.ticks(obj.Ind.wrongL), ...
            obj.HoldDuration(obj.Ind.wrongL), ...
            32, opts.color.PortL, opts.mk(2),'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % port 2 wrong
        scatter(ax, opts.ticks(obj.Ind.wrongR), ...
            obj.HoldDuration(obj.Ind.wrongR), ...
            32, opts.color.PortR, opts.mk(2), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % premature
        scatter(ax, opts.ticks(obj.Ind.premature), ...
            obj.HoldDuration(obj.Ind.premature), ...
            32, opts.color.Premature, opts.mk(2), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % late
        scatter(ax, opts.ticks(obj.Ind.late), ...
            obj.HoldDuration(obj.Ind.late), ...
            32, opts.color.Late, opts.mk(2), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        ax.XLabel.String = 'Time in session (s)';
        ax.YLabel.String = 'Hold duration (s)';
        set(ax, 'xlim', [0 opts.ticks(end)+5], 'ylim', [0 2.5], 'ticklength', [0.01 0.1]);
    end

%% ha7. plot hold duration pdf
    function plot_hold_duration_pdf_l(ax, obj, opts)

        yline(ax, obj.MixedFP+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        yline(ax, obj.MixedFP+0.8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        binEdges = obj.Bins.HoldDuration;

        for fp = 1:length(obj.MixedFP)
            yline(ax, obj.MixedFP(fp), 'Color', [.7 .7 .7], 'LineWidth', opts.lw(fp), 'LineStyle', '-');
        end

        for i = 1:length(obj.CueUncue)
            for j = 1
                plot(ax, obj.HoldDurationPDF{i, j}, binEdges, 'color', eval("opts.color.Port" + obj.Ports(j)), 'linewidth', 1.5, 'LineStyle', opts.ls(i));
            end
        end

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'xlimmode', 'auto', 'ylim', [0 2.5], 'yticklabel', [], 'ticklength', [0.02 0.1])
    end

%% ha7.1. plot hold duration cdf
    function plot_hold_duration_pdf_r(ax, obj, opts)

        yline(ax, obj.MixedFP+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        yline(ax, obj.MixedFP+0.8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        binEdges = obj.Bins.HoldDuration;

        for fp = 1:length(obj.MixedFP)
            yline(ax, obj.MixedFP(fp), 'Color', [.7 .7 .7], 'LineWidth', opts.lw(fp), 'LineStyle', '-');
        end

        for i = 1:length(obj.CueUncue)
            for j = 2
                plot(ax, obj.HoldDurationPDF{i, j}, binEdges, 'color', eval("opts.color.Port" + obj.Ports(j)), 'linewidth', 1.5, 'LineStyle', opts.ls(i));
            end
        end

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'xlimmode', 'auto', 'ylim', [0 2.5], 'yticklabel', [], 'ticklength', [0.02 0.1])
    end

%% ha8. Plot reaction time
    function plot_reaction_time(ax, obj, opts) 

        line(ax, [opts.ticks opts.ticks]', [0 0.6/20], 'color', 'K')
        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        yline(ax, .8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.ticks(obj.Ind.correctL & obj.Ind.cue), ...
            obj.RT(obj.Ind.correctL & obj.Ind.cue), ...
            24, opts.color.PortL, opts.mk(1),'Markerfacealpha', 0.8, 'linewidth', 1.5, 'MarkerFaceColor', opts.color.PortL);

        % port 2 correct
        scatter(ax, opts.ticks(obj.Ind.correctR & obj.Ind.cue), ...
            obj.RT(obj.Ind.correctR & obj.Ind.cue), ...
            24, opts.color.PortR, opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5, 'MarkerFaceColor', opts.color.PortR);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.ticks(obj.Ind.correctL & obj.Ind.uncue), ...
            obj.RT(obj.Ind.correctL & obj.Ind.uncue), ...
            24, opts.color.PortL, opts.mk(1),'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % port 2 correct
        scatter(ax, opts.ticks(obj.Ind.correctR & obj.Ind.uncue), ...
            obj.RT(obj.Ind.correctR & obj.Ind.uncue), ...
            24, opts.color.PortR, opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.ticks(obj.Ind.wrongL), ...
            obj.RT(obj.Ind.wrongL), ...
            32, opts.color.PortL, opts.mk(2),'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % port 2 wrong
        scatter(ax, opts.ticks(obj.Ind.wrongR), ...
            obj.RT(obj.Ind.wrongR), ...
            32, opts.color.PortR, opts.mk(2), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        ax.XLabel.String = 'Time in session (s)';
        ax.YLabel.String = 'Reaction time (s)';
        set(ax, 'xlim', [0 opts.ticks(end)+5], 'ylim', [0 0.8], 'ticklength', [0.01 0.1]);
    end

%% ha9. plot reaction time PDF for PortL
    function plot_reaction_time_pdf_l(ax, obj, opts)
        
        binEdges = obj.Bins.RT;

        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        yline(ax, .8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        for c = 1:length(obj.CueUncue)
            plot(ax, obj.RTPDF{c, 1}, binEdges, 'color', opts.color.PortL, 'linewidth', 1.5, 'LineStyle', opts.ls(c));
        end

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'xlimmode', 'auto', 'ylim', [0 .8], 'yticklabel', [], 'ticklength', [0.02 0.1]);
    end

%% ha9.1. plot reaction time PDF for PortR
    function plot_reaction_time_pdf_r(ax, obj, opts)
        
        binEdges = obj.Bins.RT;

        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');
        yline(ax, .8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        for c = 1:length(obj.CueUncue)
            plot(ax, obj.RTPDF{c, 2}, binEdges, 'color', opts.color.PortR, 'linewidth', 1.5, 'LineStyle', opts.ls(c));
        end

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'xlimmode', 'auto', 'ylim', [0 .8], 'yticklabel', [], 'ticklength', [0.02 0.1]);
    end

%% ha10. plot movement times
    function plot_movement_time(ax, obj, opts)
        
        line([opts.ticks opts.ticks]', [0 0.1], 'color', 'k')

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.ticks(obj.Ind.correctL & obj.Ind.cue), ...
            obj.MovementTime(obj.Ind.correctL & obj.Ind.cue), ...
            24, opts.color.PortL, opts.mk(1),'Markerfacealpha', 0.8, 'linewidth', 1.5, 'MarkerFaceColor', opts.color.PortL);

        % port 2 correct
        scatter(ax, opts.ticks(obj.Ind.correctR & obj.Ind.cue), ...
            obj.MovementTime(obj.Ind.correctR & obj.Ind.cue), ...
            24, opts.color.PortR, opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5, 'MarkerFaceColor', opts.color.PortR);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.ticks(obj.Ind.correctL & obj.Ind.uncue), ...
            obj.MovementTime(obj.Ind.correctL & obj.Ind.uncue), ...
            24, opts.color.PortL, opts.mk(1),'Markerfacealpha', 0.8, 'linewidth', 1.5);

        % port 2 correct
        scatter(ax, opts.ticks(obj.Ind.correctR & obj.Ind.uncue), ...
            obj.MovementTime(obj.Ind.correctR & obj.Ind.uncue), ...
            24, opts.color.PortR, opts.mk(1), 'Markerfacealpha', 0.8, 'linewidth', 1.5);

        ax.XLabel.String = 'Time in session (s)';
        ax.YLabel.String = 'Movement time (s)';
        set(ax, 'xlim', [0 opts.ticks(end)+5], 'ylim', [0 2], 'ticklength', [0.01 0.1]);
    end

%% ha11. Plot movement time PDF on the right for PortL
    function plot_movement_time_pdf_l(ax, obj, opts)

        binEdges = obj.Bins.MovementTime;

        for c = 1:length(obj.CueUncue)
            plot(ax, obj.MovementTimePDF{c, 1}, binEdges, 'color', opts.color.PortL, 'linewidth', 1.5, 'LineStyle', opts.ls(c))
        end

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'xlimmode', 'auto', 'ylim', [0 2], 'yticklabel', [], 'ticklength', [0.02 0.1]);
    end

%% ha11.1. Plot movement time PDF on the right for PortL
    function plot_movement_time_pdf_r(ax, obj, opts)

        binEdges = obj.Bins.MovementTime;

        for c = 1:length(obj.CueUncue)
            plot(ax, obj.MovementTimePDF{c, 2}, binEdges, 'color', opts.color.PortR, 'linewidth', 1.5, 'LineStyle', opts.ls(c))
        end

        ax.XLabel.String = 'Prob. density (1/s)';
        set(ax, 'xlimmode', 'auto', 'ylim', [0 2], 'yticklabel', [], 'ticklength', [0.02 0.1]);
    end

%% ha12. Performance of each port
    function plot_performance_collected(ax, obj, opts)

        for i = 1:length(obj.CueUncue)
            
            indL = obj.Performance.Cued_this==obj.CueUncue(i) & obj.Performance.TargetPort=="L";
            indR = obj.Performance.Cued_this==obj.CueUncue(i) & obj.Performance.TargetPort=="R";

            line(ax, [1 2.5], [obj.Performance.CorrectRatio(indL) obj.Performance.CorrectRatio(indR)], ...
                'color', opts.color.Correct, 'linewidth', 1.5, 'LineStyle', opts.ls(i));
            line(ax, [3.5 5], [obj.Performance.WrongRatio(indL) obj.Performance.WrongRatio(indR)], ...
                'color', opts.color.Wrong, 'linewidth', 1.5, 'LineStyle', opts.ls(i));
            line(ax, [6 7.5], [obj.Performance.PrematureRatio(indL) obj.Performance.PrematureRatio(indR)], ...
                'color', opts.color.Premature, 'linewidth', 1.5, 'LineStyle', opts.ls(i));
            line(ax, [8.5 10], [obj.Performance.LateRatio(indL) obj.Performance.LateRatio(indR)], ...
                'color', opts.color.Late, 'linewidth', 1.5, 'LineStyle', opts.ls(i));
        end

        ax.YLabel.String = 'Performance (%)';
        set(ax, 'xlim', [0.5 10.5], 'ylim', [0 100], 'ticklength', [0.01 0.1], ...
            'xtick', [1 2.5 3.5 5 6 7.5 8.5 10], 'xticklabel', {'L', 'R', '', '', '', '', '', ''})
    end

%  ha12.1. legend
    function legend_performance(ax, opts)
        
        line(ax, [.5 1], [4 4], 'color', opts.color.Correct,    'linewidth', 1.2);
        line(ax, [.5 1], [3 3], 'color', opts.color.Wrong,      'linewidth', 1.2);
        line(ax, [.5 1], [2 2], 'color', opts.color.Premature,  'linewidth', 1.2);
        line(ax, [.5 1], [1 1], 'color', opts.color.Late,       'linewidth', 1.2);

        text(ax, 1.2, 4, 'Correct',   'fontsize', 7);
        text(ax, 1.2, 3, 'Wrong',     'fontSize', 7);
        text(ax, 1.2, 2, 'Premature', 'fontsize', 7);
        text(ax, 1.2, 1, 'Late',      'fontsize', 7);

        set(ax, 'ylim', [0.5 4.5], 'xlim', [.5 2], 'xcolor', 'none', 'ycolor', 'none', 'color', 'none');
    end

%% ha13. Reaction time to each port
    function plot_reaction_time_collected(ax, obj, opts) 

        thisRT = nan(sum(obj.Stage==1 & obj.Ind.correct), length(obj.Ports)*length(obj.CueUncue));
        for p = 1:length(obj.Ports)
            for fp = 1:length(obj.CueUncue)
                thisRT(1:length(obj.RTSorted{fp, p}), fp+(p-1)*length(obj.CueUncue)) = obj.RTSorted{fp,p};
            end
        end
        thisRT(1, all(isnan(thisRT))) = 0;

        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');

        violinplot({thisRT(:, 3:4), thisRT(:, 1:2)}, {'Cue', 'Uncue'}, ...
            'ViolinColor', {repmat(opts.color.PortR, 2, 1), repmat(opts.color.PortL, 2, 1)}, 'ScatterSize', 6, ...
            'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false);

        scatter(ax, [1 2 1 2], median(thisRT, 'omitnan'), 24, [repmat(opts.color.PortL, 2, 1); repmat(opts.color.PortR, 2, 1)], 'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3])

        ax.YLabel.String = 'Reaction time (s)';
        set(ax, 'xlim', [0 4], 'ylim', [0 .8], 'xtick', [1 2], 'xticklabel', {'Cue', 'Uncue'}, 'ticklength', [0.01 0.1], 'box', 'off');
    end

end
