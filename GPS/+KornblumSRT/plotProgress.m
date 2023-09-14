function fig = plotProgress(obj)

%%
opts.color = GPSColor(); % Color class for GPS
opts.mk = ["o", "x"]; % Scatter marker for correct and others
opts.ls = ["-", ":"]; % Line style for [ShortFP, MedFP, LongFP]
opts.lw = [0.5, 1, 1.5]; % Line width for [ShortFP, MedFP, LongFP]
opts.trial_ticks = obj.BehavTable.TrialStartTimeProgress + obj.BehavTable.CentInTime;

opts.session_sep = zeros(1, obj.NumSessions-1);
for s = 1:obj.NumSessions-1
    ind = find(obj.BehavTable.SessionDate==obj.Sessions(s), 1, 'last');
    opts.session_sep(s) = opts.trial_ticks(ind)+5;
end

opts.session_date = char(obj.Sessions);
opts.session_date = opts.session_date(:, 5:8);

opts.plotsize1 = [8 4];
opts.plotsize2 = [8 3.5];
opts.plotsize3 = [3.5 3];
opts.plotsize4 = [2 3];

opts.sep_col = 1;
opts.sep_row = 1.5;

%%
fig = figure(23); clf(23);
set(gcf, 'unit', 'centimeters', 'position', [2 1.5 39.5 22.5], 'paperpositionmode', 'auto', 'color', 'w');

mycolormap = customcolormap_preset("red-white-blue");

uicontrol('Style', 'text', 'parent', 23, 'units', 'normalized', 'position', [0.3 0.96 0.4 0.03],...
    'string', obj.Subject+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end), 'fontsize', 11, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set axes and plot
ha_perf_track_l = axes;
set(ha_perf_track_l, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_track(ha_perf_track_l, obj, opts, "L");

ha_perf_track_r = axes;
set(ha_perf_track_r, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_track(ha_perf_track_r, obj, opts, "R");
set(ha_perf_track_r, 'yticklabel', [], 'ylabel', []);


% Performance progress
ha1 = axes;
set(ha1, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 1.5+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_progress(ha1, obj, opts, "L");

ha2 = axes;
set(ha2, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 1.5+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_progress(ha2, obj, opts, "R");
set(ha2, 'yticklabel', [], 'ylabel', []);

% Hold duration scatter
ha3 = axes;
set(ha3, 'units', 'centimeters', 'position', [1.5+2*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_scatter(ha3, obj, opts, 'cue');
set(ha3, 'xtick', [], 'xlabel', []);

ha4 = axes;
set(ha4, 'units', 'centimeters', 'position', [1.5+2*(opts.sep_row+opts.plotsize1(1)) 2.2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_scatter(ha4, obj, opts, 'uncue');

% Hold duration violin plot
ha5 = axes;
set(ha5, 'units', 'centimeters', 'position', [1+3*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_violin(ha5, obj, opts, 'cue');
set(ha5, 'xtick', [], 'xlabel', [], 'yticklabel', [], 'ylabel', []);

ha6 = axes;
set(ha6, 'units', 'centimeters', 'position', [1+3*(opts.sep_row+opts.plotsize1(1)) 2.2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_violin(ha6, obj, opts, 'uncue');
set(ha6, 'yticklabel', [], 'ylabel', []);

% Heatmap
ha7 = axes;
set(ha7, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) .7+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_heatmap(ha7, obj, opts, "L", 'cue')
set(ha7, 'xtick', [], 'xlabel', []);

ha8 = axes;
set(ha8, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 1.2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_heatmap(ha8, obj, opts, "L", 'uncue')
set(ha8, 'Title', []);

ha9 = axes;
set(ha9, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) .7+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_heatmap(ha9, obj, opts, "R", 'cue')
set(ha9, 'xtick', [], 'xlabel', [], 'yticklabel', [], 'ylabel', []);

ha10 = axes;
set(ha10, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 1.2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_heatmap(ha10, obj, opts, "R", 'uncue')
set(ha10, 'yticklabel', [], 'ylabel', [], 'Title', []);

% ha7.CLim(2) = max([ha7.CLim(2) ha8.CLim(2) ha9.CLim(2) ha10.CLim(2)]);
% ha8.CLim(2) = ha7.CLim(2);
% ha9.CLim(2) = ha7.CLim(2);
% ha10.CLim(2) = ha7.CLim(2);

%
ha121 = axes;
set(ha121, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize1(1)) 1.2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_pdf_early_late(ha121, obj, opts, "L");

ha122 = axes;
set(ha122, 'units', 'centimeters', 'position', [1.5+3*(opts.sep_row+opts.plotsize1(1)) 1.2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_pdf_early_late(ha122, obj, opts, "R");

density_max = max([ha121.YLim(2) ha122.YLim(2)]);
set(ha121, 'ylim', [0 density_max]);
set(ha122, 'ylim', [0 density_max], 'yticklabel', [], 'ylabel', []);

ha123 = axes;
set(ha123, 'units', 'centimeters', 'position', [1.2+3*(opts.sep_row+opts.plotsize1(1))+5*opts.plotsize2(1)/6 .7*(opts.sep_col+opts.plotsize1(2))+1*opts.plotsize2(2)/6, opts.plotsize2 ./ [5 3]], ...
    'nextplot', 'add', 'fontsize', 8, 'color', 'none', 'xcolor', 'none', 'ycolor', 'none', 'xlim', [0 3], 'ylim', [0 3]);
line(ha123, [1 2], [2 2], 'Color', opts.color.PhaseEarly, 'LineWidth', 1.5, 'LineStyle', '-');
line(ha123, [1 2], [1 1], 'Color', opts.color.PhaseLate, 'LineWidth', 1.5, 'LineStyle', '-');
text(ha123, 2.2, 2, "Early", 'FontSize', 8, 'Color', opts.color.PhaseEarly, 'FontWeight', 'bold');
text(ha123, 2.2, 1, "Late", 'FontSize', 8, 'Color', opts.color.PhaseLate, 'FontWeight', 'bold');

%
ha7.CLim(2) = 5;
ha8.CLim(2) = ha7.CLim(2);
ha9.CLim(2) = ha7.CLim(2);
ha10.CLim(2) = ha7.CLim(2);

cb10 = colorbar(ha10, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1))-opts.sep_row+0.2 1.2+0*(opts.sep_col+opts.plotsize1(2)), [.3 opts.plotsize1(2)]]);
cb10.Limits = [0 cb10.Limits(2)];
cb10.Label.String = "Prob. density (1/s)";
cb10.FontSize = 8;

%
ha14 = axes;
set(ha14, 'units', 'centimeters', 'position', [1.5+2*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_iqr(ha14, obj, opts);

ha15 = axes;
set(ha15, 'units', 'centimeters', 'position', [1.5+3*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_reaction_time_median(ha15, obj, opts);

ha151 = axes;
set(ha151, 'units', 'centimeters', 'position', [1.2+3*(opts.sep_row+opts.plotsize1(1))+5*opts.plotsize1(1)/6 1.7*(opts.sep_col+opts.plotsize1(2))+1*opts.plotsize1(2)/6, opts.plotsize2 ./ [5 3]], ...
    'nextplot', 'add', 'fontsize', 8, 'color', 'none', 'xcolor', 'none', 'ycolor', 'none', 'xlim', [0 3], 'ylim', [0 3]);
line(ha151, [1 2], [2 2], 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
line(ha151, [1 2], [1 1], 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
line(ha151, 1.5, 2, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-', 'marker', 'o', 'markersize', 6, 'markerfacecolor', 'k', 'markeredgecolor', 'w');
line(ha151, 1.5, 1, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':', 'marker', 'o', 'markersize', 4, 'markeredgecolor', 'k', 'markerfacecolor', 'w');
text(ha151, 2.2, 2, "Cue", 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
text(ha151, 2.2, 1, "Uncue", 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');

%% ha2. Plot every trial's hold duration during training
    function plot_hold_duration_scatter(ax, obj, opts, cued)

        fp_this = unique(obj.MixedFP);
        yline(ax, fp_this, 'Color', [.8 .8 .8], 'LineWidth', 1.2, 'LineStyle', '-', 'Alpha', 0.4);
        xline(ax, opts.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

        s_sep = [0 opts.session_sep opts.trial_ticks(end)];
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [s_sep(dcz_ind); s_sep(dcz_ind); s_sep(dcz_ind+1); s_sep(dcz_ind+1)], 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        switch cued
            case {'cue'}
                yline(ax, fp_this+.5, 'Color', [.8 .8 .8], 'LineWidth', 1.2, 'LineStyle', '--', 'Alpha', 0.8);
            case {'uncue'}
                yline(ax, fp_this+.8, 'Color', [.8 .8 .8], 'LineWidth', 1.2, 'LineStyle', ':', 'Alpha', 0.8);
        end

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.trial_ticks(obj.Ind.wrongL & obj.Ind.(cued)), ...
            obj.BehavTable.HoldDuration(obj.Ind.wrongL & obj.Ind.(cued)), ...
            32, opts.color.PortL, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 wrong
        scatter(ax, opts.trial_ticks(obj.Ind.wrongR & obj.Ind.(cued)), ...
            obj.BehavTable.HoldDuration(obj.Ind.wrongR & obj.Ind.(cued)), ...
            32, opts.color.PortR, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % premature
        scatter(ax, opts.trial_ticks(obj.Ind.premature & obj.Ind.(cued)), ...
            obj.BehavTable.HoldDuration(obj.Ind.premature & obj.Ind.(cued)), ...
            32, opts.color.Premature, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % late
        scatter(ax, opts.trial_ticks(obj.Ind.late & obj.Ind.(cued)), ...
            obj.BehavTable.HoldDuration(obj.Ind.late & obj.Ind.(cued)), ...
            32, opts.color.Late, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.trial_ticks(obj.Ind.correctL & obj.Ind.(cued)), ...
            obj.BehavTable.HoldDuration(obj.Ind.correctL & obj.Ind.(cued)), ...
            24, opts.color.PortL, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 correct
        scatter(ax, opts.trial_ticks(obj.Ind.correctR & obj.Ind.(cued)), ...
            obj.BehavTable.HoldDuration(obj.Ind.correctR & obj.Ind.(cued)), ...
            24, opts.color.PortR, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Hold duration (s)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        set(ax, 'xlim', [0 opts.trial_ticks(end)+5], 'ylim', [0 obj.Bins.HoldDuration(end)], 'ticklength', [0.01 0.1]);

        t = cued;
        t(1) = upper(t(1));
        ax.Title.String = t;
        ax.Title.FontWeight = 'Bold';

    end

%% ha5. Plot performance progress
    function plot_performance_progress(ax, obj, opts, port)

        session_id = 1:obj.NumSessions;
        xline(ax, session_id-0.5, 'Color', [.9 .9 .9], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        c = 1;
        ind_this = find(obj.Performance.Cued_this==obj.CueUncue(c) & obj.Performance.TargetPort==port);

        plot(ax, session_id, obj.Performance.PrematureRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', [opts.color.Premature .2], ...
            'markersize', 6, 'linewidth', 2, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.LateRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', [opts.color.Late .2], ...
            'markersize', 6, 'linewidth', 2, 'markerfacecolor', opts.color.Late, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.WrongRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', [opts.color.Wrong .2], ...
            'markersize', 6, 'linewidth', 2, 'markerfacecolor', opts.color.Wrong, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.CorrectRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', [opts.color.Correct .2], ...
            'markersize', 6, 'linewidth', 2, 'markerfacecolor', opts.color.Correct, 'markeredgecolor', 'w');

        c = 2;
        ind_this = find(obj.Performance.Cued_this==obj.CueUncue(c) & obj.Performance.TargetPort==port);

        plot(ax, session_id, obj.Performance.PrematureRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', .8*opts.color.Premature, ...
            'markersize', 4, 'linewidth', 2, 'markeredgecolor', .8*opts.color.Premature, 'markerfacecolor', 'w');
        plot(ax, session_id, obj.Performance.LateRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', .8*opts.color.Late, ...
            'markersize', 4, 'linewidth', 2, 'markeredgecolor', .8*opts.color.Late, 'markerfacecolor', 'w');
        plot(ax, session_id, obj.Performance.WrongRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', .8*opts.color.Wrong, ...
            'markersize', 4, 'linewidth', 2, 'markeredgecolor', .8*opts.color.Wrong, 'markerfacecolor', 'w');
        plot(ax, session_id, obj.Performance.CorrectRatio(ind_this), 'o', 'linestyle', opts.ls(c), 'color', .8*opts.color.Correct, ...
            'markersize', 4, 'linewidth', 2, 'markeredgecolor', .8*opts.color.Correct, 'markerfacecolor', 'w');

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Performance (%)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        switch port
            case {'L'}
                ax.Title.String = 'Left';
                ax.Title.Color = opts.color.PortL;
            case {'R'}
                ax.Title.String = 'Right';
                ax.Title.Color = opts.color.PortR;
        end
        ax.Title.FontWeight = 'Bold';

        set(ax, 'xlim', [0.5 obj.NumSessions+0.5], 'xtick', session_id, 'xticklabel', opts.session_date);
    end

    function plot_performance_track(ax, obj, opts, port)

        session_sep = zeros(1, obj.NumSessions-1);
        for i = 1:obj.NumSessions-1
             win_pos = obj.BehavTable.TrialStartTimeProgress(obj.BehavTable.SessionDate==obj.Sessions(i)) + 5;
             session_sep(i) = win_pos(end);
        end
        xline(ax, session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-')

        s_sep = [0 session_sep];

        for s_this = 1:obj.NumSessions
            ind_cue = find(obj.BehavTable.SessionDate==obj.Sessions(s_this) & eval("obj.Ind.port"+upper(port)) & obj.Ind.cue);
            this_track_cue = obj.getPerformanceTrack(obj.BehavTable(ind_cue, :), obj.Ind(ind_cue, :));
            plot(ax, this_track_cue.WinPos+s_sep(s_this), this_track_cue.CorrectRatio,   'linestyle', '-', 'color', [opts.color.Correct .2], 'linewidth', 1.8);
            plot(ax, this_track_cue.WinPos+s_sep(s_this), this_track_cue.WrongRatio,     'linestyle', '-', 'color', [opts.color.Wrong .2], 'linewidth', 1.8);
            plot(ax, this_track_cue.WinPos+s_sep(s_this), this_track_cue.PrematureRatio, 'linestyle', '-', 'color', [opts.color.Premature .2], 'linewidth', 1.8);
            plot(ax, this_track_cue.WinPos+s_sep(s_this), this_track_cue.LateRatio,      'linestyle', '-', 'color', [opts.color.Late .1], 'linewidth', 1.8);

            ind_uncue = find(obj.BehavTable.SessionDate==obj.Sessions(s_this) & eval("obj.Ind.port"+upper(port)) & obj.Ind.uncue);
            this_track_uncue = obj.getPerformanceTrack(obj.BehavTable(ind_uncue, :), obj.Ind(ind_uncue, :));
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.CorrectRatio,   'linestyle', ':', 'color', .8*opts.color.Correct, 'linewidth', 1.8);
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.WrongRatio,     'linestyle', ':', 'color', .8*opts.color.Wrong, 'linewidth', 1.8);
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.PrematureRatio, 'linestyle', ':', 'color', .8*opts.color.Premature, 'linewidth', 1.8);
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.LateRatio,      'linestyle', ':', 'color', .8*opts.color.Late, 'linewidth', 1.8);
        end

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Performance (%)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';

        switch port
            case {'L'}
                ax.Title.String = 'Left';
                ax.Title.Color = opts.color.PortL;
            case {'R'}
                ax.Title.String = 'Right';
                ax.Title.Color = opts.color.PortR;
        end
        ax.Title.FontWeight = 'Bold';

        set(ax, 'XLim', [0 obj.BehavTable.TrialStartTimeProgress(end)], 'YLim', [0 100]);
    end

%% ha6. Plot hold duration of each session (violin), each port
    function plot_hold_duration_violin(ax, obj, opts, cued)

        fp_this = unique(obj.MixedFP);
        yline(ax, fp_this, 'Color', [.8 .8 .8], 'LineWidth', 1.2, 'LineStyle', '-', 'Alpha', 0.4);
        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

        switch cued
            case {'cue'}
                yline(ax, fp_this+.5, 'Color', [.8 .8 .8], 'LineWidth', 1.2, 'LineStyle', '--', 'Alpha', 0.8);
            case {'uncue'}
                yline(ax, fp_this+.8, 'Color', [.8 .8 .8], 'LineWidth', 1.2, 'LineStyle', ':', 'Alpha', 0.8);
        end
        
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        thisHD = nan(height(obj.BehavTable), length(obj.Ports)*obj.NumSessions);

        HD = obj.BehavTable.HoldDuration;
        [~, ~, indrmv] = rmoutliers_custome(HD);
        HD(indrmv) = nan;

        for s_this = 1:obj.NumSessions
            session_id = obj.BehavTable.SessionDate==obj.Sessions(s_this);
            for p_this = 1:length(obj.Ports)
                ind_this = session_id & obj.BehavTable.Stage==1 & eval("obj.Ind.port"+obj.Ports(p_this)) & obj.Ind.(cued);
                thisHD(ind_this, 2*(s_this-1)+p_this) = HD(ind_this);
            end
        end
        thisHD(1, all(isnan(thisHD))) = 0;

        violinplot({thisHD(:, 2:2:end), thisHD(:, 1:2:end)}, obj.Sessions, ...
            'ViolinColor', {repmat(opts.color.PortR, obj.NumSessions, 1), repmat(opts.color.PortL, obj.NumSessions, 1)}, ...
            'ScatterSize', 4, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false, 'Bandwidth', 0.1);

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Hold duration (s)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        set(ax, 'xlim', [0 opts.trial_ticks(end)+5], 'ylim', [0 2.5], 'ticklength', [0.01 0.1]);

        t = cued;
        t(1) = upper(t(1));
        ax.Title.String = t;
        ax.Title.FontWeight = 'Bold';

        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 obj.Bins.HoldDuration(end)], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end


%%
    function plot_reaction_time_median(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.9 .9 .9], 'LineWidth', 1, 'LineStyle', '-');
        yline(ax, 0.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        yline(ax, 0.8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for p_this = 1:length(obj.Ports)
            plot(1:obj.NumSessions, ...
                obj.RTStat.Median(obj.RTStat.Port==obj.Ports(p_this) & obj.RTStat.thisCued==1), ...
                'Color', [eval("opts.color.Port"+obj.Ports(p_this)) .2], 'LineWidth', 1.5, 'LineStyle', '-', ...
                'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', eval("opts.color.Port"+obj.Ports(p_this)), 'MarkerEdgeColor', 'w');
            plot(1:obj.NumSessions, ...
                obj.RTStat.Median(obj.RTStat.Port==obj.Ports(p_this) & obj.RTStat.thisCued==0), ...
                'Color', eval("opts.color.Port"+obj.Ports(p_this)), 'LineWidth', 1.5, 'LineStyle', ':', ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', eval("opts.color.Port"+obj.Ports(p_this)), 'MarkerFaceColor', 'w');
        end

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Reaction time median (s)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 .8], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_hold_duration_iqr(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.9 .9 .9], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for p_this = 1:length(obj.Ports)
            plot(1:obj.NumSessions, ...
                obj.HDStat.IQR(obj.HDStat.Port==obj.Ports(p_this) & obj.HDStat.thisCued==1), ...
                'Color', [eval("opts.color.Port"+obj.Ports(p_this)) .2], 'LineWidth', 1.5, 'LineStyle', '-', ...
                'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', eval("opts.color.Port"+obj.Ports(p_this)), 'MarkerEdgeColor', 'w');
            plot(1:obj.NumSessions, ...
                obj.HDStat.IQR(obj.HDStat.Port==obj.Ports(p_this) & obj.HDStat.thisCued==0), ...
                'Color', eval("opts.color.Port"+obj.Ports(p_this)), 'LineWidth', 1.5, 'LineStyle', ':', ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', eval("opts.color.Port"+obj.Ports(p_this)), 'MarkerFaceColor', 'w');
        end

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Hold duration IQR (s)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 max(obj.HDStat.IQR)*1.2], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_hold_duration_heatmap(ax, obj, opts, port, cued)

        hd_l = zeros(length(obj.Bins.HoldDuration), obj.NumSessions);
        hd_r = zeros(length(obj.Bins.HoldDuration), obj.NumSessions);

        this_fp = unique(obj.MixedFP);

        for s_this = 1:obj.NumSessions
            switch cued
                case {'cue'}
                    hd_l(:, s_this) = obj.HDPDF{s_this}{1, 1};
                    hd_r(:, s_this) = obj.HDPDF{s_this}{1, 2};
                    yline(ax, (this_fp+.5)/obj.Bins.width, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
                    ax.YLabel.String = 'Hold duration (s) - Cue';
                case {'uncue'}
                    hd_l(:, s_this) = obj.HDPDF{s_this}{2, 1};
                    hd_r(:, s_this) = obj.HDPDF{s_this}{2, 2};
                    yline(ax, (this_fp+.8)/obj.Bins.width, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');
                    ax.YLabel.String = 'Hold duration (s) - Uncue';
            end
        end

        switch port
            case {"L", "l"}
                imagesc(ax, hd_l);
                clim(ax, [0 max(max([hd_l; hd_r]))]);
            case {"R", "r"}
                imagesc(ax, hd_r);
                clim(ax, [0 max(max([hd_l; hd_r]))]);
        end

        yline(ax, this_fp/obj.Bins.width, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '-');

        for s_this = 1:obj.NumSessions
            if obj.Label(s_this)=="Chemo"
                patch(ax, 'XData', s_this+[-0.5 0.5 0.5 -0.5], 'YData', [-.5 -.5 length(obj.Bins.HoldDuration)+.5 length(obj.Bins.HoldDuration)+.5], 'FaceColor', 'none', ...
                    'EdgeColor', opts.color.Treat, 'LineWidth', 2, 'LineStyle', ':', 'EdgeAlpha', 0.8);
            end
        end

        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        ax.XLabel.String = 'Sessions';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [-.5 length(obj.Bins.HoldDuration)+.5], ...
            'xtick', 1:obj.NumSessions, 'xticklabel', opts.session_date, ...
            'ytick', 0:250:length(obj.Bins.HoldDuration), 'yticklabel', string(0:0.5:obj.Bins.HoldDuration(end)), ...
            'ticklength', [0.01 0.1], 'box', 'off');
        switch port
            case {'L'}
                ax.Title.String = 'Left';
                ax.Title.Color = opts.color.PortL;
            case {'R'}
                ax.Title.String = 'Right';
                ax.Title.Color = opts.color.PortR;
        end
        ax.Title.FontWeight = 'Bold';

        colormap(ax, 'turbo');
    end

%%
    function plot_hold_duration_pdf_early_late(ax, obj, opts, port)

        p_this = obj.Ports==upper(string(port));
        fp_this = unique(obj.MixedFP);

        HDs = obj.HDSortedAll;

        xline(ax, fp_this, 'color', [.7 .7 .7], 'linewidth', 1.5, 'LineStyle', '-');
        xline(ax, fp_this+.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        xline(ax, fp_this+.8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        hd_cue   = HDs{1, p_this};
        if length(hd_cue) >= 2*obj.PhaseCount
            hd_cue_early_pdf = ksdensity(hd_cue(1:obj.PhaseCount), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_cue_late_pdf  = ksdensity(hd_cue(end-obj.PhaseCount+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');
        else
            l_cue = length(hd_cue);
            hd_cue_early_pdf = ksdensity(hd_cue(1:floor(l_cue/2)), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_cue_late_pdf  = ksdensity(hd_cue(end-floor(l_cue/2)+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');
        end
        
        hd_uncue = HDs{2, p_this};
        if length(hd_uncue) >= 2*obj.PhaseCount
            hd_uncue_early_pdf = ksdensity(hd_uncue(1:obj.PhaseCount), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_uncue_late_pdf  = ksdensity(hd_uncue(end-obj.PhaseCount+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');
        else
            l_uncue = length(hd_uncue);
            hd_uncue_early_pdf = ksdensity(hd_uncue(1:floor(l_uncue/2)), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_uncue_late_pdf  = ksdensity(hd_uncue(end-floor(l_uncue/2)+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');
        end

        plot(ax, obj.Bins.HoldDuration, hd_cue_early_pdf, ...
            'color', [opts.color.PhaseEarly .3], 'linewidth', 2, 'LineStyle', '-');
        plot(ax, obj.Bins.HoldDuration, hd_cue_late_pdf, ...
            'color', [opts.color.PhaseLate .3], 'linewidth', 2, 'LineStyle', '-');
        plot(ax, obj.Bins.HoldDuration, hd_uncue_early_pdf, ...
            'color', .9*opts.color.PhaseEarly, 'linewidth', 2, 'LineStyle', ':');
        plot(ax, obj.Bins.HoldDuration, hd_uncue_late_pdf, ...
            'color', .9*opts.color.PhaseLate, 'linewidth', 2, 'LineStyle', ':');
        
        ax.XLabel.String = "Hold duration (s)";
        ax.YLabel.String = "Prob. density (1/s)";
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        switch port
            case {'L'}
                ax.Title.String = 'Left';
                ax.Title.Color = opts.color.PortL;
            case {'R'}
                ax.Title.String = 'Right';
                ax.Title.Color = opts.color.PortR;
        end
        ax.Title.FontWeight = 'Bold';
        set(ax, 'xlim', [0 3], 'ylimmode', 'auto');
    end

end
