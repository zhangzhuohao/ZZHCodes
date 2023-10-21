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
set(gcf, 'unit', 'centimeters', 'position', [2 1.5 39.5 17.3], 'paperpositionmode', 'auto', 'color', 'w');

uicontrol('Style', 'text', 'parent', 23, 'units', 'normalized', 'position', [0.3 0.96 0.4 0.03],...
    'string', obj.Subject+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end), 'fontsize', 11, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set axes and plot
ha_perf_track_l = axes;
set(ha_perf_track_l, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_track(ha_perf_track_l, obj, opts, "L");

ha_perf_track_r = axes;
set(ha_perf_track_r, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_track(ha_perf_track_r, obj, opts, "R");
set(ha_perf_track_r, 'yticklabel', [], 'ylabel', []);

% Performance progress
ha_perf_l = axes;
set(ha_perf_l, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_progress(ha_perf_l, obj, opts, "L");
set(ha_perf_l, 'xticklabel', [], 'xlabel', [], 'Title', []);

ha_perf_r = axes;
set(ha_perf_r, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_performance_progress(ha_perf_r, obj, opts, "R");
set(ha_perf_r, 'xticklabel', [], 'xlabel', [], 'yticklabel', [], 'ylabel', [], 'Title', []);

% Hold duration scatter
ha_hd_scatter = axes;
set(ha_hd_scatter, 'units', 'centimeters', 'position', [1.5+2*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_scatter(ha_hd_scatter, obj, opts, 'uncue');

% Hold duration violin plot
ha_hd_violin = axes;
set(ha_hd_violin, 'units', 'centimeters', 'position', [1.5+3*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_violin(ha_hd_violin, obj, opts, 'uncue');

% Heatmap
ha_hd_heatmap_l = axes;
set(ha_hd_heatmap_l, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_heatmap(ha_hd_heatmap_l, obj, opts, "L", 'uncue')
set(ha_hd_heatmap_l, 'Title', []);

ha_hd_heatmap_r = axes;
set(ha_hd_heatmap_r, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_heatmap(ha_hd_heatmap_r, obj, opts, "R", 'uncue')
set(ha_hd_heatmap_r, 'yticklabel', [], 'ylabel', [], 'Title', []);

ha_hd_heatmap_l.CLim(2) = 5;
ha_hd_heatmap_r.CLim(2) = ha_hd_heatmap_l.CLim(2);

cb_heatmap_cb = colorbar(ha_hd_heatmap_r, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1))-opts.sep_row+0.2 2+0*(opts.sep_col+opts.plotsize1(2)), [.3 opts.plotsize1(2)]]);
cb_heatmap_cb.Limits = [0 cb_heatmap_cb.Limits(2)];
cb_heatmap_cb.Label.String = "Prob. density (1/s)";
cb_heatmap_cb.FontSize = 8;

% density
ha_hd_pdf_l = axes;
set(ha_hd_pdf_l, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize1(1)) 1.2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_pdf_early_late(ha_hd_pdf_l, obj, opts, "L");

ha_hd_pdf_r = axes;
set(ha_hd_pdf_r, 'units', 'centimeters', 'position', [1.5+3*(opts.sep_row+opts.plotsize1(1)) 1.2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_pdf_early_late(ha_hd_pdf_r, obj, opts, "R");

density_max = max([ha_hd_pdf_l.YLim(2) ha_hd_pdf_r.YLim(2)]);
set(ha_hd_pdf_l, 'ylim', [0 density_max]);
set(ha_hd_pdf_r, 'ylim', [0 density_max], 'yticklabel', [], 'ylabel', []);

ha_e_l_legend = axes;
set(ha_e_l_legend, 'units', 'centimeters', 'position', [1.2+3*(opts.sep_row+opts.plotsize1(1))+5*opts.plotsize2(1)/6 .7*(opts.sep_col+opts.plotsize1(2))+1*opts.plotsize2(2)/6, opts.plotsize2 ./ [5 3]], ...
    'nextplot', 'add', 'fontsize', 8, 'color', 'none', 'xcolor', 'none', 'ycolor', 'none', 'xlim', [0 3], 'ylim', [0 3]);
line(ha_e_l_legend, [1 2], [2 2], 'Color', opts.color.PhaseEarly, 'LineWidth', 1.5, 'LineStyle', '-');
line(ha_e_l_legend, [1 2], [1 1], 'Color', opts.color.PhaseLate, 'LineWidth', 1.5, 'LineStyle', '-');
text(ha_e_l_legend, 2.2, 2, "Early", 'FontSize', 8, 'Color', opts.color.PhaseEarly, 'FontWeight', 'bold');
text(ha_e_l_legend, 2.2, 1, "Late", 'FontSize', 8, 'Color', opts.color.PhaseLate, 'FontWeight', 'bold');

%
ha_hd_iqr = axes;
set(ha_hd_iqr, 'units', 'centimeters', 'position', [1.5+2*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_iqr(ha_hd_iqr, obj, opts);

ha_hd_median = axes;
set(ha_hd_median, 'units', 'centimeters', 'position', [1.5+3*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8, 'tickdir', 'out');
plot_hold_duration_median(ha_hd_median, obj, opts);

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

    end

%% ha5. Plot performance progress
    function plot_performance_progress(ax, obj, opts, port)

        session_id = 1:obj.NumSessions;
        xline(ax, session_id-0.5, 'Color', [.9 .9 .9], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        c = 2;
        ind_this = find(obj.Performance.Cued_this==obj.CueUncue(c) & obj.Performance.TargetPort==port);

        plot(ax, session_id, obj.Performance.PrematureRatio(ind_this), '-o', 'color', [opts.color.Premature .8], ...
            'markersize', 4, 'linewidth', 2, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.LateRatio(ind_this), '-o', 'color', [opts.color.Late .8], ...
            'markersize', 4, 'linewidth', 2, 'markerfacecolor', opts.color.Late, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.WrongRatio(ind_this), '-o', 'color', [opts.color.Wrong .8], ...
            'markersize', 4, 'linewidth', 2, 'markerfacecolor', opts.color.Wrong, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.CorrectRatio(ind_this), '-o', 'color', [opts.color.Correct .8], ...
            'markersize', 4, 'linewidth', 2, 'markerfacecolor', opts.color.Correct, 'markeredgecolor', 'w');

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
            ind_uncue = find(obj.BehavTable.SessionDate==obj.Sessions(s_this) & eval("obj.Ind.port"+upper(port)) & obj.Ind.uncue);
            this_track_uncue = obj.getPerformanceTrack(obj.BehavTable(ind_uncue, :), obj.Ind(ind_uncue, :));
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.CorrectRatio,   'linestyle', '-', 'color', [opts.color.Correct .7], 'linewidth', 1);
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.WrongRatio,     'linestyle', '-', 'color', [opts.color.Wrong .7], 'linewidth', 1);
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.PrematureRatio, 'linestyle', '-', 'color', [opts.color.Premature .7], 'linewidth', 1);
            plot(ax, this_track_uncue.WinPos+s_sep(s_this), this_track_uncue.LateRatio,      'linestyle', '-', 'color', [opts.color.Late .7], 'linewidth', 1);
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

        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 obj.Bins.HoldDuration(end)], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_hold_duration_median(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.9 .9 .9], 'LineWidth', 1, 'LineStyle', '-');
        yline(ax, obj.MixedFP, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');
        yline(ax, obj.MixedFP+0.8, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', ':');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for p_this = 1:length(obj.Ports)
            plot(1:obj.NumSessions, ...
                obj.HDStat.Median(obj.HDStat.Port==obj.Ports(p_this) & obj.HDStat.thisCued==0), ...
                'Color', eval("opts.color.Port"+obj.Ports(p_this)), 'LineWidth', 1.2, 'LineStyle', '-', ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', eval("opts.color.Port"+obj.Ports(p_this)), 'MarkerEdgeColor', 'w');
        end

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Reaction time median (s)';
        ax.XLabel.FontWeight = 'Bold';
        ax.YLabel.FontWeight = 'Bold';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 3], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_hold_duration_iqr(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.9 .9 .9], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);

        for p_this = 1:length(obj.Ports)
            plot(1:obj.NumSessions, ...
                obj.HDStat.IQR(obj.HDStat.Port==obj.Ports(p_this) & obj.HDStat.thisCued==0), ...
                'Color', eval("opts.color.Port"+obj.Ports(p_this)), 'LineWidth', 1.2, 'LineStyle', '-', ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', eval("opts.color.Port"+obj.Ports(p_this)), 'MarkerEdgeColor', 'w');
        end

        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

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

        if all(all(cellfun(@(x) length(x)>=100, obj.HDSortedNone)))
            HDs = obj.HDSortedNone;
        else
            HDs = cellfun(@(x, y) [x; y], obj.HDSortedNone, obj.HDSortedSaline, 'UniformOutput', false);
        end

        xline(ax, fp_this, 'color', [.7 .7 .7], 'linewidth', 1.5, 'LineStyle', '-');
        xline(ax, fp_this+.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');
        xline(ax, fp_this+.8, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', ':');

        hd_uncue = HDs{2, p_this};
        if length(hd_uncue) >= 2*obj.PhaseCount
            hd_uncue_early_pdf = ksdensity(hd_uncue(1:obj.PhaseCount), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_uncue_late_pdf  = ksdensity(hd_uncue(end-obj.PhaseCount+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');
        else
            l_uncue = length(hd_uncue);
            hd_uncue_early_pdf = ksdensity(hd_uncue(1:floor(l_uncue/2)), obj.Bins.HoldDuration, 'Function', 'pdf');
            hd_uncue_late_pdf  = ksdensity(hd_uncue(end-floor(l_uncue/2)+1:end), obj.Bins.HoldDuration, 'Function', 'pdf');
        end

        plot(ax, obj.Bins.HoldDuration, hd_uncue_early_pdf, ...
            'color', opts.color.PhaseEarly, 'linewidth', 2, 'LineStyle', '-');
        plot(ax, obj.Bins.HoldDuration, hd_uncue_late_pdf, ...
            'color', opts.color.PhaseLate, 'linewidth', 2, 'LineStyle', '-');
        
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
