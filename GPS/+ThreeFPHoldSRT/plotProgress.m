function fig = plotProgress(obj)

%%
opts.color = GPSColor(); % Color class for GPS
opts.mk = ["o", "x"]; % Scatter marker for correct and others
opts.ls = [":", "-", "-"]; % Line style for [ShortFP, MedFP, LongFP]
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
opts.plotsize2 = [9 3.5];
opts.plotsize3 = [3.5 3];
opts.plotsize4 = [2 3];

opts.sep_col = 0.5;
opts.sep_row = 1.5;

%%
fig = figure(23); clf(23);
set(gcf, 'unit', 'centimeters', 'position', [2 .7 39.5 25.4], 'paperpositionmode', 'auto', 'color', 'w');

mycolormap = customcolormap_preset("red-white-blue");
% mycolormap = "Turbo";

uicontrol('Style', 'text', 'parent', 23, 'units', 'normalized', 'position', [0.25 0.95 0.5 0.04],...
    'string', obj.Subject+" / "+obj.Task+" / "+obj.Sessions(1)+"_"+obj.Sessions(end), 'fontsize', 11, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set axes and plot
% Maze diagram
ha_diagram = axes;
set(ha_diagram, 'units', 'centimeters', 'position', [2+3*(opts.sep_row+opts.plotsize1(1))+5.5 1.35+4*(opts.sep_col+opts.plotsize1(2)), [3 6] ], 'nextplot', 'add', 'fontsize', 8);
plot_diagram(ha_diagram, opts);
set(ha_diagram, 'xlim', [0 6]);

% Hold duration scatter
ha_HD_scatter = axes;
set(ha_HD_scatter, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_scatter(ha_HD_scatter, obj, opts);
set(ha_HD_scatter, 'xtick', [], 'xlabel', []);

% Reaction time scatter
ha_RT_scatter = axes;
set(ha_RT_scatter, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_scatter(ha_RT_scatter, obj, opts);
set(ha_RT_scatter, 'xtick', [], 'xlabel', []);

% Performance track
ha_perf_track = axes;
set(ha_perf_track, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+4*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_performance_track_progress(ha_perf_track, obj, opts);
set(ha_perf_track, 'xtick', [], 'xlabel', []);

% Performance progress
ha_perf_progress = axes;
set(ha_perf_progress, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+4*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_performance_progress(ha_perf_progress, obj, opts);
set(ha_perf_progress, 'xtick', [], 'xlabel', [], 'yticklabel', [], 'ylabel', []);

% Reaction time violin plot
ha_RT_violin = axes;
set(ha_RT_violin, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_violin(ha_RT_violin, obj, opts);
set(ha_RT_violin, 'xtick', [], 'xlabel', [], 'yticklabel', [], 'ylabel', []);

%
ha_RT_median = axes;
set(ha_RT_median, 'units', 'centimeters', 'position', [.5+3*(opts.sep_row+opts.plotsize1(1))+1.5 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_median(ha_RT_median, obj, opts);
set(ha_RT_median, 'xtick', [], 'xlabel', []);

%
ha_HD_IQR = axes;
set(ha_HD_IQR, 'units', 'centimeters', 'position', [.5+3*(opts.sep_row+opts.plotsize1(1))+1.5 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_iqr(ha_HD_IQR, obj, opts);

%
ha_HD_heatmap_R = axes;
set(ha_HD_heatmap_R, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1)) 2+4*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_heatmap(ha_HD_heatmap_R, obj, "R", opts)
set(ha_HD_heatmap_R, 'xtick', [], 'xlabel', []);
clim(ha_HD_heatmap_R, [0, 10]);

colormap(ha_HD_heatmap_R, "Turbo");
cb9 = colorbar(ha_HD_heatmap_R, 'units', 'centimeters', 'position', [.5+3*(opts.sep_row+opts.plotsize1(1))-opts.sep_row+0.2 2+4*(opts.sep_col+opts.plotsize1(2)), [.3 opts.plotsize1(2)]]);
% cb9.Limits = [0 cb9.Limits(2)];
cb9.Limits = [0 10];
cb9.Label.String = "Prob. density (1/s)";
cb9.FontSize = 8;

%
ha_HD_heatmap_L = axes;
set(ha_HD_heatmap_L, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_heatmap(ha_HD_heatmap_L, obj, "L", opts)
set(ha_HD_heatmap_L, 'xtick', [], 'xlabel', []);
clim(ha_HD_heatmap_L, [0, 10]);

colormap(ha_HD_heatmap_L, "Turbo");
cb10 = colorbar(ha_HD_heatmap_L, 'units', 'centimeters', 'position', [.5+3*(opts.sep_row+opts.plotsize1(1))-opts.sep_row+0.2 2+3*(opts.sep_col+opts.plotsize1(2)), [.3 opts.plotsize1(2)]]);
% cb10.Limits = [0 cb10.Limits(2)];
cb10.Limits = [0 10];
cb10.Label.String = "Prob. density (1/s)";
cb10.FontSize = 8;

%
ha_HD_heatmap_diff = axes;
set(ha_HD_heatmap_diff, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1)) 2+2*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_heatmap(ha_HD_heatmap_diff, obj, "L-R", opts);
clim(ha_HD_heatmap_diff, [-8 8]);

colormap(ha_HD_heatmap_diff, mycolormap);
cb11 = colorbar(ha_HD_heatmap_diff, 'units', 'centimeters', 'position', [.5+3*(opts.sep_row+opts.plotsize1(1))-opts.sep_row+0.2 2+2*(opts.sep_col+opts.plotsize1(2)), [.3 4]]);
cb11.Label.String = "Δ Prob. density (1/s)";
cb11.FontSize = 8;

%
ha_HD_heatmap_colorbar = axes;
set(ha_HD_heatmap_colorbar, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_pdf_early_late(ha_HD_heatmap_colorbar, obj, "R", opts);

ha122 = axes;
set(ha122, 'units', 'centimeters', 'position', [.5+2*(opts.sep_row+opts.plotsize1(1)) 1.9+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_pdf_early_late(ha122, obj, "L", opts);

density_max = max([ha_HD_heatmap_colorbar.YLim(2) ha122.YLim(2)]);
set(ha_HD_heatmap_colorbar, 'ylim', [0 density_max], 'xticklabel', [], 'xlabel', []);
set(ha122, 'ylim', [0 density_max]);

%
ha_HD_CDF_R = axes;
set(ha_HD_CDF_R, 'units', 'centimeters', 'position', [1+3*(opts.sep_row+opts.plotsize1(1)) 1.5+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_cdf_early_late(ha_HD_CDF_R, obj, "R", opts);
set(ha_HD_CDF_R, 'xticklabel', [], 'xlabel', []);

ha_HD_CDF_L = axes;
set(ha_HD_CDF_L, 'units', 'centimeters', 'position', [1+3*(opts.sep_row+opts.plotsize1(1)) 1.9+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_cdf_early_late(ha_HD_CDF_L, obj, "L", opts);

ha_HD_legend = axes;
set(ha_HD_legend, 'units', 'centimeters', 'position', [.5+3*(opts.sep_row+opts.plotsize1(1))+5*opts.plotsize2(1)/6 1.5+0*(opts.sep_col+opts.plotsize1(2))+1*opts.plotsize2(2)/6, opts.plotsize2 ./ [5 3]], ...
    'nextplot', 'add', 'fontsize', 8, 'color', 'none', 'xcolor', 'none', 'ycolor', 'none', 'xlim', [0 3], 'ylim', [0 3]);
line(ha_HD_legend, [1 2], [2 2], 'Color', opts.color.PhaseEarly, 'LineWidth', 1, 'LineStyle', '-');
line(ha_HD_legend, [1 2], [1 1], 'Color', opts.color.PhaseLate, 'LineWidth', 1, 'LineStyle', '-');
text(ha_HD_legend, 2.2, 2, "Early", 'FontSize', 8);
text(ha_HD_legend, 2.2, 1, "Late", 'FontSize', 8);

%
ha14 = axes;
set(ha14, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time_scatter(ha14, obj, opts);
set(ha14, 'xtick', [], 'xlabel', []);

%
ha15 = axes;
set(ha15, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize1(1)) 2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_shuttle_time_scatter(ha15, obj, opts);

%
ha16 = axes;
set(ha16, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+3*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_violin(ha16, obj, opts)
set(ha16, 'xtick', [], 'xlabel', [], 'yticklabel', [], 'ylabel', []);

ha17 = axes;
set(ha17, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+1*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time_violin(ha17, obj, opts)
set(ha17, 'xtick', [], 'xlabel', [], 'yticklabel', [], 'ylabel', []);

ha18 = axes;
set(ha18, 'units', 'centimeters', 'position', [.5+1*(opts.sep_row+opts.plotsize1(1)) 2+0*(opts.sep_col+opts.plotsize1(2)), opts.plotsize1], 'nextplot', 'add', 'fontsize', 8);
plot_shuttle_time_violin(ha18, obj, opts)
set(ha18, 'yticklabel', [], 'ylabel', [])

ha19 = axes;
set(ha19, 'units', 'centimeters', 'position', [1.5+3*(opts.sep_row+opts.plotsize1(1)) 2.5+4*(opts.sep_col+opts.plotsize1(2)), opts.plotsize2 .* [3/5 1]], 'nextplot', 'add', 'fontsize', 8);
plot_interruption_early_late(ha19, obj, opts)

%% ha1. Make a diagram of the setup
    function plot_diagram(ax, opts)

        x = [1 3 5 8 9   11  11  9   8 5 3 1 1];
        y = [0 0 2 2 1.2 1.2 4.8 4.8 4 4 6 6 0];

        patch(ax, 'XData', y, 'YData', x, 'FaceColor', 'none', 'EdgeColor', 'k', 'linewidth', 2);

        viscircles(ax, [3, 9.6], 0.3, 'color', 'k', 'LineWidth', 1);
        text(ax, 3, 8.8, 'Init', 'FontWeight','bold', 'HorizontalAlignment', 'center');
        viscircles(ax, [3, 3.5], 0.3,  'color', 'k', 'LineWidth', 1);
        text(ax, 3, 4.5, 'Cent', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

        viscircles(ax, [1.5, 3], 0.3, 'color', opts.color.PortR);
        text(ax, 1.5, 2.2, 'Right', 'FontWeight','bold', 'HorizontalAlignment', 'center', 'Color', opts.color.PortR);
        viscircles(ax, [4.5, 3], 0.3,  'color', opts.color.PortL);
        text(ax, 4.5, 2.2, 'Left', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', opts.color.PortL);

        set(ax, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ylim', [0 12], 'xlim', [0 6]);
    end

%% ha2. Plot every trial's hold duration during training
    function plot_hold_duration_scatter(ax, obj, opts)

        for fp_this = 1:length(obj.MixedFP)
            yline(ax, obj.MixedFP(fp_this), 'Color', [.8 .8 .8], 'LineWidth', opts.lw(fp_this), 'LineStyle', '-', 'Alpha', 0.4);
        end
        xline(ax, opts.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

        s_sep = [0 opts.session_sep opts.trial_ticks(end)];
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [s_sep(dcz_ind); s_sep(dcz_ind); s_sep(dcz_ind+1); s_sep(dcz_ind+1)], 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.trial_ticks(obj.Ind.wrongL), ...
            obj.BehavTable.HoldDuration(obj.Ind.wrongL), ...
            22*obj.BehavTable.FP(obj.Ind.wrongL), opts.color.PortL, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 wrong
        scatter(ax, opts.trial_ticks(obj.Ind.wrongR), ...
            obj.BehavTable.HoldDuration(obj.Ind.wrongR), ...
            22*obj.BehavTable.FP(obj.Ind.wrongR), opts.color.PortR, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % premature
        scatter(ax, opts.trial_ticks(obj.Ind.premature), ...
            obj.BehavTable.HoldDuration(obj.Ind.premature), ...
            22*obj.BehavTable.FP(obj.Ind.premature), opts.color.Premature, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % late
        scatter(ax, opts.trial_ticks(obj.Ind.late), ...
            obj.BehavTable.HoldDuration(obj.Ind.late), ...
            22*obj.BehavTable.FP(obj.Ind.late), opts.color.Late, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.trial_ticks(obj.Ind.correctL), ...
            obj.BehavTable.HoldDuration(obj.Ind.correctL), ...
            18*obj.BehavTable.FP(obj.Ind.correctL), opts.color.PortL, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 correct
        scatter(ax, opts.trial_ticks(obj.Ind.correctR), ...
            obj.BehavTable.HoldDuration(obj.Ind.correctR), ...
            18*obj.BehavTable.FP(obj.Ind.correctR), opts.color.PortR, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Hold duration (s)';
        set(ax, 'xlim', [0 opts.trial_ticks(end)+5], 'ylim', [0 2.5], 'ticklength', [0.01 0.1]);

    end

%% ha3. Plot correct trial's reaction during training
    function plot_reaction_time_scatter(ax, obj, opts)

        xline(ax, opts.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-')
        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');

        s_sep = [0 opts.session_sep opts.trial_ticks(end)];
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [s_sep(dcz_ind); s_sep(dcz_ind); s_sep(dcz_ind+1); s_sep(dcz_ind+1)], 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.trial_ticks(obj.Ind.wrongL), ...
            obj.BehavTable.RT(obj.Ind.wrongL), ...
            22*obj.BehavTable.FP(obj.Ind.wrongL), opts.color.PortL, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 wrong
        scatter(ax, opts.trial_ticks(obj.Ind.wrongR), ...
            obj.BehavTable.RT(obj.Ind.wrongR), ...
            22*obj.BehavTable.FP(obj.Ind.wrongR), opts.color.PortR, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.trial_ticks(obj.Ind.correctL), ...
            obj.BehavTable.RT(obj.Ind.correctL), ...
            18*obj.BehavTable.FP(obj.Ind.correctL), opts.color.PortL, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 correct
        scatter(ax, opts.trial_ticks(obj.Ind.correctR), ...
            obj.BehavTable.RT(obj.Ind.correctR), ...
            18*obj.BehavTable.FP(obj.Ind.correctR), opts.color.PortR, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Reaction time (s)';
        set(ax, 'xlim', [0 opts.trial_ticks(end)+5], 'ylim', [0 obj.Bins.RT(end)], 'ticklength', [0.01 0.1]);
    end

%% ha3. Plot correct trial's reaction during training
    function plot_movement_time_scatter(ax, obj, opts)

        xline(ax, opts.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-')

        s_sep = [0 opts.session_sep opts.trial_ticks(end)];
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [s_sep(dcz_ind); s_sep(dcz_ind); s_sep(dcz_ind+1); s_sep(dcz_ind+1)], 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.trial_ticks(obj.Ind.wrongL), ...
            obj.BehavTable.MovementTime(obj.Ind.wrongL), ...
            22*obj.BehavTable.FP(obj.Ind.wrongL), opts.color.PortL, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 wrong
        scatter(ax, opts.trial_ticks(obj.Ind.wrongR), ...
            obj.BehavTable.MovementTime(obj.Ind.wrongR), ...
            22*obj.BehavTable.FP(obj.Ind.wrongR), opts.color.PortR, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 1 correct (defined by their action, which is also the target for correct response)
        scatter(ax, opts.trial_ticks(obj.Ind.correctL), ...
            obj.BehavTable.MovementTime(obj.Ind.correctL), ...
            18*obj.BehavTable.FP(obj.Ind.correctL), opts.color.PortL, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        % port 2 correct
        scatter(ax, opts.trial_ticks(obj.Ind.correctR), ...
            obj.BehavTable.MovementTime(obj.Ind.correctR), ...
            18*obj.BehavTable.FP(obj.Ind.correctR), opts.color.PortR, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Movement time (s)';
        set(ax, 'xlim', [0 opts.trial_ticks(end)+5], 'ylim', [0 obj.Bins.MovementTime(end)], 'ticklength', [0.01 0.1]);
    end

%% ha3. Plot correct trial's reaction during training
    function plot_shuttle_time_scatter(ax, obj, opts)

        xline(ax, opts.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-')

        s_sep = [0 opts.session_sep opts.trial_ticks(end)];
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [s_sep(dcz_ind); s_sep(dcz_ind); s_sep(dcz_ind+1); s_sep(dcz_ind+1)], 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % port 1 wrong (defined by their action, which is different from target)
        scatter(ax, opts.trial_ticks, ...
            log10(obj.BehavTable.ShuttleTime), ...
            18, [0 0 0], opts.mk(1), 'MarkerEdgeAlpha', 0.4, 'linewidth', 1);

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Log Shuttle time (s)';
        set(ax, 'xlim', [0 opts.trial_ticks(end)+5], 'ylim', [0 obj.Bins.ShuttleTimeLog(end)], 'ticklength', [0.01 0.1]);
    end

%% ha4. Plot performance track of each session
    function plot_performance_track_progress(ax, obj, opts)

        xline(ax, opts.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-')

        s_sep = [0 opts.session_sep opts.trial_ticks(end)];
        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [s_sep(dcz_ind); s_sep(dcz_ind); s_sep(dcz_ind+1); s_sep(dcz_ind+1)], 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for s_this = 1:obj.NumSessions
            ind_this = find(obj.PerformanceTrack.Sessions==obj.Sessions(s_this));
            plot(ax, obj.PerformanceTrack.WinPosProgress(ind_this), obj.PerformanceTrack.CorrectRatio(ind_this),  'linestyle', '-', 'color', opts.color.Correct, ...
                'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Correct,   'markeredgecolor', 'w');
            plot(ax, obj.PerformanceTrack.WinPosProgress(ind_this), obj.PerformanceTrack.WrongRatio(ind_this),    'linestyle', '-', 'color', opts.color.Wrong, ...
                'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Wrong,     'markeredgecolor', 'w');
            plot(ax, obj.PerformanceTrack.WinPosProgress(ind_this), obj.PerformanceTrack.PrematureRatio(ind_this), 'linestyle', '-', 'color', opts.color.Premature, ...
                'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
            plot(ax, obj.PerformanceTrack.WinPosProgress(ind_this), obj.PerformanceTrack.LateRatio(ind_this),      'linestyle', '-', 'color', opts.color.Late, ...
                'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Late,      'markeredgecolor', 'w');
        end

        ax.XLabel.String = 'Time in training (s)';
        ax.YLabel.String = 'Performance (%)';
        set(ax, 'XLim', [0 opts.trial_ticks(end)+5], 'YLim', [0 100]);
    end

%% ha5. Plot performance progress
    function plot_performance_progress(ax, obj, opts)

        session_id = 1:obj.NumSessions;
        xline(ax, session_id-0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');
        xline(ax, session_id, 'Color', [.8 .8 .8], 'LineWidth', 1, 'LineStyle', ':');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for fp_this = 1:length(obj.MixedFP)
            for p_this = 1:length(obj.Ports)

                ind_this = find(obj.Performance.Foreperiod==obj.MixedFP(fp_this) & obj.Performance.TargetPort==obj.Ports(p_this));
                
                scatter(ax, session_id+(-1.5+p_this)/3, obj.Performance.PrematureRatio(ind_this), ...
                    9*fp_this, 'Marker', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.Premature, 'MarkerFaceAlpha', .3);
                scatter(ax, session_id+(-1.5+p_this)/3, obj.Performance.LateRatio(ind_this), ...
                    9*fp_this, 'Marker', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.Late, 'MarkerFaceAlpha', .3);
                scatter(ax, session_id+(-1.5+p_this)/3, obj.Performance.WrongRatio(ind_this), ...
                    9*fp_this, 'Marker', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.Wrong, 'MarkerFaceAlpha', .3);
                scatter(ax, session_id+(-1.5+p_this)/3, obj.Performance.CorrectRatio(ind_this), ...
                    9*fp_this, 'Marker', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.Correct, 'MarkerFaceAlpha', .5);
            end
        end

        ind_this = find(obj.Performance.Foreperiod==0 & obj.Performance.TargetPort=="Both");

        plot(ax, session_id, obj.Performance.PrematureRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Premature, ...
            'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.LateRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Late, ...
            'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Late, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.WrongRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Wrong, ...
            'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Wrong, 'markeredgecolor', 'w');
        plot(ax, session_id, obj.Performance.CorrectRatio(ind_this), 'o', 'linestyle', '-', 'color', opts.color.Correct, ...
            'markersize', 5, 'linewidth', 1, 'markerfacecolor', opts.color.Correct, 'markeredgecolor', 'w');

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Performance (%)';
        set(ax, 'xlim', [0.5 obj.NumSessions+0.5], 'xtick', session_id, 'xticklabel', opts.session_date);
    end

%% ha6. Plot hold duration of each session (violin), each port
    function plot_hold_duration_violin(ax, obj, opts)

        for fp_this = 1:length(obj.MixedFP)
            yline(ax, obj.MixedFP(fp_this), 'Color', [.8 .8 .8], 'LineWidth', opts.lw(fp_this), 'LineStyle', '-', 'Alpha', 0.4);
        end
        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

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
                ind_this = session_id & obj.BehavTable.Stage==1 & eval("obj.Ind.port"+obj.Ports(p_this));
                thisHD(ind_this, 2*(s_this-1)+p_this) = HD(ind_this);
            end
        end
        thisHD(1, all(isnan(thisHD))) = 0;

        violinplot({thisHD(:, 2:2:end), thisHD(:, 1:2:end)}, obj.Sessions, ...
            'ViolinColor', {repmat(opts.color.PortR, obj.NumSessions, 1), repmat(opts.color.PortL, obj.NumSessions, 1)}, ...
            'ScatterSize', 4, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false, 'Bandwidth', 0.1);

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Hold duration (s)';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 obj.Bins.HoldDuration(end)], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%% ha6. Plot reaction time of each session (violin), each port
    function plot_reaction_time_violin(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');
        yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        thisRT = nan(height(obj.BehavTable), length(obj.Ports)*obj.NumSessions);

        RT = obj.BehavTable.RT;
        [~, ~, indrmv] = rmoutliers_custome(RT);
        RT(indrmv) = nan;

        for s_this = 1:obj.NumSessions
            session_id = obj.BehavTable.SessionDate==obj.Sessions(s_this);
            for p_this = 1:length(obj.Ports)
                ind_this = session_id & obj.BehavTable.Stage==1 & eval("obj.Ind.correct"+obj.Ports(p_this));
                thisRT(ind_this, 2*(s_this-1)+p_this) = RT(ind_this);
            end
        end

        thisRT(1, all(isnan(thisRT))) = 0;
        violinplot({thisRT(:, 2:2:end), thisRT(:, 1:2:end)}, obj.Sessions, ...
            'ViolinColor', {repmat(opts.color.PortR, obj.NumSessions, 1), repmat(opts.color.PortL, obj.NumSessions, 1)}, ...
            'ScatterSize', 4, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false);
        
        session_order = [1:obj.NumSessions; 1:obj.NumSessions];
        session_order = session_order(:);
        scatter(ax, session_order, median(thisRT, 'omitnan'), 24, ...
            repmat([opts.color.PortL; opts.color.PortR], obj.NumSessions, 1), ...
            'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3]);

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Reaction time (s)';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 obj.Bins.RT(end)], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%% ha6. Plot reaction time of each session (violin), each port
    function plot_movement_time_violin(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        thisMT = nan(height(obj.BehavTable), length(obj.Ports)*obj.NumSessions);

        MT = obj.BehavTable.MovementTime;
        [~, ~, indrmv] = rmoutliers_custome(MT);
        MT(indrmv) = nan;

        for s_this = 1:obj.NumSessions
            session_id = obj.BehavTable.SessionDate==obj.Sessions(s_this);
            for p_this = 1:length(obj.Ports)
                ind_this = session_id & obj.BehavTable.Stage==1 & eval("obj.Ind.correct"+obj.Ports(p_this));
                thisMT(ind_this, 2*(s_this-1)+p_this) = MT(ind_this);
            end
        end

        thisMT(1, all(isnan(thisMT))) = 0;
        violinplot({thisMT(:, 2:2:end), thisMT(:, 1:2:end)}, obj.Sessions, ...
            'ViolinColor', {repmat(opts.color.PortR, obj.NumSessions, 1), repmat(opts.color.PortL, obj.NumSessions, 1)}, ...
            'ScatterSize', 4, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false);

        session_order = [1:obj.NumSessions; 1:obj.NumSessions];
        session_order = session_order(:);
        scatter(ax, session_order, median(thisMT, 'omitnan'), 24, ...
            repmat([opts.color.PortL; opts.color.PortR], obj.NumSessions, 1), ...
            'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3]);

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Movement time (s)';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 obj.Bins.MovementTime(end)], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%% ha6. Plot reaction time of each session (violin), each port
    function plot_shuttle_time_violin(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        thisSTLog = nan(height(obj.BehavTable), obj.NumSessions);

        STLog = log(obj.BehavTable.ShuttleTime);
        [~, ~, indrmv] = rmoutliers_custome(STLog);
        STLog(indrmv) = nan;

        for s_this = 1:obj.NumSessions
            session_id = obj.BehavTable.SessionDate==obj.Sessions(s_this);

            ind_this = session_id & obj.BehavTable.Stage==1;
            thisSTLog(ind_this, s_this) = STLog(ind_this);
        end

        thisSTLog(1, all(isnan(thisSTLog))) = 0;
        violinplot(thisSTLog, obj.Sessions, ...
            'ViolinColor', zeros(obj.NumSessions, 3), ...
            'ScatterSize', 4, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false);

        scatter(ax, 1:obj.NumSessions, median(thisSTLog, 'omitnan'), 24, ...
            repmat([1 1 1], obj.NumSessions, 1), ...
            'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3]);

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Log Shuttle time (s)';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 obj.Bins.ShuttleTimeLog(end)], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_reaction_time_median(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');
        yline(0.5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for s_this = 1:obj.NumSessions
            for p_this = 1:2
                plot([s_this-0.3 s_this s_this+0.3], ...
                    obj.RTStat.Median(obj.RTStat.Port==obj.Ports(p_this) & obj.RTStat.Sessions==obj.Sessions(s_this) & obj.RTStat.thisFP~=0), ...
                    'Color', [eval("opts.color.Port"+obj.Ports(p_this)) 0.6], 'LineWidth', 1.5);
                scatter([s_this-0.3 s_this s_this+0.3], ...
                    obj.RTStat.Median(obj.RTStat.Port==obj.Ports(p_this) & obj.RTStat.Sessions==obj.Sessions(s_this) & obj.RTStat.thisFP~=0), ...
                    6*(obj.MixedFP+0.5), repmat(eval("opts.color.Port"+obj.Ports(p_this)), 3, 1), 'filled', 'Marker', 'o', ...
                   'MarkerEdgeAlpha', 0.8, 'MarkerFaceAlpha', 0.8);
            end
        end

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Reaction time median (s)';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 .5], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_hold_duration_iqr(ax, obj, opts)

        xline(ax, (1:obj.NumSessions-1)+0.5, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

        dcz_ind = find(obj.Treatment == "DCZ");
        num_dcz = length(dcz_ind);
        patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
            'FaceColor', opts.color.Treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        for s_this = 1:obj.NumSessions
            for p_this = 1:length(obj.Ports)
                plot([s_this-0.3 s_this s_this+0.3], ...
                    obj.HDStat.IQR(obj.HDStat.Port==obj.Ports(p_this) & obj.HDStat.Sessions==obj.Sessions(s_this) & obj.HDStat.thisFP~=0), ...
                    'Color', [eval("opts.color.Port"+obj.Ports(p_this)) 0.6], 'LineWidth', 1.5);
                scatter([s_this-0.3 s_this s_this+0.3], ...
                    obj.HDStat.IQR(obj.HDStat.Port==obj.Ports(p_this) & obj.HDStat.Sessions==obj.Sessions(s_this) & obj.HDStat.thisFP~=0), ...
                    6*(obj.MixedFP+0.5), repmat(eval("opts.color.Port"+obj.Ports(p_this)), 3, 1), 'filled', 'Marker', 'o', ...
                   'MarkerEdgeAlpha', 0.8, 'MarkerFaceAlpha', 0.8);
            end
        end

        ax.XLabel.String = 'Sessions';
        ax.YLabel.String = 'Hold duration IQR (s)';
        set(ax, 'xlim', [.5 obj.NumSessions+.5], 'ylim', [0 max(obj.HDStat.IQR)*1.5], 'xtick', 1:obj.NumSessions, ...
            'xticklabel', opts.session_date, 'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
    function plot_hold_duration_heatmap(ax, obj, port, opts)

        hd_l = zeros(length(obj.Bins.HoldDuration), obj.NumSessions*length(obj.MixedFP));
        hd_r = zeros(length(obj.Bins.HoldDuration), obj.NumSessions*length(obj.MixedFP));

        for s_this = 1:obj.NumSessions
            for fp_this = 1:length(obj.MixedFP)
                hd_l(:, s_this + obj.NumSessions*(fp_this-1)) = obj.HDPDF{s_this}{fp_this, 1};
                hd_r(:, s_this + obj.NumSessions*(fp_this-1)) = obj.HDPDF{s_this}{fp_this, 2};
            end
        end

        switch port
            case {"L", "l"}
                imagesc(ax, hd_l);
                ax.YLabel.String = 'Hold duration Left (s)';
                clim(ax, [-max(max([hd_l; hd_r])) max(max([hd_l; hd_r]))]);
            case {"R", "r"}
                imagesc(ax, hd_r);
                ax.YLabel.String = 'Hold duration Right (s)';
                clim(ax, [-max(max([hd_l; hd_r])) max(max([hd_l; hd_r]))]);
            case {"Diff", "diff", "L-R", "l-r"}
                imagesc(ax, hd_l-hd_r);
                ax.YLabel.String = 'Hold duration L - R (s)';
                if abs(min(min(hd_l-hd_r))) < abs(max(max(hd_l-hd_r)))
                    clim(ax, [-abs(max(max(hd_l-hd_r))) abs(max(max(hd_l-hd_r)))]);
                else
                    clim(ax, [-abs(min(min(hd_l-hd_r))) abs(min(min(hd_l-hd_r)))]);
                end
        end

        line(ax, [.5 obj.NumSessions+.5], [obj.MixedFP(1) obj.MixedFP(1)]/obj.Bins.width, 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', ':');
        line(ax, [obj.NumSessions+.5 2*obj.NumSessions+.5], [obj.MixedFP(2) obj.MixedFP(2)]/obj.Bins.width, 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', ':');
        line(ax, [2*obj.NumSessions+.5 3*obj.NumSessions+.5], [obj.MixedFP(3) obj.MixedFP(3)]/obj.Bins.width, 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', ':');

        xline(ax, (1:2)*obj.NumSessions + .5, 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--');

        for s_this = 1:obj.NumSessions
            if obj.Label(s_this) == "Chemo"
                for fp_this = 1:length(obj.MixedFP)
                    patch(ax, 'XData', s_this + obj.NumSessions*(fp_this-1)+[-0.5 0.5 0.5 -0.5], 'YData', [-.5 -.5 length(obj.Bins.HoldDuration)+.5 length(obj.Bins.HoldDuration)+.5], 'FaceColor', 'none', ...
                        'EdgeColor', opts.color.Treat, 'LineWidth', 0.5, 'LineStyle', ':', 'EdgeAlpha', 0.8);
                end
            end
        end

        ax.XLabel.String = 'Sessions';
        set(ax, 'xlim', [.5 3*obj.NumSessions+.5], 'ylim', [-.5 length(obj.Bins.HoldDuration)+.5], ...
            'xtick', 1:3*obj.NumSessions, 'xticklabel', repmat(opts.session_date, 3, 1), ...
            'ytick', 0:250:length(obj.Bins.HoldDuration), 'yticklabel', string(0:0.5:obj.Bins.HoldDuration(end)), ...
            'ticklength', [0.01 0.1], 'box', 'off');
    end

%%
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
            case {"R"}
                ax.YLabel.String = "Prob. density Right (1/s)";
        end
        
        ax.XLabel.String = "Hold duration (s)";
        set(ax, 'xlim', [0 obj.Bins.HoldDuration(end)], 'ylimmode', 'auto');
    end

    function plot_hold_duration_cdf_early_late(ax, obj, port, opts)

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

            hd_early_cdf = ksdensity(hd_this(1:obj.PhaseCount), obj.Bins.HoldDuration, 'Function', 'cdf');
            hd_late_cdf  = ksdensity(hd_this(end-obj.PhaseCount+1:end), obj.Bins.HoldDuration, 'Function', 'cdf');

            plot(ax, obj.Bins.HoldDuration, hd_early_cdf, ...
                'color', opts.color.PhaseEarly, 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
            plot(ax, obj.Bins.HoldDuration, hd_late_cdf, ...
                'color', opts.color.PhaseLate, 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
        end

        switch upper(string(port))
            case {"L"}
                ax.YLabel.String = "Cum. density Left";
            case {"R"}
                ax.YLabel.String = "Cum. density Right";
        end
        
        ax.XLabel.String = "Hold duration (s)";
        set(ax, 'xlim', [0 obj.Bins.HoldDuration(end)], 'ylim', [0 1]);
    end

    function plot_interruption_early_late(ax, obj, opts)

        for fp_this = 1:length(obj.MixedFP)
            xline(ax, obj.MixedFP(fp_this), 'color', [.7 .7 .7], 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
        end

        trial_ind = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="None")) & obj.BehavTable.Stage==1;

        if sum(trial_ind) < 6*obj.PhaseCount*2
            inter = [obj.InterruptionNone; obj.InterruptionSaline];
            inter = sortrows(inter, "TrialProgress");
            trial_ind = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="None" | obj.Treatment=="Saline")) & obj.BehavTable.Stage==1;
        else
            inter = obj.InterruptionNone;
        end

        if sum(trial_ind) < 6*obj.PhaseCount*2
            return;
        end

        behav_this = obj.BehavTable(trial_ind, :);

        behav_early = behav_this(1:6*obj.PhaseCount, :);
        behav_late  = behav_this(end-6*obj.PhaseCount+1:end, :);
        trial_id_early = behav_early.TrialProgress(end);
        trial_id_late  = behav_late.TrialProgress(1);

        nums_early = zeros(1, length(obj.Bins.Interruption)-1);
        nums_late  = zeros(1, length(obj.Bins.Interruption)-1);
        for i = 1:length(nums_early)
            nums_early(i) = sum(behav_early.HoldDuration >= obj.Bins.Interruption(i));
            nums_late(i)  = sum(behav_late.HoldDuration  >= obj.Bins.Interruption(i));
        end

        h_early = histcounts(inter.On(inter.TrialProgress<=trial_id_early), "BinEdges", obj.Bins.Interruption) ./ nums_early;
        bar(ax, obj.Bins.Interruption(1:end-1)+.5*obj.Bins.widthInter, h_early, 1, 'FaceColor', opts.color.PhaseEarly , 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'EdgeAlpha', 0.6);
        h_late  = histcounts(inter.On(inter.TrialProgress>=trial_id_late) , "BinEdges", obj.Bins.Interruption) ./ nums_late;
        bar(ax, obj.Bins.Interruption(1:end-1)+.5*obj.Bins.widthInter, h_late , 1, 'FaceColor', opts.color.PhaseLate  , 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'EdgeAlpha', 0.6);

        y_control = smoothdata(h_early, 'gaussian', 0.18*(length(h_early)));
        y_chemo   = smoothdata(h_late , 'gaussian', 0.18*(length(h_late)));

        plot(ax, obj.Bins.Interruption(1:end-1)+.5*obj.Bins.widthInter, y_control, 'Color', opts.color.PhaseEarly, 'LineStyle', '-', 'LineWidth', 1.5);
        plot(ax, obj.Bins.Interruption(1:end-1)+.5*obj.Bins.widthInter, y_chemo  , 'Color', opts.color.PhaseLate , 'LineStyle', '-', 'LineWidth', 1.5);

        ax.YLabel.String = "Interruptions / trial";
        ax.XLabel.String = "Hold duration (s)";
        set(ax, 'xlim', [0 1.5]);
    end

end
