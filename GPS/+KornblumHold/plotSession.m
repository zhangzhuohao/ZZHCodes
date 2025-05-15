function fig_session = plotSession(obj)

ax_sz_1 = [5 2];
ax_sz_2 = [2.4 2];
ax_sz_3 = [1.5 2];

if contains(obj.Protocol, "Self")
    ls = "-";
    lw = 1;
    sort_id = 2;
    rw = .8;
    ax_sz_s = ax_sz_1;
    cue_uncue_code = 0;
    cue_uncue_label = "Uncued";
else
    ls = ["-", ":"];
    lw = [1 1.25];
    sort_id = [1 2];
    rw = [.5 .8];
    ax_sz_s = ax_sz_2;
    cue_uncue_code = [0 1];
    cue_uncue_label = ["Uncued" "Cued"];
end

fig_session = figure(11); clf(fig_session);
set(fig_session, 'Visible', 'on', 'Units', 'centimeters', 'Position', [5 5 17 10.5], 'Color', 'w', 'toolbar', 'none');

fig_title = sprintf("%s / %s / %s / %s", obj.Subject, obj.Session, obj.Protocol, obj.Label);
set_fig_title(fig_session, fig_title);

% cumulative plot
ax_cumulative = axes(fig_session, "Units", "centimeters", "Position", [1.2 7.5 ax_sz_1], 'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');
stairs(ax_cumulative, obj.TrialCentInTime, obj.Trials, 'LineWidth', 1, 'LineStyle', '-', 'Color', [0 0 0]);
set(ax_cumulative, 'XLim', [0 1.05*max(obj.TrialCentInTime)])
xlabel(ax_cumulative, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax_cumulative, 'Trial #', 'FontWeight', 'bold');

% performance track
draw_perf_track(obj, fig_session, [1.2 4 ax_sz_1], ax_sz_2, obj.PerformanceTrack(sort_id,:), ls, lw);

% performance ratio
draw_perf_bars(obj, fig_session, [1.2 1 ax_sz_1], ax_sz_2, cue_uncue_code, cue_uncue_label);

% shuttle time
% scatter plot
ax_st_s = draw_scatter(obj, fig_session, [7.5 7.5 ax_sz_1], ax_sz_1, ...
    {log10(obj.ST)}, {obj.TrialCentInTime}, ...
    {repmat("Correct", obj.NumTrials, 1)}, {[0 0 0]}, [], []);
for i = 1:length(ax_st_s)
    set(ax_st_s{i}, 'YLim', [log10(0.5) 1], 'YTick', log10([0.1:0.1:1 2:10]), 'YTickLabel', ["0.1" repmat("",1,8) "1", repmat("", 1, 8), "10"]);
end
ylabel(ax_st_s{1}, 'ST (s)', 'FontWeight', 'bold');

% prob. density
c_first = .4 * ones(1, 3);
c_last  = 0  * ones(1, 3);
ax_st_d = draw_density(obj, fig_session, [13 7.5 3.5 2], ax_sz_3, ...
    obj.LogSTPDF([1 end])', [], [], ...
    {{c_first}, {c_last}}, {'-', '-'}, {1, 1});
for i = 1:length(ax_st_d)
    set(ax_st_d{i}, 'YLim', [log10(0.5)  1], 'YTick', log10([0.1:0.1:1 2:10]), 'YTickLabel', []);
end
text(ax_st_d{1}, ax_st_d{1}.XLim(2), .9, 'First 1/3', 'Color', c_first, 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(ax_st_d{2}, ax_st_d{2}.XLim(2), .9, 'Last 1/3' , 'Color', c_last , 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

% hold duration
% scatter plot
ax_hd_s = draw_scatter(obj, fig_session, [7.5 4 ax_sz_1], ax_sz_s, ...
    obj.HDSorted(sort_id, :), obj.TimeSorted(sort_id, :), obj.OutcomeSorted(sort_id, :), ...
    {GPSColor.PortL, GPSColor.PortR}, obj.TargetFP, rw);
for i = 1:length(ax_hd_s)
    set(ax_hd_s{i}, 'YLim', [-1 1]+obj.TargetFP);
    xlabel(ax_hd_s{i}, []);
end
ylabel(ax_hd_s{1}, 'HD (s)', 'FontWeight', 'bold');

% prob. density
ax_hd_d = draw_density(obj, fig_session, [13 4 3.5 2], ax_sz_3, ...
    obj.HDPDF(sort_id, :), obj.TargetFP, rw, ...
    {{GPSColor.PortL}, {GPSColor.PortR}}, {ls, ls}, {lw, lw});

for i = 1:length(ax_hd_d)
    set(ax_hd_d{i}, 'YLim', [-1 1]+obj.TargetFP, 'YTickLabel', []);
    xlabel(ax_hd_d{i}, []);
end

% movement time
ind = obj.Stage==1 & obj.Outcome=="Correct";
time_this = obj.TrialCentInTime(ind);
outcome_this = obj.Outcome(ind);
refs_this = cellfun(@(x) x(ind), obj.SortRefs, 'UniformOutput', false);
time_sorted = obj.sort_data(time_this, refs_this, obj.SortCodes);
outcome_sorted = obj.sort_data(outcome_this, refs_this, obj.SortCodes);

% scatter plot
ax_mt_s = draw_scatter(obj, fig_session, [7.5 1 ax_sz_1], ax_sz_s, ...
    obj.MTSorted(sort_id,:), time_sorted(sort_id,:), outcome_sorted(sort_id,:), ...
    {GPSColor.PortL, GPSColor.PortR}, [], []);
for i = 1:length(ax_mt_s)
    title(ax_mt_s{i}, []);
    set(ax_mt_s{i}, 'YLim', [.2 1]);
end
ylabel(ax_mt_s{1}, 'MT (s)', 'FontWeight', 'bold');

% prob. density
ax_mt_d = draw_density(obj, fig_session, [13 1 3.5 2], ax_sz_3, ...
    obj.MTPDF(sort_id,:), [], [], ...
    {{GPSColor.PortL}, {GPSColor.PortR}}, {ls, ls}, {lw, lw});
for i = 1:length(ax_mt_d)
    set(ax_mt_d{i}, 'YLim', [.2 1], 'YTickLabel', []);
end

end

%% Functions
% performance track
function ax = draw_perf_track(obj, fig, pos, ax_sz, perf_track, ls, lw)
draw_sz = size(perf_track);
if draw_sz(2)==1
    draw_sz = [1 draw_sz(1)];
    perf_track = reshape(perf_track, draw_sz);
end
ax = cell(1, draw_sz(2));
[dist_w, ~] = obj.get_plot_dist(ax, pos, ax_sz);

x_lim = max(max(cellfun(@(x) max(x.Pos), perf_track)));
for j = 1:draw_sz(2)
    x_now = pos(1) + dist_w * (j-1);
    y_now = pos(2);

    ax{j} = axes(fig, "Units", "centimeters", "Position", [x_now y_now ax_sz], ...
        'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');
    for i = 1:draw_sz(1)
        for k = 1:length(obj.PerformanceType)
            perf_this = obj.PerformanceType(k);
            plot(ax{j}, perf_track{i,j}.Pos, 100*perf_track{i,j}.(perf_this), ls(i), 'LineWidth', lw(i), 'Color', GPSColor.(perf_this));
        end
        set(ax{j}, 'XLim', [0 1.05*x_lim], 'YLim', [0 100], 'YTick', [0 50 100]);

        if j==1 && i==1
            ylabel(ax{j}, 'Perfomance %', 'FontWeight', 'bold');
            xlabel(ax{j}, 'Time (s)', 'FontWeight', 'bold');
        elseif j~=1
            set(ax{j}, 'YTickLabel', []);
        end
        if i==1
            switch j
                case 1
                    title(ax{j}, 'Left', 'FontWeight', 'bold', 'Color', GPSColor.PortL);
                case 2
                    title(ax{j}, 'Right', 'FontWeight', 'bold', 'Color', GPSColor.PortR);
            end
        end
    end
end
drawnow();
end % draw_perf_track

function ax = draw_perf_bars(obj, fig, pos, ax_sz, code, lb)
ax = cell(1, 2);
[dist_w, ~] = obj.get_plot_dist(ax, pos, ax_sz);

for j = 1:2
    x_now = pos(1) + dist_w * (j-1);
    y_now = pos(2);

    ax{j} = axes(fig, "Units", "centimeters", "Position", [x_now y_now ax_sz], ...
        'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');
    id = zeros(1, length(code));
    for i = 1:length(code)
        id(i) = find(obj.Performance.Cued==code(i) & obj.Performance.PortCorrect==j);
    end
    b = bar(ax{j}, code, 100*table2array(obj.Performance(id, end-3:end)), 'EdgeColor', 'none');
    for i = 1:4
        b(i).FaceColor = GPSColor.(obj.PerformanceType(i));
    end
    set(ax{j}, 'XLim', code+[-.5 .5], 'XTick', code, 'XTickLabel', lb, 'XTickLabelRotation', 0, 'YLim', [0 100]);
    if j==1
        ylabel(ax{j}, 'Perfomance %', 'FontWeight', 'bold');
        xlabel(ax{j}, 'Condition', 'FontWeight', 'bold');
    elseif j~=1
        set(ax{j}, 'YTickLabel', []);
    end
end
drawnow();
end % draw_perf_bars

% scatter plot
function ax = draw_scatter(obj, fig, pos, ax_sz, data_sorted, time_sorted, perf_sorted, color, fp, rw)
draw_sz = size(data_sorted, 1);
ax = cell(1, draw_sz);
dist_w = obj.get_plot_dist(ax, pos, ax_sz);

x_lim = max(obj.TrialCentInTime);
for i = 1:draw_sz
    x_now = pos(1) + dist_w * (i-1);
    y_now = pos(2);

    ax{i} = axes(fig, "Units", "centimeters", "Position", [x_now y_now ax_sz], ...
        'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');

    switch length(data_sorted(1,:))
        case 1
            data_this = data_sorted{i,1};
            time_this = time_sorted{i,1};
            color_this = repmat(color{1}, length(data_sorted{i,1}), 1);
            perf_this = perf_sorted{i,1};
        case 2
            data_this = [data_sorted{i,1}; data_sorted{i,2}];
            time_this = [time_sorted{i,1}; time_sorted{i,2}];
            color_this = [repmat(color{1}, length(data_sorted{i,1}), 1); repmat(color{2}, length(data_sorted{i,2}), 1)];
            perf_this = [perf_sorted{i,1}; perf_sorted{i,2}];
    end

    [time_this, id_rank] = sort(time_this);
    data_this = data_this(id_rank);
    color_this = color_this(id_rank, :);
    perf_this = perf_this(id_rank);

    mk_this = repmat("x", length(data_this), 1);
    mk_this(perf_this=="Correct") = "o";

    for j = 1:length(data_this)
        scatter(ax{i}, time_this(j), data_this(j), 16, color_this(j,:), 'Marker', mk_this(j), 'LineWidth', 1);
    end
    if ~isempty(fp)
        yline(ax{i}, fp, 'LineWidth', .5, 'LineStyle', '-');
        if ~isempty(rw)
            yline(ax{i}, fp+rw(i), 'LineWidth', .5, 'LineStyle', ':');
        end
    end
    set(ax{i}, 'XLim', [0 1.05*x_lim]);

    if i==1
        xlabel(ax{i}, 'Time (s)', 'FontWeight', 'bold');
    else
        set(ax{i}, 'YTickLabel', []);
    end
    if draw_sz>1
        switch i
            case 1
                title(ax{i}, 'Cued', 'FontWeight', 'bold');
            case 2
                title(ax{i}, 'Uncued', 'FontWeight', 'bold');
        end
    end
    add_time_tick(fig, ax{i}, time_this);
end
drawnow();
end % draw_scatter

% density plot
function ax = draw_density(obj, fig, pos, ax_sz, data_kde, fp, rw, color, ls, lw)
draw_sz = size(data_kde, 2);
ax = cell(1, draw_sz);
dist_w = obj.get_plot_dist(ax, pos, ax_sz);

x_lim = max(max(cellfun(@(x) max(x.f), data_kde)));
x_lim = (ceil(x_lim) + round(x_lim)) / 2;
if x_lim==0
    x_lim = 1;
end
for i = 1:draw_sz
    x_now = pos(1) + dist_w * (i-1);
    y_now = pos(2);

    ax{i} = axes(fig, "Units", "centimeters", "Position", [x_now y_now ax_sz], ...
        'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out');

    obj.plot_distr(ax{i}, data_kde(:,i), 'Color', color{:,i}, 'LineStyle', ls{:,i}, 'LineWidth', lw{:,i}, 'Reversed', 1);
    if ~isempty(fp)
        yline(ax{i}, fp, 'LineWidth', .5, 'LineStyle', '-');
        if ~isempty(rw)
            yline(ax{i}, fp+rw, 'LineWidth', .5, 'LineStyle', ':');
        end
    end

    set(ax{i}, 'XLim', [0 x_lim]);

    if i==1
        xlabel(ax{i}, 'Density (1/s)', 'FontWeight', 'bold');
    else
        set(ax{i}, 'YTickLabel', []);
    end
end
drawnow();
end % draw_density

% time ticks
function add_time_tick(fig, ax, x_val)
ax_tick = axes(fig, "Units", ax.Units, "Position", ax.Position, 'NextPlot', 'add', 'FontSize', 8, 'TickDir', 'out', ...
    'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', ax.XLim, 'YLim', [0 1]);
for i = 1:length(x_val)
    plot(ax_tick, x_val([i i]), [0 .05], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '-');
end
end % add_time_tick
