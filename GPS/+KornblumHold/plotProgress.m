function fig_progress = plotProgress(obj)
%%
beh = obj.BehavTable;
beh_sorted = obj.sort_table(beh(beh.Stage==1, :), obj.SortVars, obj.SortCodes);
beh_sorted_st = obj.sort_table(beh(beh.Stage==1, :), [], []);
beh_sorted_mt = obj.sort_table(beh(beh.Stage==1 & beh.Outcome=="Correct", :), obj.SortVars, obj.SortCodes);
C = GPSColor();

ax_sz_1 = [6 2];    
ax_sz_2 = [2.9 2];
ax_sz_3 = [2 2];

if contains(obj.Protocol, "Self")
    ls = "-";
    lw = 1;
    cued = 0;
    sort_id = 2;
    rw = .8;
    ax_sz_s = ax_sz_1;
    n_s = 1;
    cue_uncue_code = 0;
    cue_uncue_label = "Uncued";
else
    ls = ["-", ":"];
    lw = [1 1.25];
    cued = [1 0];
    sort_id = [1 2];
    rw = [.5 .8];
    ax_sz_s = ax_sz_2;
    n_s = 2;
    cue_uncue_code = [0 1];
    cue_uncue_label = ["Uncued" "Cued"];
end

%%
fig_progress = figure(12); clf(fig_progress);
set(fig_progress, 'Visible', 'on', 'Units', 'centimeters', 'Position', [5 5 29 13.2], 'Color', 'w', 'toolbar', 'none');

fig_title = sprintf("%s / %s / %s - %s", obj.Subject, obj.Protocol, obj.Sessions(1), obj.Sessions(end));
set_fig_title(fig_progress, fig_title);

% performance track
ax_perf_track = obj.assign_ax_to_fig(fig_progress, 1, 2, [1.2 10 12.5 2], ax_sz_1);
perf_track = obj.splice_data(obj.PerformanceTrack.Session);
data_perf_track = obj.assign_data_to_ax(ax_perf_track, perf_track(sort_id, :));
draw_perf_track(obj, ax_perf_track, data_perf_track, repmat({ls}, 1, 2), repmat({lw}, 1, 2));
title(ax_perf_track{1,1}, 'Left' , 'FontWeight', 'bold', 'Color', C.PortL);
title(ax_perf_track{1,2}, 'Right', 'FontWeight', 'bold', 'Color', C.PortR);

% performance progress
ax_perf = obj.assign_ax_to_fig(fig_progress, 1, 2, [1.2 7 12.5 2], ax_sz_1);
text(ax_perf{1}, .5, .5, "performance progress", 'HorizontalAlignment', 'center');

% pdf heatmap of hold duration
ax_heat = obj.assign_ax_to_fig(fig_progress, 1, 2, [1.2 4 12.5 2], ax_sz_1);
text(ax_heat{1}, .5, .5, "heat map", 'HorizontalAlignment', 'center');

% pdf hold duration in early and late training phase
ax_e_l = obj.assign_ax_to_fig(fig_progress, 1, 2, [1.2 1 12.5 2], ax_sz_1);
text(ax_e_l{1}, .5, .5, "early late", 'HorizontalAlignment', 'center');

% Shuttle time
% scatter
ax_st_s = obj.assign_ax_to_fig(fig_progress, 1, 1, [15 10 5 2], ax_sz_1);
data_st = obj.assign_data_to_ax(ax_st_s, beh_sorted_st);

text(ax_st_s{1}, .5, .5, "LogST scatter", 'HorizontalAlignment', 'center');

% violin
ax_st_v = obj.assign_ax_to_fig(fig_progress, 1, 1, [22.2 10 5 2], ax_sz_1);
text(ax_st_v{1}, .5, .5, "LogST violin", 'HorizontalAlignment', 'center');

% Movement time
% scatter
ax_mt_s = obj.assign_ax_to_fig(fig_progress, 1, n_s, [15 7 5 2], ax_sz_s);
data_mt = obj.assign_data_to_ax(ax_mt_s, beh_sorted_mt(sort_id, :)');

text(ax_mt_s{1}, .5, .5, "MT scatter", 'HorizontalAlignment', 'center');

% violin
ax_mt_v = obj.assign_ax_to_fig(fig_progress, 1, n_s, [22.2 7 5 2], ax_sz_s);
text(ax_mt_v{1}, .5, .5, "MT violin", 'HorizontalAlignment', 'center');

% Hold duration
% scatter
ax_hd_s = obj.assign_ax_to_fig(fig_progress, 1, n_s, [15 4 5 2], ax_sz_s);
data_hd = obj.assign_data_to_ax(ax_hd_s, beh_sorted(sort_id, :)');
c_hd  = repmat({[C.PortL; C.PortR]}, 1, n_s);
ms_hd = repmat({[8; 8]}, 1, n_s);
for i = 1:n_s
    set(ax_hd_s{i}, "YLim", [-1 1]+obj.TargetFP);
end
draw_scatter(obj, ax_hd_s, data_hd, "TrialCentInTimeProgress", "HD", c_hd, ms_hd, obj.TargetFP, rw);

% violin
ax_hd_v = obj.assign_ax_to_fig(fig_progress, 1, n_s, [22.2 4 5 2], ax_sz_s);
text(ax_hd_v{1}, .5, .5, "HD violin", 'HorizontalAlignment', 'center');

% Hold duration statistics
% median
ax_median = obj.assign_ax_to_fig(fig_progress, 1, 1, [15 1 5 2], ax_sz_1);
text(ax_median{1}, .5, .5, "HD median", 'HorizontalAlignment', 'center');

% iqr
ax_iqr = obj.assign_ax_to_fig(fig_progress, 1, 1, [22.2 1 5 2], ax_sz_1);
text(ax_iqr{1}, .5, .5, "HD IQR", 'HorizontalAlignment', 'center');

%%
end

%%
function ax_cell = draw_perf_track(obj, ax_cell, perf_track_cell, ls, lw)
session_info  = obj.DurationSession;
session_sep   = [0; cumsum(session_info)];
session_chemo = find(obj.Label=="Chemo");
sep_lesion    = find(obj.Label=="Lesion", 1);
sz_draw = size(ax_cell);
for i = 1:sz_draw(1)
    for j = 1:sz_draw(2)

        ax = ax_cell{i,j};
        mark_sessions(ax, session_sep, session_chemo, sep_lesion);
        for k = 1:length(perf_track_cell{i,j})
            perf_track = perf_track_cell{i,j}{k};
            perf_track = obj.add_progress_info(perf_track, "Pos", session_info);
            plot_perf_track(obj, ax, perf_track, ls{i,j}(k), lw{i,j}(k))
        end
        set(ax, 'XLim', [0 obj.BehavTable.TrialCentInTimeProgress(end)], 'YLim', [0 100]);

        if i==sz_draw(1) && j==1
            xlabel(ax, 'Time (s)', 'FontWeight', 'bold');
            ylabel(ax, 'Perfomance %', 'FontWeight', 'bold');
        end
        if i~=sz_draw(i)
            set(ax, 'XTickLabel', []);
        end
        if j~=1
            set(ax, 'YTickLabel', []);
        end
    end % for j = 1:sz_draw(2)
end % for i = 1:sz_draw(1)
drawnow();
end % draw_perf_track

function draw_scatter(obj, ax_cell, beh_cell, var_x, var_y, c, mk_sz, fp, rw)
session_info  = obj.DurationSession;
session_sep   = [0; cumsum(session_info)];
session_chemo = find(obj.Label=="Chemo");
sep_lesion    = find(obj.Label=="Lesion", 1);
sz_draw = size(ax_cell);
for i = 1:sz_draw(1)
    for j = 1:sz_draw(2)

        ax = ax_cell{i,j};
        mark_sessions(ax, session_sep, session_chemo, sep_lesion);
        if ~isempty(fp)
            yline(ax, fp, 'LineWidth', .5, 'LineStyle', '-', 'Alpha', .5);
            yline(ax, fp+rw(j), 'LineWidth', .5, 'LineStyle', ':', 'Alpha', .5);
        end
        plot_scatter(ax, beh_cell{i,j}, var_x, var_y, c{i,j}, mk_sz{i,j});
        set(ax, 'XLim', [0 obj.BehavTable.TrialCentInTimeProgress(end)]);

        if i==sz_draw(1) && j==1
            xlabel(ax, 'Time (s)', 'FontWeight', 'bold');
        end
        if i~=sz_draw(i)
            set(ax, 'XTickLabel', []);
        end
        if j~=1
            set(ax, 'YTickLabel', []);
        end
    end % for j = 1:sz_draw(2)
end % for i = 1:sz_draw(1)
drawnow();
end % draw_scatter

%% Plot functions
function plot_perf_track(obj, ax, perf_track, ls, lw)
sessions = unique(perf_track.Session, "stable");
for i = 1:length(sessions)
    ind_this = perf_track.Session == sessions(i);
    track_this = perf_track(ind_this, :);
    for j = 1:length(obj.PerformanceType)
        perf_this = obj.PerformanceType(j);
        plot(ax, track_this.PosProgress, 100*track_this.(perf_this), ...
            'LineStyle', ls, 'LineWidth', lw, 'Color', GPSColor.(perf_this));
    end
end
end % plot_perf_track

function plot_scatter(ax, beh, var_x, var_y, c, mk_sz)
n_beh = cellfun(@height, beh);
beh_this = vertcat(beh{:});

x_this = beh_this.(var_x);
y_this = beh_this.(var_y);

c_this  = zeros(sum(n_beh), 3);
ms_this = zeros(sum(n_beh), 1);
count_1 = 1;
for i = 1:length(n_beh)
    count_2 = count_1+n_beh(i)-1;
    c_this(count_1:count_2, :) = repmat(c(i, :), n_beh(i), 1);
    ms_this(count_1:count_2)   = repmat(mk_sz(i), n_beh(i), 1);
    count_1 = count_2 + 1;
end
[x, id] = sort(x_this);
y = y_this(id);
c = c_this(id,:);
ms = ms_this(id);

scatter(ax, x, y, ms, c, 'filled', 'MarkerFaceAlpha', .6);

end

%% Marking sessions
function mark_sessions(ax, session_sep, session_chemo, sep_lesion)
xline(ax, session_sep(2:end-1), 'LineWidth', .5, 'LineStyle', '-', 'Color', [.7 .7 .7], 'Alpha', .7);

if ~isempty(session_chemo)
    for k = 1:length(session_chemo)
        s_this = session_sep([0 1] + session_chemo(k));
        fill(ax, [s_this(1) s_this(2) s_this(2) s_this(1)], [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)], 'r', ...
            'FaceColor', GPSColor.Treat, 'FaceAlpha', .2, 'EdgeColor', 'none');
    end
end
if ~isempty(sep_lesion)
    xline(ax, session_sep(sep_lesion), 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0 0 0], 'Alpha', .7);
end
end