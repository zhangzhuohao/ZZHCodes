clear;

TargetDir = uigetdir('E:\YuLab\Work\GPS\Data\');

[~, Pair] = fileparts(TargetDir);
Pair = string(split(Pair, '_'));

%%
ProtocolDir = get_folders(TargetDir, "FolderType", 'Protocol');
[~, Tasks] = arrayfun(@(x) fileparts(x), ProtocolDir);

check = any(contains(Tasks, Pair(1)), 2);
if ~check(1)
    ProtocolDir = flipud(ProtocolDir);
end

%%
OBJs = cell(length(ProtocolDir), 1);

for d = 1:length(ProtocolDir)
    load(get_mat_files(ProtocolDir(d), "FileType", 'ProgressClass'));
    OBJs{d} = obj;
end

%%
opts.color = struct(GPSColor()); % Color class for GPS
opts.color.(Pair(1)) = [40 120 181] / 255;
opts.color.(Pair(2)) = [224 72 75] / 255;

opts.mk = ["o", "x"]; % Scatter marker for correct and others
opts.ls = [":", "-.", "-"]; % Line style for [ShortFP, MedFP, LongFP]
opts.lw = [1, 1.5, 2]; % Line width for [ShortFP, MedFP, LongFP]

opts.plotsize = [8   4;
                8.5  3.5;
                8    2;
                8.5  5;
                5    5;
                4    4];

opts.sep_col = 1.2;
opts.sep_row = 1.5;

opts.MixedFP = OBJs{1}.MixedFP;
opts.Ports = OBJs{1}.Ports;
opts.Pair = Pair;

%%
cp.session_1 = OBJs{1}.Sessions;
cp.session_2 = OBJs{2}.Sessions;

cp.date_1    = char(cp.session_1);
cp.date_1    = string(cp.date_1(:, 5:8));
cp.date_2    = char(cp.session_2);
cp.date_2    = string(cp.date_2(:, 5:8));

if length(cp.session_2) ~= length(cp.session_1)
    fprintf("The numbers of 1 and 2 sessions dont match.\n");
    return
end

cp.behav_1   = OBJs{1}.BehavTable;
cp.behav_2   = OBJs{2}.BehavTable;
cp.ind_1     = OBJs{1}.Ind;
cp.ind_2     = OBJs{2}.Ind;

cp.trial_1   = cp.behav_1.TrialStartTime + cp.behav_1.CentInTime;
cp.trial_2   = cp.behav_2.TrialStartTime + cp.behav_2.CentInTime;

cp.num_session_pairs = length(cp.session_1);
cp.session_sep       = zeros(1, cp.num_session_pairs-1);
for s = 1:cp.num_session_pairs-1
    i_1   = find(cp.behav_1.SessionDate==cp.session_1(s), 1, 'last');
    i_2   = find(cp.behav_2.SessionDate==cp.session_2(s), 1, 'last');
    end_1 = cp.trial_1(i_1);
    end_2 = cp.trial_2(i_2);

    cp.session_sep(s) = max([end_1, end_2]) + 5;

    if s==1
        cp.trial_1(i_1+1:end) = cp.trial_1(i_1+1:end) + cp.session_sep(s);
        cp.trial_2(i_2+1:end)     = cp.trial_2(i_2+1:end) + cp.session_sep(s);
    else
        cp.trial_1(i_1+1:end) = cp.trial_1(i_1+1:end) + cp.session_sep(s) - cp.session_sep(s-1);
        cp.trial_2(i_2+1:end)     = cp.trial_2(i_2+1:end) + cp.session_sep(s) - cp.session_sep(s-1);
    end
end

%%
fig = figure(24); clf(24);
set(gcf, 'unit', 'centimeters', 'position', [2 1.5 35.5 23.4], 'paperpositionmode', 'auto', 'color', 'w');

uicontrol('Style', 'text', 'parent', 24, 'units', 'normalized', 'position', [0.3 0.95 0.4 0.04],...
    'string', OBJs{1}.Subject+" / "+Pair(1)+"_"+Pair(2), 'fontsize', 11, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set axes and plot
%% Maze diagram
ha1 = axes;
set(ha1, 'units', 'centimeters', 'position', [2.8+2*(opts.sep_row+opts.plotsize(1, 1))+5.5 1.5+3*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(1, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_diagram(ha1, opts);

%% Performance
% Hold duration scatter
% Control
ha21 = axes;
set(ha21, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize(1, 1)) 1.5+3*(opts.sep_col+opts.plotsize(1, 2))+2.7, opts.plotsize(3, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_scatter(ha21, cp, 1, opts);
set(ha21, 'xtick', [], 'xlabel', []);

% Chemo
ha22 = axes;
set(ha22, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize(1, 1)) 1.5+3*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(3, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_scatter(ha22, cp, 2, opts);

% Performance track
ha3 = axes;
set(ha3, 'units', 'centimeters', 'position', [1.5+1*(opts.sep_row+opts.plotsize(1, 1)) 1.5+3*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(4, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_performance_track(ha3, OBJs, cp, opts);

% Performance compare
ha4 = axes;
set(ha4, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1)) 1.5+3*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(5, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_performance_compare(ha4, OBJs, opts);

ha41 = axes;
set(ha41, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1))+3 1.5+3*(opts.sep_col+opts.plotsize(1, 2))+.5, 0.4*opts.plotsize(5, :) ], 'nextplot', 'add', 'fontsize', 7);
plot_performance_compare(ha41, OBJs, opts);
set(ha41, 'xlim', [0 20], 'ylim', [0 20], 'xlabel', [], 'ylabel', [], 'xtick', 0:5:20, 'ytick', 0:5:20, 'title', []);

% legend
ha42 = axes;
set(ha42, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1))+5.1 1.5+3*(opts.sep_col+opts.plotsize(1, 2))+2, [1.5 3] ], 'nextplot', 'add', 'fontsize', 7);
line(ha42, [0 .6], [0 0], 'Color', opts.color.Correct, 'LineWidth', 1);
text(ha42, .7, 0, 'Correct', 'FontSize', 8, 'Color', opts.color.Correct);
line(ha42, [0 .6], [1 1], 'Color', opts.color.Premature, 'LineWidth', 1);
text(ha42, .7, 1, 'Premature', 'FontSize', 8, 'Color', opts.color.Premature);
line(ha42, [0 .6], [2 2], 'Color', opts.color.Wrong, 'LineWidth', 1);
text(ha42, .7, 2, 'Wrong', 'FontSize', 8, 'Color', opts.color.Wrong);
line(ha42, [0 .6], [3 3], 'Color', opts.color.Late, 'LineWidth', 1);
text(ha42, .7, 3, 'Late', 'FontSize', 8, 'Color', opts.color.Late);

line(ha42, [0 .6], [5 5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
text(ha42, .7, 5, Pair(1), 'FontSize', 8, 'Color', 'k');
line(ha42, [0 .6], [6 6], 'Color', 'k', 'LineWidth', 1, 'LineStyle', ':')
text(ha42, .7, 6, Pair(2), 'FontSize', 8, 'Color', 'k');

set(ha42, 'xlim', [0 2], 'ylim', [0 7], 'xcolor', 'none', 'ycolor', 'none', 'ydir', 'reverse', 'color', 'none');

%% Hold duration
% Hold duration violin plot
ha5 = axes;
set(ha5, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize(1, 1)) 1.5+2*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(1, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_violin(ha5, OBJs, opts);

% Hold duration statistic
% Median
ha61 = axes;
set(ha61, 'units', 'centimeters', 'position', [1.5+1*(opts.sep_row+opts.plotsize(1, 1)) 1.5+2*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(6, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_median_compare(ha61, OBJs, opts);

% IQR
ha62 = axes;
set(ha62, 'units', 'centimeters', 'position', [1.8+1*(opts.sep_row+opts.plotsize(1, 1))+4.5 1.5+2*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(6, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_iqr_compare(ha62, OBJs, opts);
set(ha62, 'ylabel', []);

% Hold duration PDF
% Right
ha71 = axes;
set(ha71, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize(1, 1)) 1.5+1*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(2, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_pdf(ha71, OBJs, "R", opts);

% Left
ha72 = axes;
set(ha72, 'units', 'centimeters', 'position', [1.5+0*(opts.sep_row+opts.plotsize(1, 1)) 1.5+0*(opts.sep_col+opts.plotsize(1, 2))+0.5, opts.plotsize(2, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_pdf(ha72, OBJs, "L", opts);

ha71.YLim(2) = max([ha71.YLim(2) ha72.YLim(2)]);
ha72.YLim(2) = ha71.YLim(2);

% Hold duration CDF
% Right
ha81 = axes;
set(ha81, 'units', 'centimeters', 'position', [1.5+1*(opts.sep_row+opts.plotsize(1, 1)) 1.5+1*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(2, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_cdf(ha81, OBJs, "R", opts);

% Left
ha82 = axes;
set(ha82, 'units', 'centimeters', 'position', [1.5+1*(opts.sep_row+opts.plotsize(1, 1)) 1.5+0*(opts.sep_col+opts.plotsize(1, 2))+0.5, opts.plotsize(2, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_hold_duration_cdf(ha82, OBJs, "L", opts);

%% Reaction time
% Reaction time violin
ha9 = axes;
set(ha9, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1)) 1.5+2*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(1, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_violin(ha9, OBJs, opts);

% reaction time compare
ha10 = axes;
set(ha10, 'units', 'centimeters', 'position', [2+3*(opts.sep_row+opts.plotsize(1, 1)) 1.5+2*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(6, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_reaction_time_median_compare(ha10, OBJs, opts);

%% Movement time
% movement time violin
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1)) 1.5+1*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(1, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time_violin(ha11, OBJs, opts);

% movement time median
ha12 = axes;
set(ha12, 'units', 'centimeters', 'position', [2+3*(opts.sep_row+opts.plotsize(1, 1)) 1.5+1*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(6, :) ], 'nextplot', 'add', 'fontsize', 8);
plot_movement_time_median_compare(ha12, OBJs, opts);

%% Shuttle time
% log shuttle time violin
ha13 = axes;
set(ha13, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1)) 1.5+0*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(1, :) .* [1/3, 1] ], 'nextplot', 'add', 'fontsize', 8);
plot_shuttle_time_violin(ha13, OBJs, cp, opts);

%% Interruption
% Interruption histogram
ha14 = axes;
set(ha14, 'units', 'centimeters', 'position', [2+2*(opts.sep_row+opts.plotsize(1, 1))+4.5 2+0*(opts.sep_col+opts.plotsize(1, 2)), opts.plotsize(2, :) .* [3/5 1] ], 'nextplot', 'add', 'fontsize', 8);
plot_interruption(ha14, OBJs, cp, opts);

%%
savename = fullfile(TargetDir, "Comparison_"+OBJs{1}.Subject+"_"+Pair(1)+"_"+Pair(2));
print(fig, '-dpdf', savename, '-bestfit')
print(fig, '-dpng', savename)
saveas(fig, savename, 'fig')

%% ha1. Make a diagram of the setup
function plot_diagram(ax, opts)

x = [1 3 5 8 9   11  11  9   8 5 3 1 1];
y = [0 0 2 2 1.2 1.2 4.8 4.8 4 4 6 6 0];

patch(ax, 'XData', x, 'YData', y, 'FaceColor', 'none', 'EdgeColor', 'k', 'linewidth', 2);

viscircles(ax, [9.6, 3], 0.3, 'color', 'k', 'LineWidth', 1);
text(ax, 8.6, 3, 'Init', 'FontWeight','bold', 'HorizontalAlignment', 'center');
viscircles(ax, [3.5, 3], 0.3,  'color', 'k', 'LineWidth', 1);
text(ax, 3.9, 3, 'Cent', 'FontWeight', 'bold');

viscircles(ax, [3, 1.5], 0.3, 'color', opts.color.PortL);
text(ax, 1.2, 1.5, 'Left', 'FontWeight','bold', 'HorizontalAlignment', 'left', 'Color', opts.color.PortL);
viscircles(ax, [3, 4.5], 0.3,  'color', opts.color.PortR);
text(ax, 1.2, 4.5, 'Right', 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Color', opts.color.PortR);

set(ax, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ylim', [0 6], 'xlim', [0 12]);
end

%% ha2. Plot every trial's hold duration during training
function plot_hold_duration_scatter(ax, cp, lb, opts)

for fp_this = 1:length(opts.MixedFP)
    yline(ax, opts.MixedFP(fp_this), 'Color', [.8 .8 .8], 'LineWidth', opts.lw(fp_this), 'LineStyle', ':', 'Alpha', 0.4);
end
xline(ax, cp.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-');

switch lb
    case 1
        hd  = cp.behav_1.HoldDuration;
        fp  = cp.behav_1.FP;
        ind = cp.ind_1;
        tk  = cp.trial_1;
    case 2
        hd  = cp.behav_2.HoldDuration;
        fp  = cp.behav_2.FP;
        ind = cp.ind_2;
        tk  = cp.trial_2;
end

% port 1 wrong (defined by their action, which is different from target)
scatter(ax, tk(ind.wrongL), hd(ind.wrongL), ...
    22*fp(ind.wrongL), opts.color.PortL, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

% port 2 wrong
scatter(ax, tk(ind.wrongR), hd(ind.wrongR), ...
    22*fp(ind.wrongR), opts.color.PortR, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

% premature
scatter(ax, tk(ind.premature), hd(ind.premature), ...
    22*fp(ind.premature), opts.color.Premature, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

% late
scatter(ax, tk(ind.late), hd(ind.late), ...
    22*fp(ind.late), opts.color.Late, opts.mk(2), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

% port 1 correct (defined by their action, which is also the target for correct response)
scatter(ax, tk(ind.correctL), hd(ind.correctL), ...
    18*fp(ind.correctL), opts.color.PortL, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

% port 2 correct
scatter(ax, tk(ind.correctR), hd(ind.correctR), ...
    18*fp(ind.correctR), opts.color.PortR, opts.mk(1), 'MarkerEdgeAlpha', 0.6, 'linewidth', 1);

switch lb
    case 1
        t = "\it{" + opts.Pair(1) + " ";
        for s_this = 1:cp.num_session_pairs
            if s_this < cp.num_session_pairs
                t = t + cp.date_1(s_this) + "/";
            else
                t = t + cp.date_1(s_this) + "}";
            end
        end
        ax.Title.String = t;
        ax.Title.Color  = opts.color.(opts.Pair(1));
    case 2
        t = "\it{" + opts.Pair(2) + " ";
        for s_this = 1:cp.num_session_pairs
            if s_this < cp.num_session_pairs
                t = t + cp.date_2(s_this) + "/";
            else
                t = t + cp.date_2(s_this) + "}";
            end
        end
        ax.Title.String = t;
        ax.Title.Color  = opts.color.(opts.Pair(2));
end
ax.Title.FontWeight  = 'bold';
ax.Title.FontSize    = 9;
ax.Title.Interpreter = 'tex';

ax.XLabel.String = 'Time in training (s)';
ax.YLabel.String = 'Hold duration (s)';
set(ax, 'xlim', [0 max([cp.trial_1(end) cp.trial_2(end)])+5], 'ylim', [0 2.5], 'ticklength', [0.01 0.1]);
end

%% ha3. Plot performance track of each session
function plot_performance_track(ax, OBJs, cp, opts)

xline(ax, cp.session_sep, 'Color', [.7 .7 .7], 'LineWidth', 1, 'LineStyle', '-')

s_sep = [0 cp.session_sep];

for s_this = 1:cp.num_session_pairs
    ind_1 = find(OBJs{1}.PerformanceTrack.Sessions==cp.session_1(s_this));
    plot(ax, OBJs{1}.PerformanceTrack.WinPos(ind_1)+s_sep(s_this), OBJs{1}.PerformanceTrack.CorrectRatio(ind_1),   'linestyle', '-', 'color', opts.color.Correct, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Correct,   'markeredgecolor', 'w');
    plot(ax, OBJs{1}.PerformanceTrack.WinPos(ind_1)+s_sep(s_this), OBJs{1}.PerformanceTrack.WrongRatio(ind_1),     'linestyle', '-', 'color', opts.color.Wrong, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Wrong,     'markeredgecolor', 'w');
    plot(ax, OBJs{1}.PerformanceTrack.WinPos(ind_1)+s_sep(s_this), OBJs{1}.PerformanceTrack.PrematureRatio(ind_1), 'linestyle', '-', 'color', opts.color.Premature, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
    plot(ax, OBJs{1}.PerformanceTrack.WinPos(ind_1)+s_sep(s_this), OBJs{1}.PerformanceTrack.LateRatio(ind_1),      'linestyle', '-', 'color', opts.color.Late, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Late,      'markeredgecolor', 'w');

    ind_2 = find(OBJs{2}.PerformanceTrack.Sessions==cp.session_2(s_this));
    plot(ax, OBJs{2}.PerformanceTrack.WinPos(ind_2)+s_sep(s_this), OBJs{2}.PerformanceTrack.CorrectRatio(ind_2),   'linestyle', ':', 'color', opts.color.Correct, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Correct,   'markeredgecolor', 'w');
    plot(ax, OBJs{2}.PerformanceTrack.WinPos(ind_2)+s_sep(s_this), OBJs{2}.PerformanceTrack.WrongRatio(ind_2),     'linestyle', ':', 'color', opts.color.Wrong, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Wrong,     'markeredgecolor', 'w');
    plot(ax, OBJs{2}.PerformanceTrack.WinPos(ind_2)+s_sep(s_this), OBJs{2}.PerformanceTrack.PrematureRatio(ind_2), 'linestyle', ':', 'color', opts.color.Premature, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Premature, 'markeredgecolor', 'w');
    plot(ax, OBJs{2}.PerformanceTrack.WinPos(ind_2)+s_sep(s_this), OBJs{2}.PerformanceTrack.LateRatio(ind_2),      'linestyle', ':', 'color', opts.color.Late, ...
        'markersize', 5, 'linewidth', 1.5, 'markerfacecolor', opts.color.Late,      'markeredgecolor', 'w');
end

ax.XLabel.String = 'Time in training (s)';
ax.YLabel.String = 'Performance (%)';
set(ax, 'XLim', [0 max([cp.trial_1(end) cp.trial_2(end)])+5], 'YLim', [0 100]);
end

%% ha4. Plot performance comparation
function plot_performance_compare(ax, OBJs, opts)

line(ax, [0 100], [0 100], 'LineStyle', ':', 'LineWidth', 1, 'Color', 'k');

perf_1 = OBJs{1}.PerformanceAll;
perf_2 = OBJs{2}.PerformanceAll;

r_short = OBJs{1}.PerformanceAll.Foreperiod==0.5 & OBJs{1}.PerformanceAll.TargetPort=="R";
l_short = OBJs{1}.PerformanceAll.Foreperiod==0.5 & OBJs{1}.PerformanceAll.TargetPort=="L";
r_med   = OBJs{1}.PerformanceAll.Foreperiod==1   & OBJs{1}.PerformanceAll.TargetPort=="R";
l_med   = OBJs{1}.PerformanceAll.Foreperiod==1   & OBJs{1}.PerformanceAll.TargetPort=="L";
r_long  = OBJs{1}.PerformanceAll.Foreperiod==1.5 & OBJs{1}.PerformanceAll.TargetPort=="R";
l_long  = OBJs{1}.PerformanceAll.Foreperiod==1.5 & OBJs{1}.PerformanceAll.TargetPort=="L";

r_all   = OBJs{1}.PerformanceAll.Foreperiod==0   & OBJs{1}.PerformanceAll.TargetPort=="R";
l_all   = OBJs{1}.PerformanceAll.Foreperiod==0   & OBJs{1}.PerformanceAll.TargetPort=="L";

scatter(ax, perf_1.PrematureRatio(r_short | r_med | r_long), perf_2.PrematureRatio(r_short | r_med | r_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Premature, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
scatter(ax, perf_1.PrematureRatio(l_short | l_med | l_long), perf_2.PrematureRatio(l_short | l_med | l_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Premature, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
%         plot(ax, [perf_1.PrematureRatio(r_all) perf_2.PrematureRatio(r_all)], [perf_1.PrematureRatio(l_all) perf_2.PrematureRatio(l_all)], ...
%             'Color', [opts.color.Premature 0.6], 'LineWidth', 1.5);

scatter(ax, perf_1.WrongRatio(r_short | r_med | r_long), perf_2.WrongRatio(r_short | r_med | r_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Wrong, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
scatter(ax, perf_1.WrongRatio(l_short | l_med | l_long), perf_2.WrongRatio(l_short | l_med | l_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Wrong, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
%         plot(ax, [perf_1.WrongRatio(r_all) perf_2.WrongRatio(r_all)], [perf_1.WrongRatio(l_all) perf_2.WrongRatio(l_all)], ...
%             'Color', [opts.color.Wrong 0.6], 'LineWidth', 1.5);

scatter(ax, perf_1.LateRatio(r_short | r_med | r_long), perf_2.LateRatio(r_short | r_med | r_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Late, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
scatter(ax, perf_1.LateRatio(l_short | l_med | l_long), perf_2.LateRatio(l_short | l_med | l_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Late, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
%         plot(ax, [perf_1.LateRatio(r_all) perf_2.LateRatio(r_all)], [perf_1.LateRatio(l_all) perf_2.LateRatio(l_all)], ...
%             'Color', [opts.color.Late 0.6], 'LineWidth', 1.5);

scatter(ax, perf_1.CorrectRatio(r_short | r_med | r_long), perf_2.CorrectRatio(r_short | r_med | r_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Correct, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
scatter(ax, perf_1.CorrectRatio(l_short | l_med | l_long), perf_2.CorrectRatio(l_short | l_med | l_long), ...
    [.5 1 1.5]*24, 'Marker', '^', 'MarkerFaceColor', opts.color.Correct, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.5, ...
    'LineWidth', 1);
%         plot(ax, [perf_1.CorrectRatio(r_all) perf_2.CorrectRatio(r_all)], [perf_1.CorrectRatio(l_all) perf_2.CorrectRatio(l_all)], ...
%             'Color', [opts.color.Correct 0.6], 'LineWidth', 1.5);

scatter(ax, perf_1.PrematureRatio(r_all), perf_2.PrematureRatio(r_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Premature, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, perf_1.PrematureRatio(l_all), perf_2.PrematureRatio(l_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Premature, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

scatter(ax, perf_1.WrongRatio(r_all), perf_2.WrongRatio(r_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Wrong, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, perf_1.WrongRatio(l_all), perf_2.WrongRatio(l_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Wrong, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

scatter(ax, perf_1.LateRatio(r_all), perf_2.LateRatio(r_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Late, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, perf_1.LateRatio(l_all), perf_2.LateRatio(l_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Late, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

scatter(ax, perf_1.CorrectRatio(r_all), perf_2.CorrectRatio(r_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Correct, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, perf_1.CorrectRatio(l_all), perf_2.CorrectRatio(l_all), 48, ...
    'Marker', 'o', 'MarkerFaceColor', opts.color.Correct, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

text(ax, 50, 100, "Performance (%)", 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

ax.XLabel.String     = opts.Pair(1);
ax.XLabel.Color      = opts.color.(opts.Pair(1));
ax.XLabel.FontWeight = "Bold";

ax.YLabel.String     = opts.Pair(2);
ax.YLabel.Color      = opts.color.(opts.Pair(2));
ax.YLabel.FontWeight = "Bold";

set(ax, 'xlim', [0 100], 'ylim', [0 100]);
end

%% Hold duration violin
function plot_hold_duration_violin(ax, OBJs, opts)

dcz_ind = 2 * (1:length(opts.MixedFP))';
num_dcz = length(opts.MixedFP);
patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
    'FaceColor', opts.color.(opts.Pair(2)), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

for fp_this = 1:length(opts.MixedFP)
    line(ax, -.5 + 2*[fp_this-0.45 fp_this+0.45], [opts.MixedFP(fp_this) opts.MixedFP(fp_this)], 'Color', 'k', 'LineWidth', .5, 'LineStyle', '-');
end

num_violin = length(opts.Ports)*2*length(opts.MixedFP);
thisHD = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);

for fp_this = 1:length(opts.MixedFP)
    for p_this = 1:length(opts.Ports)
        thisHD(1:length(OBJs{1}.HDSortedAll{fp_this, p_this}), 4*(fp_this-1)+p_this)   = OBJs{1}.HDSortedAll{fp_this, p_this};
        thisHD(1:length(OBJs{2}.HDSortedAll{fp_this, p_this}), 4*(fp_this-1)+p_this+2) = OBJs{2}.HDSortedAll{fp_this, p_this};
    end
end

thisHD(all(isnan(thisHD), 2), :) = [];

violinplot({thisHD(:, 2:2:end), thisHD(:, 1:2:end)}, 1:num_violin, ...
    'ViolinColor', {repmat(opts.color.PortR, num_violin/2, 1), repmat(opts.color.PortL, num_violin/2, 1)}, ...
    'ScatterSize', 8, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false, 'Bandwidth', 0.1);

scatter(ax, ceil((1:num_violin)/2) + repmat([-.02 .02], 1, num_violin/2), median(thisHD, 'omitnan'), 32, ...
    repmat([opts.color.PortL; opts.color.PortR], num_violin/2, 1), ...
    'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3], 'MarkerFaceAlpha', 0.6);

ax.XLabel.String = 'Foreperiod (s)';
ax.YLabel.String = 'Hold duration (s)';
set(ax, 'xlim', [.5 2*length(opts.MixedFP)+.5], 'ylim', [0 OBJs{1}.Bins.HoldDuration(end)], 'xtick', 1.5:2:2*length(opts.MixedFP), ...
    'xticklabel', opts.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
end

%% Hold duration comparation
% median
function plot_hold_duration_median_compare(ax, OBJs, opts)

ind_l = find(ismember(OBJs{1}.HDStatAll.thisFP, opts.MixedFP) & OBJs{1}.HDStatAll.Port=="L");
ind_r = find(ismember(OBJs{1}.HDStatAll.thisFP, opts.MixedFP) & OBJs{1}.HDStatAll.Port=="R");

%         for i = 1:length(ind_l)
%             plot(ax, [OBJs{1}.HDStatAll.Median(ind_l(i)) OBJs.HDStatChemo.Median(ind_l(i))], [OBJs{1}.HDStatAll.Median(ind_r(i)) OBJs.HDStatChemo.Median(ind_r(i))], ...
%                 'Color', [0 0 0 0.6], 'LineWidth', 1.5);
%         end

scatter(ax, OBJs{1}.HDStatAll.Median(ind_l), OBJs{2}.HDStatAll.Median(ind_l), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortL, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, OBJs{1}.HDStatAll.Median(ind_r), OBJs{2}.HDStatAll.Median(ind_r), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortR, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

ax.YLabel.String     = opts.Pair(2);
ax.YLabel.Color      = opts.color.(opts.Pair(2));
ax.YLabel.FontWeight = "Bold";

ax.XLabel.String     = opts.Pair(1);
ax.XLabel.Color      = opts.color.(opts.Pair(1));
ax.XLabel.FontWeight = "Bold";

set(ax, 'xlimmode', 'auto', 'ylimmode', 'auto');
ax.XLim = [min([ax.XLim(1) ax.YLim(1)]) max([ax.XLim(2) ax.YLim(2)])] .* [.95 1.05];
ax.YLim = ax.XLim;
ax.XTick = ax.YTick;

text(ax, mean(ax.XLim), ax.YLim(2), "HD median (s)", 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

line(ax, [0 100], [0 100], 'Color', 'k', 'LineStyle', ':', 'LineWidth', .5);
end

% IQR
function plot_hold_duration_iqr_compare(ax, OBJs, opts)

ind_l = find(ismember(OBJs{1}.HDStatAll.thisFP, opts.MixedFP) & OBJs{1}.HDStatAll.Port=="L");
ind_r = find(ismember(OBJs{1}.HDStatAll.thisFP, opts.MixedFP) & OBJs{1}.HDStatAll.Port=="R");

%         for i = 1:length(ind_l)
%             plot(ax, [OBJs{1}.HDStatAll.IQR(ind_l(i)) OBJs.HDStatChemo.IQR(ind_l(i))], [OBJs{1}.HDStatAll.IQR(ind_r(i)) OBJs.HDStatChemo.IQR(ind_r(i))], ...
%                 'Color', [0 0 0 0.6], 'LineWidth', 1.5);
%         end

ax.YLabel.String     = opts.Pair(2);
ax.YLabel.Color      = opts.color.(opts.Pair(2));
ax.YLabel.FontWeight = "Bold";

ax.XLabel.String     = opts.Pair(1);
ax.XLabel.Color      = opts.color.(opts.Pair(1));
ax.XLabel.FontWeight = "Bold";

scatter(ax, OBJs{1}.HDStatAll.IQR(ind_l), OBJs{2}.HDStatAll.IQR(ind_l), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortL, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, OBJs{1}.HDStatAll.IQR(ind_r), OBJs{2}.HDStatAll.IQR(ind_r), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortR, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

set(ax, 'xlimmode', 'auto', 'ylimmode', 'auto');
ax.XLim = [min([ax.XLim(1) ax.YLim(1)]) max([ax.XLim(2) ax.YLim(2)])] .* [.95 1.05];
ax.YLim = ax.XLim;
ax.XTick = ax.YTick;

text(ax, mean(ax.XLim), ax.YLim(2), "HD IQR (s)", 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

line(ax, [0 100], [0 100], 'Color', 'k', 'LineStyle', ':', 'LineWidth', .5);
end

%% Hold duration density
% Prob. density
function plot_hold_duration_pdf(ax, OBJs, port, opts)

p_this = opts.Ports==upper(string(port));

for fp_this = 1:length(opts.MixedFP)

    xline(ax, opts.MixedFP(fp_this), 'color', [.7 .7 .7], 'linewidth', opts.lw(fp_this), 'LineStyle', '-');

    hd_1 = OBJs{1}.HDSortedAll{fp_this, p_this};
    hd_2 = OBJs{2}.HDSortedAll{fp_this, p_this};

    hd_1_pdf = ksdensity(hd_1, OBJs{1}.Bins.HoldDuration, 'Function', 'pdf');
    hd_2_pdf = ksdensity(hd_2, OBJs{1}.Bins.HoldDuration, 'Function', 'pdf');

    plot(ax, OBJs{1}.Bins.HoldDuration, hd_1_pdf, ...
        'color', opts.color.(opts.Pair(1)), 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
    plot(ax, OBJs{1}.Bins.HoldDuration, hd_2_pdf, ...
        'color', opts.color.(opts.Pair(2)), 'linewidth', opts.lw(fp_this), 'LineStyle', '-');

end

switch upper(string(port))
    case {"L"}
        ax.XLabel.String = "Hold duration Left (s)";
        ax.XLabel.Color  = opts.color.PortL;
    case {"R"}
        ax.XLabel.String = "Hold duration Right (s)";
        ax.XLabel.Color  = opts.color.PortR;
end
ax.XLabel.FontWeight = "Bold";

ax.YLabel.String = "Prob. density (1/s)";
set(ax, 'xlim', [0 OBJs{1}.Bins.HoldDuration(end)], 'ylimmode', 'auto');
end

% Cum. density
function plot_hold_duration_cdf(ax, OBJs, port, opts)

p_this = opts.Ports==upper(string(port));

for fp_this = 1:length(opts.MixedFP)

    xline(ax, opts.MixedFP(fp_this), 'color', [.7 .7 .7], 'linewidth', opts.lw(fp_this), 'LineStyle', '-');

    hd_1 = OBJs{1}.HDSortedAll{fp_this, p_this};
    hd_2 = OBJs{2}.HDSortedAll{fp_this, p_this};

    hd_1_pdf = ksdensity(hd_1, OBJs{1}.Bins.HoldDuration, 'Function', 'cdf');
    hd_2_pdf = ksdensity(hd_2, OBJs{1}.Bins.HoldDuration, 'Function', 'cdf');

    plot(ax, OBJs{1}.Bins.HoldDuration, hd_1_pdf, ...
        'color', opts.color.(opts.Pair(1)), 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
    plot(ax, OBJs{1}.Bins.HoldDuration, hd_2_pdf, ...
        'color', opts.color.(opts.Pair(2))  , 'linewidth', opts.lw(fp_this), 'LineStyle', '-');

end

switch upper(string(port))
    case {"L"}
        ax.XLabel.String = "Hold duration Left (s)";
        ax.XLabel.Color  = opts.color.PortL;
    case {"R"}
        ax.XLabel.String = "Hold duration Right (s)";
        ax.XLabel.Color  = opts.color.PortR;
end
ax.XLabel.FontWeight = "Bold";

ax.YLabel.String = "Cum. density";
set(ax, 'xlim', [0 OBJs{1}.Bins.HoldDuration(end)], 'ylim', [0 1]);
end

%% Reaction time violin
function plot_reaction_time_violin(ax, OBJs, opts)

dcz_ind = 2 * (1:length(opts.MixedFP))';
num_dcz = length(opts.MixedFP);
patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
    'FaceColor', opts.color.(opts.Pair(2)), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--')

num_violin = length(opts.Ports)*2*length(opts.MixedFP);
thisRT = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);

for fp_this = 1:length(opts.MixedFP)
    for p_this = 1:length(opts.Ports)
        thisRT(1:length(OBJs{1}.RTSortedAll{fp_this, p_this}), 4*(fp_this-1)+p_this)   = OBJs{1}.RTSortedAll{fp_this, p_this};
        thisRT(1:length(OBJs{2}.RTSortedAll{fp_this, p_this})  , 4*(fp_this-1)+p_this+2) = OBJs{2}.RTSortedAll{fp_this, p_this};
    end
end

thisRT(all(isnan(thisRT), 2), :) = [];

violinplot({thisRT(:, 2:2:end), thisRT(:, 1:2:end)}, 1:num_violin, ...
    'ViolinColor', {repmat(opts.color.PortR, num_violin/2, 1), repmat(opts.color.PortL, num_violin/2, 1)}, ...
    'ScatterSize', 8, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false, 'Bandwidth', 0.1);

scatter(ax, ceil((1:num_violin)/2) + repmat([-.02 .02], 1, num_violin/2), median(thisRT, 'omitnan'), 32, ...
    repmat([opts.color.PortL; opts.color.PortR], num_violin/2, 1), ...
    'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3], 'MarkerFaceAlpha', 0.6);

ax.XLabel.String = 'Foreperiod (s)';
ax.YLabel.String = 'Reaction time (s)';
set(ax, 'xlim', [.5 2*length(opts.MixedFP)+.5], 'ylim', [0 OBJs{1}.Bins.RT(end)], 'xtick', 1.5:2:2*length(opts.MixedFP), ...
    'xticklabel', opts.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
end

%% Reaction time comparation
% median
function plot_reaction_time_median_compare(ax, OBJs, opts)

ind_l = find(ismember(OBJs{1}.RTStatAll.thisFP, opts.MixedFP) & OBJs{1}.RTStatAll.Port=="L");
ind_r = find(ismember(OBJs{1}.RTStatAll.thisFP, opts.MixedFP) & OBJs{1}.RTStatAll.Port=="R");

%         for i = 1:length(ind_l)
%             plot(ax, [OBJs{1}.RTStatAll.Median(ind_l(i)) OBJs{2}.RTStatAll.Median(ind_l(i))], [OBJs{1}.RTStatAll.Median(ind_r(i)) OBJs{2}.RTStatAll.Median(ind_r(i))], ...
%                 'Color', [0 0 0 0.6], 'LineWidth', 1.5);
%         end

scatter(ax, OBJs{1}.RTStatAll.Median(ind_l), OBJs{2}.RTStatAll.Median(ind_l), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortL, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, OBJs{1}.RTStatAll.Median(ind_r), OBJs{2}.RTStatAll.Median(ind_r), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortR, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

ax.YLabel.String     = opts.Pair(2);
ax.YLabel.Color      = opts.color.(opts.Pair(2));
ax.YLabel.FontWeight = "Bold";

ax.XLabel.String     = opts.Pair(1);
ax.XLabel.Color      = opts.color.(opts.Pair(1));
ax.XLabel.FontWeight = "Bold";

set(ax, 'xlimmode', 'auto', 'ylimmode', 'auto');
ax.XLim = [min([ax.XLim(1) ax.YLim(1)]) max([ax.XLim(2) ax.YLim(2)])] .* [.95 1.05];
ax.YLim = ax.XLim;
ax.YTick = ax.XTick;

text(ax, mean(ax.XLim), ax.YLim(2), "RT median (s)", 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

line(ax, [0 100], [0 100], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
end

%% Movement time violin
function plot_movement_time_violin(ax, OBJs, opts)

dcz_ind = 2 * (1:length(opts.MixedFP))';
num_dcz = length(opts.MixedFP);
patch(ax, 'XData', [dcz_ind-.5 dcz_ind-.5 dcz_ind+.5 dcz_ind+.5]', 'YData', repmat([0; 100; 100; 0], 1, num_dcz), ...
    'FaceColor', opts.color.(opts.Pair(2)), 'FaceAlpha', 0.2, 'EdgeColor', 'none');


num_violin = length(opts.Ports)*2*length(opts.MixedFP);
thisMT = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);

for fp_this = 1:length(opts.MixedFP)
    for p_this = 1:length(opts.Ports)
        thisMT(1:length(OBJs{1}.MTSortedAll{fp_this, p_this}), 4*(fp_this-1)+p_this)   = OBJs{1}.MTSortedAll{fp_this, p_this};
        thisMT(1:length(OBJs{2}.MTSortedAll{fp_this, p_this})  , 4*(fp_this-1)+p_this+2) = OBJs{2}.MTSortedAll{fp_this, p_this};
    end
end

thisMT(all(isnan(thisMT), 2), :) = [];

violinplot({thisMT(:, 2:2:end), thisMT(:, 1:2:end)}, 1:num_violin, ...
    'ViolinColor', {repmat(opts.color.PortR, num_violin/2, 1), repmat(opts.color.PortL, num_violin/2, 1)}, ...
    'ScatterSize', 8, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false, 'Bandwidth', 0.1);

scatter(ax, ceil((1:num_violin)/2) + repmat([-.02 .02], 1, num_violin/2), median(thisMT, 'omitnan'), 32, ...
    repmat([opts.color.PortL; opts.color.PortR], num_violin/2, 1), ...
    'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3], 'MarkerFaceAlpha', 0.6);

ax.XLabel.String = 'Foreperiod (s)';
ax.YLabel.String = 'Movement time (s)';
set(ax, 'xlim', [.5 2*length(opts.MixedFP)+.5], 'ylim', [0 OBJs{1}.Bins.MovementTime(end)], 'xtick', 1.5:2:2*length(opts.MixedFP), ...
    'xticklabel', opts.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
end

%% Movement time comparation
function plot_movement_time_median_compare(ax, OBJs, opts)

ind_l = find(ismember(OBJs{1}.MTStatAll.thisFP, opts.MixedFP) & OBJs{1}.MTStatAll.Port=="L");
ind_r = find(ismember(OBJs{1}.MTStatAll.thisFP, opts.MixedFP) & OBJs{1}.MTStatAll.Port=="R");

%         for i = 1:length(ind_l)
%             plot(ax, [OBJs{1}.MTStatAll.Median(ind_l(i)) OBJs{2}.MTStatAll.Median(ind_l(i))], [OBJs{1}.MTStatAll.Median(ind_r(i)) OBJs{2}.MTStatAll.Median(ind_r(i))], ...
%                 'Color', [0 0 0 0.6], 'LineWidth', 1.5);
%         end

scatter(ax, OBJs{1}.MTStatAll.Median(ind_l), OBJs{2}.MTStatAll.Median(ind_l), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortL, 'MarkerEdgeColor', opts.color.PortL, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);
scatter(ax, OBJs{1}.MTStatAll.Median(ind_r), OBJs{2}.MTStatAll.Median(ind_r), ...
    [.5 1 1.5]*32, 'Marker', 'o', 'MarkerFaceColor', opts.color.PortR, 'MarkerEdgeColor', opts.color.PortR, ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5);

ax.YLabel.String     = opts.Pair(2);
ax.YLabel.Color      = opts.color.(opts.Pair(2));
ax.YLabel.FontWeight = "Bold";

ax.XLabel.String     = opts.Pair(1);
ax.XLabel.Color      = opts.color.(opts.Pair(1));
ax.XLabel.FontWeight = "Bold";

set(ax, 'xlimmode', 'auto', 'ylimmode', 'auto');
ax.XLim = [min([ax.XLim(1) ax.YLim(1)]) max([ax.XLim(2) ax.YLim(2)])] .* [.95 1.05];
ax.YLim = ax.XLim;
ax.YTick = ax.XTick;

text(ax, mean(ax.XLim), ax.YLim(2), "MT median (s)", 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

line(ax, [0 100], [0 100], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
end

%% log Shuttle time violin
function plot_shuttle_time_violin(ax, OBJs, cp, opts)

patch(ax, 'XData', [1.5 1.5 2.5 2.5]', 'YData', [-1; 100; 100; -1], ...
    'FaceColor', opts.color.(opts.Pair(2)), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

thisSTLog = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), 2);

STLog_1 = log(cp.behav_1.ShuttleTime);
[~, ~, indrmv] = rmoutliers_custome(STLog_1);
STLog_1(indrmv) = nan;

STLog_2 = log(cp.behav_2.ShuttleTime);
[~, ~, indrmv] = rmoutliers_custome(STLog_2);
STLog_2(indrmv) = nan;

thisSTLog(1:length(STLog_1), 1) = STLog_1;
thisSTLog(1:length(STLog_2), 2)  = STLog_2;

violinplot(thisSTLog, [opts.Pair(1), opts.Pair(2)], ...
    'ViolinColor', [opts.color.(opts.Pair(1)); opts.color.(opts.Pair(2))], ...
    'ScatterSize', 2, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false);

scatter(ax, 1:2, median(thisSTLog, 'omitnan'), 40, ...
    ones(2, 3), ...
    'filled', 'LineWidth', 1, 'MarkerEdgeColor', [.3 .3 .3]);

ax.XLabel.String = 'Treatment';
ax.YLabel.String = 'Shuttle time (s)';
set(ax, 'xlim', [.5 2.5], 'ylim', [-1 OBJs{1}.Bins.ShuttleTimeLog(end)], 'xtick', 1:2, ...
    'xticklabel', [opts.Pair(1), opts.Pair(2)], ...
    'ytick', -1:3, 'yticklabel', ["10^-1", "10^0", "10^1", "10^2", "10^3"], ...
    'ticklength', [0.01 0.1], 'box', 'off');
end

%% Interruption
function plot_interruption(ax, OBJs, cp, opts)

for fp_this = 1:length(opts.MixedFP)
    xline(ax, opts.MixedFP(fp_this), 'color', [.7 .7 .7], 'linewidth', opts.lw(fp_this), 'LineStyle', '-');
end

nums_1 = zeros(1, length(OBJs{1}.Bins.Interruption)-1);
nums_2   = zeros(1, length(OBJs{1}.Bins.Interruption)-1);
for i = 1:length(nums_1)
    nums_1(i) = sum(cp.behav_1.HoldDuration(cp.behav_1.Stage==1) >= OBJs{1}.Bins.Interruption(i));
    nums_2(i) = sum(cp.behav_2.HoldDuration(cp.behav_2.Stage==1) >= OBJs{1}.Bins.Interruption(i));
end

h_1 = histcounts(OBJs{1}.InterruptionNone.On, "BinEdges", OBJs{1}.Bins.Interruption) ./ nums_1;
bar(ax, OBJs{1}.Bins.Interruption(1:end-1)+.5*OBJs{1}.Bins.widthInter, h_1, 1, 'FaceColor', opts.color.(opts.Pair(1)) , 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'EdgeAlpha', 0.6);
h_2   = histcounts(OBJs{2}.InterruptionNone.On  , "BinEdges", OBJs{1}.Bins.Interruption) ./ nums_2;
bar(ax, OBJs{1}.Bins.Interruption(1:end-1)+.5*OBJs{1}.Bins.widthInter, h_2  , 1, 'FaceColor', opts.color.(opts.Pair(2))   , 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'EdgeAlpha', 0.6);

y_1 = smoothdata(h_1, 'gaussian', 0.18*(length(h_1)));
y_2   = smoothdata(h_2  , 'gaussian', 0.18*(length(h_2)));

plot(ax, OBJs{1}.Bins.Interruption(1:end-1)+.5*OBJs{1}.Bins.widthInter, y_1, 'Color', opts.color.(opts.Pair(1)), 'LineStyle', '-', 'LineWidth', 1.5);
plot(ax, OBJs{1}.Bins.Interruption(1:end-1)+.5*OBJs{1}.Bins.widthInter, y_2  , 'Color', opts.color.(opts.Pair(2))  , 'LineStyle', '-', 'LineWidth', 1.5);

ax.YLabel.String = "Interruptions / trial";
ax.XLabel.String = "Hold duration (s)";
set(ax, 'xlim', [0 1.5]);
end





