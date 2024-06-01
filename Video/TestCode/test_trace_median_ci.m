clear; close all;

%%
Treatment = "ChemoBiM2";
Protocol  = "ThreeFPHoldSRT";
Subject   = "Kennard";

Drives = string(char('A':'Z')');
for i = 1:26
    folder_i = Drives(i)+":\OneDrive";
    if isfolder(folder_i)
        WorkFolder = fullfile(folder_i, "YuLab\GPS_Project\Scripts");
    end
end

DataFolder = fullfile(WorkFolder, "Data", Treatment, Protocol, Subject);

FigFolder  = fullfile(WorkFolder, "Figures", Treatment, Protocol, Subject);
if ~isfolder(FigFolder)
    mkdir(FigFolder);
end

DataFile = dir(fullfile(DataFolder, "*Traj*.mat"));
if isempty(DataFile)
    error('No traj data files here');
elseif length(DataFile)>1
    error('More than one traj data file')
end

load(fullfile(DataFolder, DataFile.name));

%%
ci_alpha = 0.01; % 99ci

Trace    = struct();
Info     = obj.TrialInfo;
TargetFP = sort(obj.TargetFP, 'descend');

Labels   = ["Control", "Chemo"];
Features = ["AngleHead", "PosXHead", "PosYHead"];

%% calculate trace median with 99ci in fore-period
TimePoints = obj.TimeMatIn;
TimeBounds = num2cell(repmat(1000*TargetFP', 1, 2));
TimeBounds = cellfun(@(x) [0 x], TimeBounds, 'UniformOutput', false);

for lb = 1:length(Labels)
    lb_this   = Labels(lb);
    info_this = Info(Info.Label==lb_this, :);
    ind_this  = info_this.Stage==1 & info_this.Outcome=="Correct";

    sort_refs  = {info_this.FP(ind_this), info_this.PortCorrect(ind_this)};
    sort_codes = {TargetFP, obj.LeftRight};

    for f = 1:length(Features)
        feature_this = Features(f);
        trace_matrix = obj.TraceMatrix.(lb_this).In.(feature_this);
        trace_sorted = obj.sort_data(trace_matrix(ind_this, :), sort_refs, sort_codes);
        trace_median = cellfun(@(x, t_bound) obj.cal_trace_median(x, TimePoints, t_bound, ci_alpha), trace_sorted, TimeBounds, 'UniformOutput', false);

        Trace.(lb_this).FP.(feature_this) = trace_median;
    end
end

%% plot
plt = GPSPlot();

ax_sz = [2 1.8];
lw = 1.5;
ls = "-";
alpha = 0.35;

color_cell = repmat({{GPSColor.PortL, GPSColor.PortR}}, 3, 1);
lw_cell = repmat({{lw, lw}}, 3, 1);
ls_cell = repmat({{ls, ls}}, 3, 1);

%
fig = figure(12); clf(fig);
set(fig, 'Visible', 'on', 'Units', 'centimeters', 'Position', [5 5 18.5 8], 'Color', 'w', 'toolbar', 'none');

fig_title = sprintf("%s / %s / %s", Treatment, Protocol, Subject);
set_fig_title(fig, fig_title);

% Head angle
% Control
ax_ang_control   = plt.assign_ax_to_fig(fig, 3, 1, [1.5 1 2 6], ax_sz);
data_ang_control = plt.assign_data_to_ax(ax_ang_control, Trace.Control.FP.AngleHead);
ax_ang_control   = draw_trace_median(ax_ang_control, data_ang_control, color_cell, lw_cell, ls_cell, alpha);

for i = 1:2
    set(ax_ang_control{i}, 'YColor', 'none', 'XColor', 'none');
end
xlabel(ax_ang_control{3}, "Time from cent-in (ms)", 'FontWeight', 'bold');
ylabel(ax_ang_control{3}, "Head angle (°)", 'FontWeight', 'bold');

% Chemo
ax_ang_chemo   = plt.assign_ax_to_fig(fig, 3, 1, [4 1 2 6], ax_sz);
data_ang_chemo = plt.assign_data_to_ax(ax_ang_chemo, Trace.Chemo.FP.AngleHead);
ax_ang_chemo   = draw_trace_median(ax_ang_chemo, data_ang_chemo, color_cell, lw_cell, ls_cell, alpha);

y_lim_ang_1 = floor(min(cellfun(@(ax1, ax2) min([ax1.YLim(1) ax2.YLim(1)]), ax_ang_control, ax_ang_chemo)));
y_lim_ang_2 = ceil(max(cellfun(@(ax1, ax2) max([ax1.YLim(2) ax2.YLim(2)]), ax_ang_control, ax_ang_chemo)));

cellfun(@(x) set(x, 'XLim', [0 1500], 'YLim', [y_lim_ang_1 y_lim_ang_2], 'YColor', 'none', 'XColor', 'none'), ax_ang_chemo);
cellfun(@(x, fp) plt.add_shade(x, fp, 'alpha', .1, 'color', GPSColor.Treat), ax_ang_chemo, TimeBounds(:,1));
cellfun(@(x) set(x, 'XLim', [0 1500], 'YLim', [y_lim_ang_1 y_lim_ang_2]), ax_ang_control);
cellfun(@(x, fp) plt.add_shade(x, fp, 'alpha', .1, 'color', GPSColor.Control), ax_ang_control, TimeBounds(:,1));

% Head position x
% Control
ax_posx_control   = plt.assign_ax_to_fig(fig, 3, 1, [7.5 1 2 6], ax_sz);
data_posx_control = plt.assign_data_to_ax(ax_posx_control, Trace.Control.FP.PosXHead);
ax_posx_control   = draw_trace_median(ax_posx_control, data_posx_control, color_cell, lw_cell, ls_cell, alpha);

for i = 1:2
    set(ax_posx_control{i}, 'YColor', 'none', 'XColor', 'none');
end
ylabel(ax_posx_control{3}, "Head pos-x (pix)", 'FontWeight', 'bold');

% Chemo
ax_posx_chemo   = plt.assign_ax_to_fig(fig, 3, 1, [10 1 2 6], ax_sz);
data_posx_chemo = plt.assign_data_to_ax(ax_posx_chemo, Trace.Chemo.FP.PosXHead);
ax_posx_chemo   = draw_trace_median(ax_posx_chemo, data_posx_chemo, color_cell, lw_cell, ls_cell, alpha);

y_lim_posx_1 = floor(min(cellfun(@(ax1, ax2) min([ax1.YLim(1) ax2.YLim(1)]), ax_posx_control, ax_posx_chemo)));
y_lim_posx_2 = ceil(max(cellfun(@(ax1, ax2) max([ax1.YLim(2) ax2.YLim(2)]), ax_posx_control, ax_posx_chemo)));

cellfun(@(x) set(x, 'XLim', [0 1500], 'YLim', [y_lim_posx_1 y_lim_posx_2], 'YColor', 'none', 'XColor', 'none'), ax_posx_chemo);
cellfun(@(x, fp) plt.add_shade(x, fp, 'alpha', .1, 'color', GPSColor.Treat), ax_posx_chemo, TimeBounds(:,1));
cellfun(@(x) set(x, 'XLim', [0 1500], 'YLim', [y_lim_posx_1 y_lim_posx_2]), ax_posx_control);
cellfun(@(x, fp) plt.add_shade(x, fp, 'alpha', .1, 'color', GPSColor.Control), ax_posx_control, TimeBounds(:,1));

% Head position y
% Control
ax_posy_control   = plt.assign_ax_to_fig(fig, 3, 1, [13.5 1 2 6], ax_sz);
data_posy_control = plt.assign_data_to_ax(ax_posy_control, Trace.Control.FP.PosYHead);
ax_posy_control   = draw_trace_median(ax_posy_control, data_posy_control, color_cell, lw_cell, ls_cell, alpha);

for i = 1:2
    set(ax_posy_control{i}, 'YColor', 'none', 'XColor', 'none');
end
ylabel(ax_posy_control{3}, "Head pos-y (pix)", 'FontWeight', 'bold');

% Chemo
ax_posy_chemo   = plt.assign_ax_to_fig(fig, 3, 1, [16 1 2 6], ax_sz);
data_posy_chemo = plt.assign_data_to_ax(ax_posy_chemo, Trace.Chemo.FP.PosYHead);
ax_posy_chemo   = draw_trace_median(ax_posy_chemo, data_posy_chemo, color_cell, lw_cell, ls_cell, alpha);

y_lim_posy_1 = floor(min(cellfun(@(ax1, ax2) min([ax1.YLim(1) ax2.YLim(1)]), ax_posy_control, ax_posy_chemo)));
y_lim_posy_2 = ceil(max(cellfun(@(ax1, ax2) max([ax1.YLim(2) ax2.YLim(2)]), ax_posy_control, ax_posy_chemo)));

cellfun(@(x) set(x, 'XLim', [0 1500], 'YLim', [y_lim_posy_1 y_lim_posy_2], 'YColor', 'none', 'XColor', 'none'), ax_posy_chemo);
cellfun(@(x, fp) plt.add_shade(x, fp, 'alpha', .1, 'color', GPSColor.Treat), ax_posy_chemo, TimeBounds(:,1));
cellfun(@(x) set(x, 'XLim', [0 1500], 'YLim', [y_lim_posy_1 y_lim_posy_2]), ax_posy_control);
cellfun(@(x, fp) plt.add_shade(x, fp, 'alpha', .1, 'color', GPSColor.Control), ax_posy_control, TimeBounds(:,1));

fig_name = sprintf("Fig_Trace_ci_%s_%s_%s", Treatment, Protocol, Subject);
fig_path = fullfile(FigFolder, fig_name);
exportgraphics(fig, fig_path+".jpg", 'Resolution', 600);
exportgraphics(fig, fig_path+".pdf", 'ContentType', 'vector');

%%
function ax_cell = draw_trace_median(ax_cell, trace_cell, color_cell, lw_cell, ls_cell, alpha)
[n_row, n_col] = size(ax_cell);
for i = 1:n_row
    for j = 1:n_col
        ax_now = ax_cell{i,j};
        trace_now = trace_cell{i,j};
        color_now = color_cell{i,j};
        lw_now = lw_cell{i,j};
        ls_now = ls_cell{i,j};

        ax_cell{i,j} = plot_trace_median(ax_now, trace_now, color_now, lw_now, ls_now, alpha);
    end
end
end

function ax = plot_trace_median(ax, trace, color, lw, ls, alpha)

    for i = 1:length(trace)
        fill(ax, [trace{i}.time flip(trace{i}.time)], [trace{i}.ci(1,:) flip(trace{i}.ci(2,:))], 'r', 'FaceColor', color{i}, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    end
    for i = 1:length(trace)
        plot(ax, trace{i}.time, trace{i}.trace, 'Color', color{i}, 'LineWidth', lw{i}, 'LineStyle', ls{i});
    end

end

