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

opts.plotsize = [8    4;
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
fig = figure(25); clf(25);
set(gcf, 'unit', 'centimeters', 'position', [2 2 21 12], 'paperpositionmode', 'auto', 'color', 'w');

uicontrol('Style', 'text', 'parent', fig, 'units', 'normalized', 'position', [0.3 0.95 0.4 0.04],...
    'string', OBJs{1}.Subject+" / "+Pair(1)+"_"+Pair(2), 'fontsize', 11, 'fontweight', 'bold', 'backgroundcolor', 'w');

%% Set axes and plot
%
ha_maze = axes();
set(ha_maze, 'units', 'centimeters', 'position', [1.5 1, 8 4], ...
    'nextplot', 'add', 'fontsize', 7, 'tickdir', 'out');
plot_diagram(ha_maze, opts);

ha_rt = axes;
set(ha_rt, 'units', 'centimeters', 'position', [1.5 6.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_reaction_time(ha_rt, OBJs, opts, 'left');

ha_rt1 = axes;
set(ha_rt1, 'units', 'centimeters', 'position', [6 6.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_reaction_time(ha_rt1, OBJs, opts, 'right');
set(ha_rt1, 'yticklabel', [], 'ylabel', []);

ha_rt.YLim(2) = 0.3; % max([ha2.YLim(2) ha21.YLim(2)]);
ha_rt1.YLim(2) = 0.3; % ha2.YLim(2);

%
ha_ct = axes;
set(ha_ct, 'units', 'centimeters', 'position', [12 1.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_choice_time(ha_ct, OBJs, opts, 'left');

ha_ct1 = axes;
set(ha_ct1, 'units', 'centimeters', 'position', [16.5 1.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_choice_time(ha_ct1, OBJs, opts, 'right');
set(ha_ct1, 'yticklabel', [], 'ylabel', []);

ha_ct.YLim(2) = 0.8; % max([ha_ct.YLim(2) ha_ct1.YLim(2)]);
ha_ct1.YLim(2) = 0.8; % ha_ct.YLim(2);

%
ha_mt = axes;
set(ha_mt, 'units', 'centimeters', 'position', [12 6.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_movement_time(ha_mt, OBJs, opts, 'left');
set(ha_mt, 'xlabel', [], 'xticklabel', []);

ha_mt1 = axes;
set(ha_mt1, 'units', 'centimeters', 'position', [16.5 6.5, 4 4], ...
    'nextplot', 'add', 'fontsize', 9, 'tickdir', 'out');
plot_movement_time(ha_mt1, OBJs, opts, 'right');
set(ha_mt1, 'yticklabel', [], 'ylabel', [], 'xlabel', [], 'xticklabel', []);

ha_mt.YLim(2) = 0.8; % max([ha_mt.YLim(2) ha_mt1.YLim(2)]);
ha_mt1.YLim(2) = 0.8; % ha_mt.YLim(2);

%%
savename = fullfile(TargetDir, "Comparison_Show_"+OBJs{1}.Subject+"_"+Pair(1)+"_"+Pair(2));
print(fig, '-dpdf', savename, '-bestfit')
print(fig, '-dpng', savename)
saveas(fig, savename, 'fig')

%%
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


%% Reaction time violin
function plot_reaction_time(ax, OBJs, opts, port)

sig_ax = axes('Units', ax.Units, 'Position', ax.Position, 'nextplot', 'add', ...
    'xlim', [.5 length(opts.MixedFP)+.5], 'YLim', [0 1], ...
    'Color', 'none', 'XColor', 'none', 'YColor', 'none');

yline(ax, .5, 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'LineStyle', '--')

num_violin = length(opts.Ports)*length(opts.MixedFP);
thisRT = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);
thisFP = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);
thisTask = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);

switch port
    case {'left'}
        p = 1;
    case {'right'}
        p = 2;
end

for fp_this = 1:length(opts.MixedFP)
    for p_this = 1:length(opts.Pair)
        thisFP(:, 2*(fp_this-1)+p_this) = fp_this;
        thisTask(:, 2*(fp_this-1)+p_this) = p_this;
        
        thisRT(1:length(OBJs{p_this}.RTSortedAll{fp_this, p}), 2*(fp_this-1)+p_this) = OBJs{p_this}.RTSortedAll{fp_this, p};
    end
end

switch port
    case {'left'}
        ax.Title.String = "Left";
        ax.Title.Color  = opts.color.PortL;
    case {'right'}
        ax.Title.String = "Right";
        ax.Title.Color  = opts.color.PortR;
end

means = mean(thisRT, 'omitnan');
sems  = std(thisRT, 'omitnan') ./ sqrt(sum(~isnan(thisRT), 'omitnan'));

pval = anovan(thisRT(:), {thisFP(:), thisTask(:)}, 'model', 'interaction', 'varnames', {'FP', 'Port'}, 'display', 'off');
%         if pval(3)<.05 % interaction
if pval(1)<.05 % control Port L, test FP
    [pval_1,~,stats_1] = anovan(thisRT(thisTask==1), thisFP(thisTask==1), 'varnames', {'FP'}, 'display', 'off');
    if pval_1<.05
        results_1 = multcompare(stats_1, 'Display', 'off');
        tbl_1 = array2table(results_1, ...
            "VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        sig_1 = find(tbl_1.("P-value")<.05);
        for s_this = 1:length(sig_1)
            sig_this = sig_1(s_this);
            if mean(means(2:2:end)) < mean(means(1:2:end)) % R < L
                y_this   = 1.01 - 0.08 * sig_this;
            else
                y_this   = -.05 + 0.08 * sig_this;
            end
            line(sig_ax, results_1(sig_this, 1:2)-0.08, [y_this y_this], 'Color', opts.color.(opts.Pair(1)), 'LineWidth', 1);
            if tbl_1.("P-value")(sig_this)>0.01
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "*", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            elseif tbl_1.("P-value")(sig_this)>0.001
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "**", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            else
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "***", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            end
        end
    end

    [pval_right,~,stats_right] = anovan(thisRT(thisTask==2), thisFP(thisTask==2), 'varnames', {'FP'}, 'display', 'off');
    if pval_right<.05
        results_right = multcompare(stats_right, 'Display', 'off');
        tbl_right = array2table(results_right, ...
            "VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        sig_right = find(tbl_right.("P-value")<.05);
        for s_this = 1:length(sig_right)
            sig_this = sig_right(s_this);
            if mean(means(2:2:end)) > mean(means(1:2:end)) % R > L
                y_this   = 1.01 - 0.08 * sig_this;
            else
                y_this   = -.05 + 0.08 * sig_this;
            end
            line(sig_ax, results_right(sig_this, 1:2)+0.08, [y_this y_this], 'Color', opts.color.(opts.Pair(2)), 'LineWidth', 1);
            if tbl_right.("P-value")(sig_this)>0.01
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "*", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            elseif tbl_right.("P-value")(sig_this)>0.001
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "**", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            else
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "***", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            end
        end
    end
end

medians = median(thisRT, 'omitnan');
iqr_pos = prctile(thisRT, 75)-medians;
iqr_neg = medians-prctile(thisRT, 25);

errorbar(ax, (1:3)+0.08, means(2:2:end), sems(2:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
    'Color', opts.color.(opts.Pair(2)), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.(opts.Pair(2)), 'MarkerSize', 5);
errorbar(ax, (1:3)-0.08, means(1:2:end), sems(1:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
    'Color', opts.color.(opts.Pair(1)), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.(opts.Pair(1)), 'MarkerSize', 5);

ax.Title.FontWeight = "Bold";
ax.XLabel.FontWeight = "Bold";
ax.YLabel.FontWeight = "Bold";
ax.XLabel.String = 'Foreperiod (s)';
ax.YLabel.String = 'Reaction time (s)';
set(ax, 'xlim', [.5 length(opts.MixedFP)+.5], 'ylim', [0 max(means+sems)+0.08], 'xtick', 1:length(opts.MixedFP), ...
    'xticklabel', opts.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
end

%% Choice time
function plot_choice_time(ax, OBJs, opts, port)

sig_ax = axes('Units', ax.Units, 'Position', ax.Position, 'nextplot', 'add', ...
    'xlim', [.5 length(opts.MixedFP)+.5], 'YLim', [0 1], ...
    'Color', 'none', 'XColor', 'none', 'YColor', 'none');

num_violin = length(opts.Ports)*length(opts.MixedFP);
thisCT = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);
thisFP = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);
thisTask = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);

switch port
    case {'left'}
        p = 1;
    case {'right'}
        p = 2;
end

for fp_this = 1:length(opts.MixedFP)
    for p_this = 1:length(opts.Pair)
        thisFP(:, 2*(fp_this-1)+p_this) = fp_this;
        thisTask(:, 2*(fp_this-1)+p_this) = p_this;
        
        this_ct = rmoutliers_custome(OBJs{p_this}.CTSortedAll{fp_this, p});

        thisCT(1:length(this_ct), 2*(fp_this-1)+p_this) = this_ct;
    end
end

switch port
    case {'left'}
        ax.Title.String = "Left";
        ax.Title.Color  = opts.color.PortL;
    case {'right'}
        ax.Title.String = "Right";
        ax.Title.Color  = opts.color.PortR;
end

means = mean(thisCT, 'omitnan');
sems  = std(thisCT, 'omitnan') ./ sqrt(sum(~isnan(thisCT), 'omitnan'));

pval = anovan(thisCT(:), {thisFP(:), thisTask(:)}, 'model', 'interaction', 'varnames', {'FP', 'Port'}, 'display', 'off');
%         if pval(3)<.05 % interaction
if pval(1)<.05 % control Port L, test FP
    [pval_1,~,stats_1] = anovan(thisCT(thisTask==1), thisFP(thisTask==1), 'varnames', {'FP'}, 'display', 'off');
    if pval_1<.05
        results_1 = multcompare(stats_1, 'Display', 'off');
        tbl_1 = array2table(results_1, ...
            "VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        sig_1 = find(tbl_1.("P-value")<.05);
        for s_this = 1:length(sig_1)
            sig_this = sig_1(s_this);
            if mean(means(2:2:end)) < mean(means(1:2:end)) % R < L
                y_this   = 1.01 - 0.08 * sig_this;
            else
                y_this   = -.05 + 0.08 * sig_this;
            end
            line(sig_ax, results_1(sig_this, 1:2)-0.08, [y_this y_this], 'Color', opts.color.(opts.Pair(1)), 'LineWidth', 1);
            if tbl_1.("P-value")(sig_this)>0.01
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "*", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            elseif tbl_1.("P-value")(sig_this)>0.001
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "**", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            else
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "***", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            end
        end
    end

    [pval_right,~,stats_right] = anovan(thisCT(thisTask==2), thisFP(thisTask==2), 'varnames', {'FP'}, 'display', 'off');
    if pval_right<.05
        results_right = multcompare(stats_right, 'Display', 'off');
        tbl_right = array2table(results_right, ...
            "VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        sig_right = find(tbl_right.("P-value")<.05);
        for s_this = 1:length(sig_right)
            sig_this = sig_right(s_this);
            if mean(means(2:2:end)) > mean(means(1:2:end)) % R > L
                y_this   = 1.01 - 0.08 * sig_this;
            else
                y_this   = -.05 + 0.08 * sig_this;
            end
            line(sig_ax, results_right(sig_this, 1:2)+0.08, [y_this y_this], 'Color', opts.color.(opts.Pair(2)), 'LineWidth', 1);
            if tbl_right.("P-value")(sig_this)>0.01
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "*", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            elseif tbl_right.("P-value")(sig_this)>0.001
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "**", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            else
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "***", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            end
        end
    end
end

medians = median(thisCT, 'omitnan');
iqr_pos = prctile(thisCT, 75)-medians;
iqr_neg = medians-prctile(thisCT, 25);

errorbar(ax, (1:3)+0.08, means(2:2:end), sems(2:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
    'Color', opts.color.(opts.Pair(2)), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.(opts.Pair(2)), 'MarkerSize', 5);
errorbar(ax, (1:3)-0.08, means(1:2:end), sems(1:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
    'Color', opts.color.(opts.Pair(1)), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.(opts.Pair(1)), 'MarkerSize', 5);

ax.Title.FontWeight = "Bold";
ax.XLabel.FontWeight = "Bold";
ax.YLabel.FontWeight = "Bold";
ax.XLabel.String = 'Foreperiod (s)';
ax.YLabel.String = 'Choice time (s)';
set(ax, 'xlim', [.5 length(opts.MixedFP)+.5], 'ylim', [0 max(means+sems)+0.08], 'xtick', 1:length(opts.MixedFP), ...
    'xticklabel', opts.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
end

%% Choice time
function plot_movement_time(ax, OBJs, opts, port)

sig_ax = axes('Units', ax.Units, 'Position', ax.Position, 'nextplot', 'add', ...
    'xlim', [.5 length(opts.MixedFP)+.5], 'YLim', [0 1], ...
    'Color', 'none', 'XColor', 'none', 'YColor', 'none');

num_violin = length(opts.Ports)*length(opts.MixedFP);
thisMT = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);
thisFP = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);
thisTask = nan(max([height(OBJs{1}.BehavTable), height(OBJs{2}.BehavTable)]), num_violin);

switch port
    case {'left'}
        p = 1;
    case {'right'}
        p = 2;
end

for fp_this = 1:length(opts.MixedFP)
    for p_this = 1:length(opts.Pair)
        thisFP(:, 2*(fp_this-1)+p_this) = fp_this;
        thisTask(:, 2*(fp_this-1)+p_this) = p_this;
        
        this_mt = rmoutliers_custome(OBJs{p_this}.MTSortedAll{fp_this, p});

        thisMT(1:length(this_mt), 2*(fp_this-1)+p_this) = this_mt;
    end
end
% thisMT = rmoutliers(thisMT, 1);

switch port
    case {'left'}
        ax.Title.String = "Left";
        ax.Title.Color  = opts.color.PortL;
    case {'right'}
        ax.Title.String = "Right";
        ax.Title.Color  = opts.color.PortR;
end

means = mean(thisMT, 'omitnan');
sems  = std(thisMT, 'omitnan') ./ sqrt(sum(~isnan(thisMT), 'omitnan'));

pval = anovan(thisMT(:), {thisFP(:), thisTask(:)}, 'model', 'interaction', 'varnames', {'FP', 'Port'}, 'display', 'off');
%         if pval(3)<.05 % interaction
if pval(1)<.05 % control Port L, test FP
    [pval_1,~,stats_1] = anovan(thisMT(thisTask==1), thisFP(thisTask==1), 'varnames', {'FP'}, 'display', 'off');
    if pval_1<.05
        results_1 = multcompare(stats_1, 'Display', 'off');
        tbl_1 = array2table(results_1, ...
            "VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        sig_1 = find(tbl_1.("P-value")<.05);
        for s_this = 1:length(sig_1)
            sig_this = sig_1(s_this);
            if mean(means(2:2:end)) < mean(means(1:2:end)) % R < L
                y_this   = 1.01 - 0.08 * sig_this;
            else
                y_this   = -.05 + 0.08 * sig_this;
            end
            line(sig_ax, results_1(sig_this, 1:2)-0.08, [y_this y_this], 'Color', opts.color.(opts.Pair(1)), 'LineWidth', 1);
            if tbl_1.("P-value")(sig_this)>0.01
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "*", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            elseif tbl_1.("P-value")(sig_this)>0.001
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "**", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            else
                text(sig_ax, mean(results_1(sig_this, 1:2)-0.08), y_this+.01, "***", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(1)));
            end
        end
    end

    [pval_right,~,stats_right] = anovan(thisMT(thisTask==2), thisFP(thisTask==2), 'varnames', {'FP'}, 'display', 'off');
    if pval_right<.05
        results_right = multcompare(stats_right, 'Display', 'off');
        tbl_right = array2table(results_right, ...
            "VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        sig_right = find(tbl_right.("P-value")<.05);
        for s_this = 1:length(sig_right)
            sig_this = sig_right(s_this);
            if mean(means(2:2:end)) > mean(means(1:2:end)) % R > L
                y_this   = 1.01 - 0.08 * sig_this;
            else
                y_this   = -.05 + 0.08 * sig_this;
            end
            line(sig_ax, results_right(sig_this, 1:2)+0.08, [y_this y_this], 'Color', opts.color.(opts.Pair(2)), 'LineWidth', 1);
            if tbl_right.("P-value")(sig_this)>0.01
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "*", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            elseif tbl_right.("P-value")(sig_this)>0.001
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "**", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            else
                text(sig_ax, mean(results_right(sig_this, 1:2)+0.08), y_this+.01, "***", ...
                    "FontSize", 14, 'HorizontalAlignment', 'center', 'Color', opts.color.(opts.Pair(2)));
            end
        end
    end
end


medians = median(thisMT, 'omitnan');
iqr_pos = prctile(thisMT, 75)-medians;
iqr_neg = medians-prctile(thisMT, 25);

errorbar(ax, (1:3)+0.08, means(2:2:end), sems(2:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
    'Color', opts.color.(opts.Pair(2)), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.(opts.Pair(2)), 'MarkerSize', 5);
errorbar(ax, (1:3)-0.08, means(1:2:end), sems(1:2:end), 'LineStyle', '-', 'CapSize', 2, 'LineWidth', 1.5, ...
    'Color', opts.color.(opts.Pair(1)), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color.(opts.Pair(1)), 'MarkerSize', 5);

ax.Title.FontWeight = "Bold";
ax.XLabel.FontWeight = "Bold";
ax.YLabel.FontWeight = "Bold";
ax.XLabel.String = 'Foreperiod (s)';
ax.YLabel.String = 'Movement time (s)';
set(ax, 'xlim', [.5 length(opts.MixedFP)+.5], 'ylim', [0 max(means+sems)+0.08], 'xtick', 1:length(opts.MixedFP), ...
    'xticklabel', opts.MixedFP, 'ticklength', [0.01 0.1], 'box', 'off');
end

