function heatmapAngleHeadKb(ax, obj, opts)

arguments
    ax  (1, 1) handle
    obj (1, 1) 
    opts.color (1, 1) GPSColor = GPSColor()
    opts.PortCorrect (1, 1) {mustBeMember(opts.PortCorrect, ["L", "R"])} = "L"
    opts.Performance (1, 1) {mustBeMember(opts.Performance, ["Correct", "Wrong", "Late", "Premature", "Valid", "All"])} = "All"
    opts.Cued (1, 1) {mustBeMember(opts.Cued, ["Cue", "Uncue", "All"])} = "Uncue"
    opts.SortBy (1, 1) {mustBeMember(opts.SortBy, ["HD", "MT", "FP"])} = "HD"
    opts.AlignTo (1, 1) {mustBeMember(opts.AlignTo, ["In", "Out"])} = "In"
    opts.Label (1, 1) {mustBeMember(opts.Label, ["All", "None", "Saline", "Control", "Chemo"])} = "All"
end

ind = obj.Ind;

%% extract by performance
switch opts.Performance
    case {'All'}

        these_trials = find(ind.("Target"+opts.PortCorrect));

        ax.Title.String = "Target Port"+opts.PortCorrect;
        ax.Title.Color  = opts.color.("Port"+opts.PortCorrect);

    case {'Premature', 'Late'}

        these_trials = find(ind.("Target"+opts.PortCorrect) & ind.(opts.Performance));
        ax.Title.String = opts.Performance;

    case {'Valid'}

        these_trials = find(ind.("Target"+opts.PortCorrect) & (ind.Wrong | ind.Correct));

        yline(ax, sum(ind.("Target"+opts.PortCorrect) & ind.Correct) + 0.5, ':', 'LineWidth', 1);
        ax.Title.String = "Target Port"+opts.PortCorrect;
        ax.Title.Color  = opts.color.("Port"+opts.PortCorrect);

    otherwise

        these_trials = find(ind.("Target"+opts.PortCorrect) & ind.(opts.Performance));

        ax.Title.String = "Target Port"+opts.PortCorrect;
        ax.Title.Color  = opts.color.("Port"+opts.PortCorrect);

end

ax.Title.FontWeight = 'bold';

%% extract by experiment label
if ~strcmp(opts.Label, "All")
    these_trials = these_trials(obj.TrialInfo.Label(these_trials)==opts.Label);
end

if length(these_trials) > 100
    these_trials = these_trials(randperm(length(these_trials), 100));
end

%% sort
switch opts.Performance

    case {'All', 'Correct', 'Wrong'}
        [~, sort_id] = sort(obj.(opts.SortBy)(these_trials));
        these_trials = these_trials(sort_id);

    case {'Valid'}
        [~, sort_id] = sort(obj.(opts.SortBy)(these_trials));
        these_trials = these_trials(sort_id);

        [~, sort_id] = sort(obj.Performance(these_trials));
        these_trials = these_trials(sort_id);

    case {'Premature', 'Late'}
        [~, sort_id] = sort(obj.(opts.SortBy)(these_trials));
        these_trials = these_trials(sort_id);

        [~, sort_id] = sort(cellfun(@(x) x(end)>0, obj.AngleHead(these_trials)));
        these_trials = these_trials(sort_id);
end

%%
M = obj.("AngleHeadMat"+opts.AlignTo)(these_trials, :);

if isempty(M)
    set(ax, 'ylim', [0.5 1.5], 'ytick', [], 'xticklabel', [], 'ycolor', 'none');
    ax.Position(4) = 0.05*ax.Position(4);
    switch opts.AlignTo
        case "In"
            point_in = find(obj.TimePointsIn==0);
            set(ax, 'xtick', point_in + (-500:500:2500), 'xticklabel', ["-500", "0", "500", "1000", "1500", "2000", "2500"]);

        case "Out"
            point_out = find(obj.TimePointsOut==0);
            set(ax, 'xtick', point_out + (-2500:500:600), 'xticklabel', ["-2500", "-2000", "-1500", "-1000", "-500", "0", "500"]);
    end
    ax.XLabel.String     = "Time from poke "+opts.AlignTo+" (ms)";
    ax.XLabel.FontWeight = "Bold";

    set(ax, 'ydir', 'reverse', 'xlim', [0 length(obj.("TimePoints"+opts.AlignTo))]);

    return
end

switch opts.AlignTo
    case "In"
        point_in = find(obj.TimePointsIn==0);
        xline(ax, point_in, '-', 'LineWidth', 1.2);

        imagesc(ax, M, "AlphaData", ~isnan(M));

        for i = 1:length(these_trials)
            point_out    = point_in + 1000*obj.HD(these_trials(i));
            point_choose = point_in + 1000*(obj.HD(these_trials(i)) + obj.MT(these_trials(i)));
            point_fp     = point_in + 1000*obj.FP(these_trials(i));
            line(ax, [point_out point_out], [i-.5 i+.5], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', 'k');
            line(ax, [point_choose point_choose], [i-.5 i+.5], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', opts.color.(obj.Performance(these_trials(i))));
            line(ax, [point_fp point_fp], [i-.4 i+.4], 'LineStyle', '-', 'LineWidth', 1.2, 'Color', [.3 .3 .3]);
        end

        set(ax, 'xtick', point_in + (-500:500:2500), 'xticklabel', ["-500", "0", "500", "1000", "1500", "2000", "2500"]);

    case "Out"
        point_out = find(obj.TimePointsOut==0);
        xline(ax, point_out, '-', 'LineWidth', 1.2);

        imagesc(ax, M, "AlphaData", ~isnan(M));

        for i = 1:length(these_trials)
            point_in     = point_out - 1000*obj.HD(these_trials(i));
            point_choose = point_out + 1000*obj.MT(these_trials(i));
            point_fp     = point_in + 1000*obj.FP(these_trials(i));
            line(ax, [point_in point_in], [i-.5 i+.5], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', 'k');
            line(ax, [point_choose point_choose], [i-.5 i+.5], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', opts.color.(obj.Performance(these_trials(i))));
            line(ax, [point_fp point_fp], [i-.4 i+.4], 'LineStyle', '-', 'LineWidth', 1.2, 'Color', [.3 .3 .3]);
        end

        set(ax, 'xtick', point_out + (-2500:500:600), 'xticklabel', ["-2500", "-2000", "-1500", "-1000", "-500", "0", "500"]);
end

clim(ax, [-90 90]);

ax.Position(4) = min([ax.Position(4) size(M, 1)*0.08]);

ax.XLabel.String     = "Time from poke "+opts.AlignTo+" (ms)";
ax.XLabel.FontWeight = "Bold";
ax.YLabel.String     = "Trial #";
ax.YLabel.FontWeight = "Bold";

if length(these_trials)<10
    set(ax, 'ytick', 1:5:10);
    ax.YLabel.String = "";
elseif length(these_trials)<40
    set(ax, 'ytick', 0:5:40);
elseif length(these_trials)<80
    set(ax, 'ytick', 0:10:80);
elseif length(these_trials)<160
    set(ax, 'ytick', 0:20:160);
end
set(ax, 'ydir', 'reverse', 'ylim', [0.5 length(these_trials)+0.5], ...
    'xlim', [0 length(obj.("TimePoints"+opts.AlignTo))]);
end