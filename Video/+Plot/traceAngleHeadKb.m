function traceAngleHeadKb(ax, obj, AlignTo, opts)

arguments
    ax  (1, 1) handle
    obj (1, 1)
    AlignTo (1, 1) {mustBeMember(AlignTo, ["In", "Out"])} = "In"
    opts.color (1, 1) GPSColor = GPSColor()
    opts.Condition (1, 1) {mustBeMember(opts.Condition, ["Cue", "Uncue"])} = "Uncue"
    opts.PortCorrect (1, 1) {mustBeMember(opts.PortCorrect, ["L", "R", "Both"])} = "Both"
    opts.Performance (1, 1) {mustBeMember(opts.Performance, ["Correct", "Wrong", "Valid", "All"])} = "All"
    opts.Label (1, 1) {mustBeMember(opts.Label, ["All", "None", "Saline", "Control", "Chemo"])} = "All"
end

switch opts.PortCorrect
    case {'Both'}
        PortCorrect = ["L", "R"];
    otherwise
        PortCorrect = opts.PortCorrect;
end

switch opts.Performance
    case {'All'}
        Perfs = ["Premature", "Wrong", "Correct", "Late"];
    case {'Valid'}
        Perfs = ["Wrong", "Correct"];
    otherwise
        Perfs = opts.Performance;
end

try
    fp = obj.SessionFP;
catch
    fp = obj.TaskFP;
end

time_bin_center = obj.("AngleHeadTrace"+AlignTo).TimePoints;

ind = obj.Ind;
these_trials = [];
for i = 1:length(Perfs)
    for j = 1:length(PortCorrect)
        N = sum(ind.("Target"+PortCorrect(j)) & ind.(Perfs(i)) & ind.(opts.Condition));
        ind_these = find(ind.("Target"+PortCorrect(j)) & ind.(Perfs(i)) & ind.(opts.Condition));
        if N < 5
            continue;
        elseif N > 50
            ind_these = ind_these(randperm(N, 50));
        end
        these_trials = [these_trials ind_these];
    end
end

these_trials = these_trials(randperm(length(these_trials)));
ind_time_bin_center = ismember(obj.("TimePoints"+AlignTo), time_bin_center);

for i = 1:length(these_trials)
    port_target = obj.PortCorrect(these_trials(i));
    this_color = GPSColor().("Port"+port_target);

    plot(ax, obj.("TimeFrom"+AlignTo){these_trials(i)}, obj.("AngleHead"){these_trials(i)}, ...
        'Color', this_color, 'LineStyle', '-', 'LineWidth', 1.5);
end

patch(ax, 'XData', [-2000 2000 2000 -2000], 'YData', [-100 -100 100 100], ...
    'FaceColor', 'w', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

switch AlignTo
    case {'In'}
        plot(ax, [-2000 2000], [0 0], ':', 'LineWidth', 1.2, 'Color', [.2 .2 .2 .8]);

        switch obj.Task
            case {'ThreeFPHoldWM'}
                patch(ax, 'XData', [0 250 250 0], 'YData', [-3 -3 3 3], 'FaceColor', GPSColor.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
            case {'ThreeFPHoldSRT'}
                patch(ax, 'XData', [0 fp*1000+300 fp*1000+300 0], 'YData', [-3 -3 3 3], 'FaceColor', GPSColor.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
            case {'ThreeFPHoldCRT'}
                patch(ax, 'XData', [fp*1000 fp*1000+200 fp*1000+200 fp*1000], 'YData', [-3 -3 3 3], 'FaceColor', GPSColor.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        end

        patch(ax, 'XData', [0 1000*fp 1000*fp 0], 'YData', [-100 -100 100 100], ...
            'FaceColor', 'm', 'FaceAlpha', 0.05, 'EdgeColor', 'none');

        set(ax, 'xlim', [-100 fp*1000+300], 'xtick', 0:500:2000, 'ylim', [-90 90], 'ytick', -90:30:90);

    case {'Out'}
        plot(ax, [-2000 2000], [0 0], ':', 'LineWidth', 1.2, 'Color', [.2 .2 .2 .8]);
        switch obj.Task
            case {'ThreeFPHoldSRT'}
                patch(ax, 'XData', [-fp*1000-100 200 200 -fp*1000-100], 'YData', [-3 -3 3 3], 'FaceColor', GPSColor.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
            case {'ThreeFPHoldCRT'}
                patch(ax, 'XData', [-200 0 0 -200], 'YData', [-3 -3 3 3], 'FaceColor', GPSColor.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        end

        plot(ax, [0 0], [-90 90], '-', 'LineWidth', 2, 'Color', [0 0 0 .7]);

        set(ax, 'xlim', [-fp*1000-100 300], 'xtick', -1500:500:500, 'ylim', [-90 90], 'ytick', -90:30:90);
end

if opts.Performance=="All"
    for j = 1:length(PortCorrect)

        if ~isfield(obj.("AngleHeadTrace"+AlignTo), PortCorrect(j)+"_"+opts.Condition)
            continue
        end

        M = obj.("AngleHeadTrace"+AlignTo).(PortCorrect(j)+"_"+opts.Condition);

        N = sum(ind.("Target"+PortCorrect(j)) & ind.(opts.Condition));
        if N < 5
            continue;
        end
        AngleHeadMean = cellfun(@(x) mean(x.ang, 'omitnan'), M);
        AngleHeadSTD  = cellfun(@(x) std(x.ang, 'omitnan'),  M);
        AngleHeadSEM  = AngleHeadSTD ./ sqrt(cellfun(@(x) length(x.ang), M));

        CI95 = tinv([0.025 0.975], N-1);
        AngleHeadCI95 = bsxfun(@times, AngleHeadSEM, CI95(:));

        AngleHeadMed  = cellfun(@(x) median(x.ang, 'omitnan'), M);
        AngleHeadQ1   = cellfun(@(x) prctile(x.ang, 25),  M);
        AngleHeadQ3   = cellfun(@(x) prctile(x.ang, 75),  M);

        AngleHeadMean = smoothdata(AngleHeadMean, "gaussian", 5);
        AngleHeadMed  = smoothdata(AngleHeadMed, "gaussian", 7);
        AngleHeadQ1   = smoothdata(AngleHeadQ1, "gaussian", 7);
        AngleHeadQ3   = smoothdata(AngleHeadQ3, "gaussian", 7);

        this_color = GPSColor().("Port"+PortCorrect(j));
        ls_this = "-";

        plot(ax, time_bin_center, AngleHeadMed, 'Color', this_color, 'LineStyle', ls_this, 'LineWidth', 3);
    end
else
    for i = 1:length(Perfs)
        for j = 1:length(PortCorrect)

            if ~isfield(obj.("AngleHeadTrace"+AlignTo), PortCorrect(j)+"_"+opts.Condition+"_"+Perfs(i))
                continue
            end

            M = obj.("AngleHeadTrace"+AlignTo).(PortCorrect(j)+"_"+opts.Condition+"_"+Perfs(i));

            N = sum(ind.("Target"+PortCorrect(j)) & ind.(opts.Condition) & ind.(Perfs(i)));
            if N < 5
                continue;
            end
            AngleHeadMean = cellfun(@(x) mean(x.ang, 'omitnan'), M);
            AngleHeadSTD  = cellfun(@(x) std(x.ang, 'omitnan'),  M);
            AngleHeadSEM  = AngleHeadSTD ./ sqrt(cellfun(@(x) length(x.ang), M));

            CI95 = tinv([0.025 0.975], N-1);
            AngleHeadCI95 = bsxfun(@times, AngleHeadSEM, CI95(:));

            AngleHeadMed  = cellfun(@(x) median(x.ang, 'omitnan'), M);
            AngleHeadQ1   = cellfun(@(x) prctile(x.ang, 25),  M);
            AngleHeadQ3   = cellfun(@(x) prctile(x.ang, 75),  M);

            AngleHeadMean = smoothdata(AngleHeadMean, "gaussian", 5);
            AngleHeadMed  = smoothdata(AngleHeadMed, "gaussian", 7);
            AngleHeadQ1   = smoothdata(AngleHeadQ1, "gaussian", 7);
            AngleHeadQ3   = smoothdata(AngleHeadQ3, "gaussian", 7);

            switch Perfs(i)
                case {'Premature'}
                    continue;
                case {'Correct'}
                    this_color = GPSColor().("Port"+PortCorrect(j));
                    ls_this = "-";
                case {'Wrong'}
                    this_color = GPSColor().("Port"+setdiff(obj.Ports, PortCorrect(j)));
                    ls_this = ":";
                case {'Late'}
                    continue;
            end
            %         patch(ax, 'XData', [time_bin_center flip(time_bin_center)], 'YData', [AngleHeadQ1 flip(AngleHeadQ3)], ...
            %             'FaceColor', this_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

            plot(ax, time_bin_center, AngleHeadMed, 'Color', this_color, 'LineStyle', ls_this, 'LineWidth', 3);
        end
    end
end

ax.XLabel.String     = "Time from poke "+AlignTo+" (ms)";
ax.XLabel.FontWeight = "Bold";
ax.YLabel.String     = "Head angle (°)";
ax.YLabel.FontWeight = "Bold";

end

