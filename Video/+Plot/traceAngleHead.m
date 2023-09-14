function traceAngleHead(ax, obj, ForePeriod, AlignTo, opts)

arguments
    ax  (1, 1) handle
    obj (1, 1) 
    ForePeriod (1, 1) {mustBeMember(ForePeriod, ["Short", "Med", "Long"])} = "Med"
    AlignTo (1, 1) {mustBeMember(AlignTo, ["In", "Out"])} = "In"
    opts.color (1, 1) GPSColor = GPSColor()
    opts.PortChosen (1, 1) {mustBeMember(opts.PortChosen, ["L", "R", "Both"])} = "Both"
    opts.Performance (1, 1) {mustBeMember(opts.Performance, ["Correct", "Wrong", "Valid"])} = "Correct"
    opts.Label (1, 1) {mustBeMember(opts.Label, ["All", "None", "Saline", "Control", "Chemo"])} = "All"
end

switch opts.PortChosen
    case {'Both'}
        PortChosen = ["L", "R"];
    otherwise
        PortChosen = opts.PortChosen;
end

switch opts.Performance
    case {'Valid'}
        Perfs = ["Wrong", "Correct"];
    otherwise
        Perfs = opts.Performance;
end

switch ForePeriod
    case {'Short'}
        fp = 0.5;
    case {'Med'}
        fp = 1.0;
    case {'Long'}
        fp = 1.5;
end

time_bin_center = obj.("AngleHeadTrace"+AlignTo).("TimePoints_"+ForePeriod);

ind = obj.Ind;
these_trials = [];
for i = 1:length(Perfs)
    for j = 1:length(PortChosen)
        N = sum(ind.("Choose"+PortChosen(j)) & ind.(Perfs(i)) & ind.(ForePeriod));
        ind_these = find(ind.("Choose"+PortChosen(j)) & ind.(Perfs(i)) & ind.(ForePeriod));
        if N < 5
            continue;
        elseif N > 30
            ind_these = ind_these(randperm(N, 30));
        end
        these_trials = [these_trials ind_these];
    end
end
these_trials = these_trials(randperm(length(these_trials)));
ind_time_bin_center = ismember(obj.("TimePoints"+AlignTo), time_bin_center);
for i = 1:length(these_trials)
    perf_this = obj.Performance(these_trials(i));
    port_chosen = obj.PortChosen(these_trials(i));
    switch perf_this
        case {'Correct'}
            this_color = GPSColor().("Port"+port_chosen);
            ls_this = "-";
        case {'Wrong'}
            this_color = GPSColor().("Port"+setdiff(obj.Ports, port_chosen));
            ls_this = ":";
    end
    plot(ax, obj.("TimeFrom"+AlignTo){these_trials(i)}, obj.("AngleHead"){these_trials(i)}, ...
        'Color', this_color, 'LineStyle', ls_this, 'LineWidth', 1.5);
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

for i = 1:length(Perfs)

    for j = 1:length(PortChosen)

        M = obj.("AngleHeadTrace"+AlignTo).(PortChosen(j)+"_"+ForePeriod+"_"+Perfs(i));

        N = sum(ind.("Choose"+PortChosen(j)) & ind.(Perfs(i)) & ind.(ForePeriod));
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
            case {'Correct'}
                this_color = GPSColor().("Port"+PortChosen(j));
                ls_this = "-";
            case {'Wrong'}
                this_color = GPSColor().("Port"+setdiff(obj.Ports, PortChosen(j)));
                ls_this = ":";
        end
%         patch(ax, 'XData', [time_bin_center flip(time_bin_center)], 'YData', [AngleHeadQ1 flip(AngleHeadQ3)], ...
%             'FaceColor', this_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        plot(ax, time_bin_center, AngleHeadMed, 'Color', this_color, 'LineStyle', ls_this, 'LineWidth', 3);
    end
end

if ~isempty(obj.("AngleHeadTrace"+AlignTo+"Test"))

    % 1. to discriminant between chosen port,
    % between correct trials (chose L) and correct trials (chose R)
    cc_table = obj.("AngleHeadTrace"+AlignTo+"Test").(ForePeriod+"_CC");
    if ~isempty(cc_table)
        sig_cc_beg = cc_table.time_points(find(diff([0; smoothdata(cc_table.sig, "movmean", 3)==1])==1));
        sig_cc_end = cc_table.time_points(find(diff([smoothdata(cc_table.sig, "movmean", 3)==1; 0])==-1));

        for s = 1:length(sig_cc_beg)
            patch(ax, 'XData', [sig_cc_beg(s) sig_cc_end(s) sig_cc_end(s) sig_cc_beg(s)], 'YData', [58.5 58.5 60 60], ...
                'FaceColor', 'k', 'EdgeColor', 'k');
        end
    end

%     % 2. to discriminant between chosen port,
%     % 2.1.
%     cw_lr_table = obj.("AngleHeadTrace"+AlignTo+"Test").(ForePeriod+"_CW_LR");
%     if ~isempty(cw_lr_table)
%         sig_cw_lr_beg = cw_lr_table.time_points(find(diff([0; smoothdata(cw_lr_table.sig, "movmean", 3)==1])==1));
%         sig_cw_lr_end = cw_lr_table.time_points(find(diff([smoothdata(cw_lr_table.sig, "movmean", 3)==1; 0])==-1));
% 
%         for s = 1:length(sig_cw_lr_beg)
%             patch(ax, 'XData', [sig_cw_lr_beg(s) sig_cw_lr_end(s) sig_cw_lr_end(s) sig_cw_lr_beg(s)], 'YData', [56 56 57.5 57.5], ...
%                 'FaceColor', 'w', 'EdgeColor', GPSColor.PortL);
%         end
%     end
% 
%     % 2.2. between correct trials (chose R) and wrong trials (chose L)
%     cw_rl_table = obj.("AngleHeadTrace"+AlignTo+"Test").(ForePeriod+"_CW_RL");
%     if ~isempty(cw_rl_table)
%         sig_cw_rl_beg = cw_rl_table.time_points(find(diff([0; smoothdata(cw_rl_table.sig, "movmean", 3)==1])==1));
%         sig_cw_rl_end = cw_rl_table.time_points(find(diff([smoothdata(cw_rl_table.sig, "movmean", 3)==1; 0])==-1));
% 
%         for s = 1:length(sig_cw_rl_beg)
%             patch(ax, 'XData', [sig_cw_rl_beg(s) sig_cw_rl_end(s) sig_cw_rl_end(s) sig_cw_rl_beg(s)], 'YData', [53.5 53.5 55 55], ...
%                 'FaceColor', 'w', 'EdgeColor', GPSColor.PortR, 'LineWidth', 1);
%         end
%     end
% 
%     % 3. to discriminant between correct port,
%     % 3.1. between correct trials (chose L) and wrong trials (chose L)
%     cw_ll_table = obj.("AngleHeadTrace"+AlignTo+"Test").(ForePeriod+"_CW_LL");
%     if ~isempty(cw_ll_table)
%         sig_cw_ll_beg = cw_ll_table.time_points(find(diff([0; smoothdata(cw_ll_table.sig, "movmean", 3)==1])==1));
%         sig_cw_ll_end = cw_ll_table.time_points(find(diff([smoothdata(cw_ll_table.sig, "movmean", 3)==1; 0])==-1));
% 
%         for s = 1:length(sig_cw_ll_beg)
%             patch(ax, 'XData', [sig_cw_ll_beg(s) sig_cw_ll_end(s) sig_cw_ll_end(s) sig_cw_ll_beg(s)], 'YData', -[56 56 57.5 57.5], ...
%                 'FaceColor', GPSColor.PortL, 'EdgeColor', GPSColor.PortL, 'LineWidth', 1);
%         end
%     end
% 
%     % 3.2. between correct trials (chose R) and wrong trials (chose R)
%     cw_rr_table = obj.("AngleHeadTrace"+AlignTo+"Test").(ForePeriod+"_CW_RR");
%     if ~isempty(cw_rr_table)
%         sig_cw_rr_beg = cw_rr_table.time_points(find(diff([0; smoothdata(cw_rr_table.sig, "movmean", 3)==1])==1));
%         sig_cw_rr_end = cw_rr_table.time_points(find(diff([smoothdata(cw_rr_table.sig, "movmean", 3)==1; 0])==-1));
% 
%         for s = 1:length(sig_cw_rr_beg)
%             patch(ax, 'XData', [sig_cw_rr_beg(s) sig_cw_rr_end(s) sig_cw_rr_end(s) sig_cw_rr_beg(s)], 'YData', -[53.5 53.5 55 55], ...
%                 'FaceColor', GPSColor.PortR, 'EdgeColor', GPSColor.PortR, 'LineWidth', 1);
%         end
%     end

end

ax.XLabel.String     = "Time from poke "+AlignTo+" (ms)";
ax.XLabel.FontWeight = "Bold";
ax.YLabel.String     = "Head angle (°)";
ax.YLabel.FontWeight = "Bold";

end

