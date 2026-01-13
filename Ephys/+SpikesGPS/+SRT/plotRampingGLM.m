function fig_ramp = plotRampingGLM(r, id, StatGLM, order)

rs = struct();
rs.Units = r.Units;
rs.EphysTable = r.EphysTable;
rs.PopPSTH = r.PopPSTH;

switch length(id)
    case 1
        ku = id;
        uid = rs.Units.SpikeNotes(ku, [1 2]);
    case 2
        uid = id;
        ku = find(rs.Units.SpikeNotes(:,1)==id(1) & rs.Units.SpikeNotes(:,2)==id(2));
end

p_all = [StatGLM.ramping.pval_pre; StatGLM.ramping.pval_post];
alpha_adj = BH_method(p_all, .05);

stat_tbl = StatGLM.ramping(StatGLM.ramping.unit_id==ku, :);

SpikeTimes = rs.Units.SpikeTimes(ku);
EphysTable = rs.EphysTable;
PopPSTH = rs.PopPSTH;

EphysTable = EphysTable(ismember(EphysTable.Outcome, "Correct"), :);
num_trial = height(EphysTable);

FPs = [.5 1 1.5];
Ports = [1 2];
num_fp = length(FPs);
num_port = length(Ports);

plt = GPSPlot();
c = GPSColor();


fig_ramp = figure(68); clf(fig_ramp);
set_fig_default(fig_ramp);
set(fig_ramp, 'Name', 'Ramping', 'Position', [5 5 8.5 9.25]);

% PSTH
ax_psth_bg = plt.assign_ax_to_fig(fig_ramp, 3, 1, [1 1 2.5 3.2], [2.5 1]);
cellfun(@(x) set(x, 'XLim', [-500 2000], 'YLim', [0 1], 'XColor', 'none', 'YColor', 'none'), ax_psth_bg);
ax_psth = plt.assign_ax_to_fig(fig_ramp, 3, 1, [1 1 2.5 3.2], [2.5 1]);
cellfun(@(x) set(x, 'XLim', [-500 2000], 'Color', 'none'), ax_psth);
cellfun(@(x) set(x, 'XTickLabel', [], 'YTickLabel', []), ax_psth(1:2));

for fp = 1:num_fp
    ax = ax_psth_bg{fp};
    fill(ax, 1000*[0 FPs(fp) FPs(fp) 0], [0 0 1 1], 'r', 'FaceColor', c.FP, 'FaceAlpha', .15, 'EdgeColor', 'none');
    ax = ax_psth{fp};
    for p = 1:num_port
        t_psth = PopPSTH.CentIn{fp,p}(1,:);
        psth   = PopPSTH.CentIn{fp,p}(ku+1,:);

        id_this = find(t_psth>=-500 & t_psth<FPs(fp)*1000+500);
        t_psth = t_psth(id_this);
        psth   = psth(id_this);

        plot(ax, t_psth, psth, 'Color', c.EphysDir{3,p}, 'LineWidth', 1);
    end
end
y_lim = cellfun(@(x) x.YLim(2), ax_psth);
y_lim = max(y_lim, [], 'all');
cellfun(@(x) set(x, 'YLim', [0 y_lim]), ax_psth);

xlabel(ax, 'Time to Cent-In (ms)');
ylabel(ax, 'Firing rate (Hz)');

% raster
h_raster = 0.03 * num_trial;
ax_raster = plt.assign_ax_to_fig(fig_ramp, 1, 1, [1 4.3 2.5 h_raster], [2.5 h_raster]);
ax = ax_raster{1};
set(ax, 'YDir', 'reverse', 'XLim', [-500 2000], 'YLim', [.5 num_trial+.5], 'XColor', 'none', 'YColor', 'none');

t_raster = -500:1:2000;
raster_mat = nan(num_trial, length(t_raster));
trial_count = 0;

for fp = 1:num_fp
    for p = 1:num_port
        e_tbl = EphysTable(EphysTable.FP==FPs(fp) & EphysTable.PortCorrect==Ports(p), :);
        switch order
            case "Sorted"
                e_tbl = sortrows(e_tbl, "RT");
            case "ByTime"
        end
        xx_all = []; % spike times
        yy_all = [];
        xxrt_all = []; % reaction (cent-out) times
        yyrt_all = [];
        for i = 1:height(e_tbl)
            yy = [-.4 .4] + i + trial_count;
            t_spk = SpikeTimes.timings - e_tbl.tCentIn(i);
            t_spk = t_spk(t_spk>=-500 & t_spk<FPs(fp)*1000+500);
            if ~isempty(t_spk)
                for i_xx = 1:length(t_spk)
                    xx_all = [xx_all, t_spk(i_xx), t_spk(i_xx), nan];
                    yy_all = [yy_all, yy, nan];
                end
                xxrt_all = [xxrt_all; e_tbl.tCentOut(i) - e_tbl.tCentIn(i)];
                yyrt_all = [yyrt_all; i + trial_count];
            end
        end
        fill(ax, 1000*[0 FPs(fp) FPs(fp) 0], trial_count+.5+[0 0 1 1]*height(e_tbl), 'r', 'FaceColor', c.FP, 'FaceAlpha', .15, 'EdgeColor', 'none');
        line(ax, xx_all, yy_all, 'color', c.EphysDir{3,p}, 'linewidth', .2);
        scatter(ax, xxrt_all, yyrt_all, 4, c.EphysCentOut, 'filled', 'o', 'MarkerFaceAlpha', .5);

        trial_count = trial_count + height(e_tbl);
    end
end

% GLM
ax_ramp = plt.assign_ax_to_fig(fig_ramp, 1, 1, [5 1 2 1.5], [2.5 1.5]);
ax = ax_ramp{1};
set(ax, 'XLim', [-50 1050]);

bin_width = 0.1;
time_bin = 0:bin_width:1;
time_pos = movmean(time_bin, 2, 'Endpoints', 'discard');
time_sep = 0.5;

t_pred = table(time_pos', categorical(time_pos'<time_sep, [true false], ["Pre", "Post"]), 'VariableNames', {'Time', 'Period'});

for p = 1:num_port
    stat_this = stat_tbl(stat_tbl.port==p, :);
    if stat_this.status~="glm_fitted"
        continue
    end
    if stat_this.pval_pre < alpha_adj
        ls_pre = '-';
    else
        ls_pre = '--';
    end
    if stat_this.pval_post < alpha_adj
        ls_post = '-';
    else
        ls_post = '--';
    end
    [c_pred, c_pred_ci] = stat_this.mdl{1}.predict(t_pred);

    id_pre = t_pred.Period=="Pre";
    fill(ax, [t_pred.Time(id_pre); flip(t_pred.Time(id_pre))]*1000, [c_pred_ci(id_pre, 1); flip(c_pred_ci(id_pre, 2))], 'r', 'FaceColor', c.EphysDir{3,p}, 'FaceAlpha', .2, 'EdgeColor', 'none');
    id_post = t_pred.Period=="Post";
    fill(ax, [t_pred.Time(id_post); flip(t_pred.Time(id_post))]*1000, [c_pred_ci(id_post, 1); flip(c_pred_ci(id_post, 2))], 'r', 'FaceColor', c.EphysDir{3,p}, 'FaceAlpha', .2, 'EdgeColor', 'none');
    
    plot(ax, t_pred.Time(t_pred.Period=="Pre")*1000, c_pred(t_pred.Period=="Pre"), 'LineStyle', ls_pre, 'Color', c.EphysDir{3,p}, 'LineWidth', 1);
    plot(ax, t_pred.Time(t_pred.Period=="Post")*1000, c_pred(t_pred.Period=="Post"), 'LineStyle', ls_post, 'Color', c.EphysDir{3,p}, 'LineWidth', 1);
end
xlabel(ax, 'Time to Cent-In (ms)');
ylabel(ax, 'Pred. firing rate (Hz)');

ax_ramp_lb = plt.assign_ax_to_fig(fig_ramp, 1, 1, [5 2.75 2.5 1], [2.5 1]);
ax = ax_ramp_lb{1};
set(ax, 'XLim', [-50 1050], 'YLim', [0 10], 'XColor', 'none', 'YColor', 'none');

txt_f = [7.5 2];
txt_p = [5 0];
for p = 1:num_port
    stat_this = stat_tbl(stat_tbl.port==p, :);
    if stat_this.status~="glm_fitted"
        continue
    end
    text(ax, 250, txt_f(p), sprintf('fold = %.2f', stat_this.fold_pre), 'Color', c.EphysDir{3,p}, 'FontSize', 6, 'HorizontalAlignment', 'center');
    if stat_this.pval_pre < .001
        text(ax, 250, txt_p(p), 'p < .001', 'Color', c.EphysDir{3,p}, 'FontSize', 6, 'HorizontalAlignment', 'center');
    else
        text(ax, 250, txt_p(p), sprintf('p = %.3f', stat_this.pval_pre), 'Color', c.EphysDir{3,p}, 'FontSize', 6, 'HorizontalAlignment', 'center');
    end

    text(ax, 750, txt_f(p), sprintf('fold = %.2f', stat_this.fold_post), 'Color', c.EphysDir{3,p}, 'FontSize', 6, 'HorizontalAlignment', 'center');
    if stat_this.pval_post < .001
        text(ax, 750, txt_p(p), 'p < .001', 'Color', c.EphysDir{3,p}, 'FontSize', 6, 'HorizontalAlignment', 'center');
    else
        text(ax, 750, txt_p(p), sprintf('p = %.3f', stat_this.pval_post), 'Color', c.EphysDir{3,p}, 'FontSize', 6, 'HorizontalAlignment', 'center');
    end
end

% auto-correlation
ax_ac = plt.assign_ax_to_fig(fig_ramp, 1, 1, [5.5 4.75 1.5 1.5], [1.5 1.5]);
ax = ax_ac{1};

[auto_cor, lags] = computeAutoCorr(SpikeTimes.timings, 50, 1);
bar(ax, lags, auto_cor, 'FaceColor', 'k');
xlabel(ax, 'Lags (ms)');
ylabel(ax, 'Count');

% waveform
ax_wave = plt.assign_ax_to_fig(fig_ramp, 1, 1, [5.5 6.5 1.5 1.5], [1.5 1.5]);
ax = ax_wave{1};
set(ax, 'XColor', 'none', 'YColor', 'none', 'XLim', [-32 32]);

t_wave = -31:1:32;
wave_m = mean(SpikeTimes.wave, 1);
% wave_ci = bootci(1000, {@mean, SpikeTimes.wave}, 'Alpha', .01);
% fill(ax, [t_wave'; flip(t_wave')], [wave_ci(1,:)'; flip(wave_ci(2,:)')], 'k', 'FaceColor', 'k', 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(ax, t_wave, wave_m, '-k', 'LineWidth', 1);

switch rs.Units.SpikeNotes(ku, 3)
    case 1
        title(sprintf('#%d (Ch %d | Unit %d | SU)', ku, uid(1), uid(2)), 'fontsize', 7);
    case 2
        title(sprintf('#%d (Ch %d | Unit %d | MU)', ku, uid(1), uid(2)), 'fontsize', 7);
    otherwise
end

h_ax = sum(ax_raster{1}.Position([2 4]));
h_fig = max([h_ax+.75, 9.25]);
fig_ramp.Position(4) = h_fig;

%
set_fig_title(fig_ramp, sprintf('%s | %s | %s', EphysTable.Subject(1), EphysTable.Session(1), order));

fig_name = sprintf('Ramping_%s_%s_Ch%d_Unit%d_%s', EphysTable.Subject(1), EphysTable.Session(1), uid(1), uid(2), order);
fig_folder = fullfile('Figure', 'Ramping');
if ~isfolder(fig_folder)
    mkdir(fig_folder);
end
fig_path = fullfile(fig_folder, fig_name);
exportgraphics(fig_ramp, sprintf('%s.png', fig_path), 'Resolution', 300);
exportgraphics(fig_ramp, sprintf('%s.pdf', fig_path), 'ContentType', 'vector');

end % plotRampingGLM

