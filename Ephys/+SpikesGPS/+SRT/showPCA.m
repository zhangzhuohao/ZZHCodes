function [fig_pc, Coef, PCs] = showPCA(r)

c = cell(3,2);
c{1,1} = [246 178 147] / 255;
c{2,1} = [220 109 87]  / 255;
c{3,1} = [183 34 48]   / 255;
c{1,2} = [182 215 232] / 255;
c{2,2} = [109 173 209] / 255;
c{3,2} = [49 124 183]  / 255;

plt = GPSPlot();

% cluster spikes by waveform, get regular good spikes
regular_file = dir('RegularIndex*.mat');
if ~isempty(regular_file)
    load(regular_file.name, 'RegularIndex');
else
    RegularIndex = clusterSpikes(r, 'PeakTrough');
end
% get units location
location_file = dir('UnitLocation*.mat');
if ~isempty(location_file)
    load(location_file.name, 'unitLocation');
else
    unitLocation = getUnitLocation(r);
end
loc_unit = unitLocation.where(:, RegularIndex.Good);
k_unit   = r.ChanMap.kcoords(r.Units.SpikeNotes(RegularIndex.Good, 1));
x_unit   = r.ChanMap.xcoords(r.Units.SpikeNotes(RegularIndex.Good, 1));

shanks  = unique(k_unit);
n_shank = length(shanks);
x_k = cell(1,n_shank);
y_k = cell(1,n_shank);
for i = 1:n_shank
   x_k{i} = unique(x_unit(k_unit==shanks(i)));
   y_k{i} = [min(r.ChanMap.ycoords(r.ChanMap.kcoords==shanks(i))) max(r.ChanMap.ycoords(r.ChanMap.kcoords==shanks(i)))];
end

k_sep = 150;
y_range = [min(cell2mat(y_k)) max(cell2mat(y_k))];


%% show warped population activity (FP pooled)
fr_pool = [r.popSDFWarped.popWarpPooled.sdf.centin; r.popSDFWarped.popWarpPooled.sdf.trigger];
t_pool  = [r.popSDFWarped.popWarpPooled.t_warp.centin; r.popSDFWarped.popWarpPooled.t_warp.trigger];
fr_pool = cellfun(@(x) x(RegularIndex.Good,:), fr_pool, 'UniformOutput', false);

n_u = size(fr_pool{1}, 1); % number of spikes in total
n_t_pool = cellfun(@(x) size(x, 2), fr_pool); % time length in each conditions

% Do soft normalization
fr_pool_c = cell2mat(fr_pool(:)'); % concatenate all conditions
fr_pool_c_norm = nan(size(fr_pool_c)); % initialize soft-norm output
for i = 1:n_u
    fr_pool_c_norm(i,:) = cal_soft_norm(fr_pool_c(i,:));
end
fr_pool_norm = mat2cell(fr_pool_c_norm, n_u, n_t_pool(:)'); % split conditions
fr_pool_norm = reshape(fr_pool_norm, size(fr_pool));


%% plot heatmap
% regular units only
ind_sort = cellfun(@(x) x(ismember(x, RegularIndex.Good)), r.popSDFWarped.popWarpPooled.IndSort, 'UniformOutput', false);
[~,~,ind_sort] = cellfun(@(x) unique(x, 'sorted'), ind_sort, 'UniformOutput', false);

fig_pc = figure(88); clf(fig_pc);
set_fig_default(fig_pc);
set(fig_pc, 'Name', 'PCA', 'Position', [5 5 18 14.75]);

% 1250 1250
ax_pop_centin = plt.assign_ax_to_fig(fig_pc, 1, 2, [1 10.5 7 3], [2 3]);
cellfun(@(x) set(x, 'XLim', [-1000 250], 'YLim', [.5 n_u+0.5], 'YDir', 'reverse', 'CLim', [-.2 .6]), ax_pop_centin);
cellfun(@(x) xline(x, 0, '-w', LineWidth=1), ax_pop_centin);
set(ax_pop_centin{2}, 'YTickLabel', []);
cellfun(@(x) colormap(x, 'parula'), ax_pop_centin);

ax_pop_trig = plt.assign_ax_to_fig(fig_pc, 1, 2, [3.1 10.5 7 3], [2 3]);
cellfun(@(x) set(x, 'XLim', [-250 1000], 'YLim', [.5 n_u+0.5], 'YDir', 'reverse', 'CLim', [-.2 .6], 'YTick', []), ax_pop_trig);
cellfun(@(x, t) xline(x, t(1:3), '-w', LineWidth=1), ax_pop_trig, r.popSDFWarped.popWarpPooled.t_points.trigger);
cellfun(@(x) colormap(x, 'parula'), ax_pop_trig);

% sorted by contra
sortby = 1;
side_label = {'Contra', 'Ipsi'};
for p = 1:2
    title(ax_pop_centin{p}, side_label{p}, 'Color', c{3,p});
    imagesc(ax_pop_centin{p}, t_pool{1,p}, 1:n_u, fr_pool_norm{1,p}(ind_sort{sortby}, :));
    imagesc(ax_pop_trig{p}, t_pool{2,p}, 1:n_u, fr_pool_norm{2,p}(ind_sort{sortby}, :));

    if p==1
        xlabel(ax_pop_centin{p}, 'Time to Cent-In (ms)');
        ylabel(ax_pop_centin{p}, 'Unit');
        xlabel(ax_pop_trig{p}, ' / Trigger (ms)');
    end
end
cb = colorbar(ax_pop_trig{2}, 'Units', 'centimeters', 'Position', [10.6 10.5 .25 3], 'TickDirection', 'out', 'FontSize', 7);
cb.Label.String = "Soft-norm. activity";
cb.Label.FontSize = 7;


%% PCA on different phases
% Do soft normalization
fr = r.popSDFWarped.sdf;
fr = cellfun(@(x) x(RegularIndex.Good,:), fr, 'UniformOutput', false);

fr_c = cell2mat(fr(:)'); % concatenate all conditions
fr_c_norm = nan(size(fr_c)); % initialize soft-norm output
for i = 1:n_u
    fr_c_norm(i,:) = cal_soft_norm(fr_c(i,:));
end

n_t = cellfun(@(x) size(x, 2), fr); % time length in each conditions

fr_norm = mat2cell(fr_c_norm, n_u, n_t(:)'); % split conditions
fr_norm = reshape(fr_norm, 3, 2);

% time information
t_warp = r.popSDFWarped.t_warp;
t_point = r.popSDFWarped.t_points;

% analyze from 1500ms pre cent-in, to 200ms post choice-in
fr_norm = cellfun(@(x, t_w, t_p) x(:, t_w>-1000 & t_w<=t_p(4)+200), fr_norm, t_warp, t_point, 'UniformOutput', false);
t_warp = cellfun(@(t_w, t_p) t_w(t_w>-1000 & t_w<=t_p(4)+200), t_warp, t_point, 'UniformOutput', false);

fr_c_norm = cell2mat(fr_norm(:)');
n_u = size(fr_norm{1}, 1); % number of units in total
n_t = cellfun(@(x) size(x, 2), fr_norm); % time length in each conditions

% covariance matrxi and total var
cov_tot = cov(fr_c_norm');
var_tot = trace(cov_tot);

%% Perform PCA using activity in different phases

% coef axis
coef_thres = 0.3;
ax_coef_1 = plt.assign_ax_to_fig(fig_pc, 4, 1, [1 1 3 8], [3 1.5]);
cellfun(@(x) set(x, 'XLim', [50 k_sep*n_shank+100], 'YLim', y_range+[-25 25], 'XColor', 'none', 'YColor', 'none'), ax_coef_1);
ax_coef_2 = plt.assign_ax_to_fig(fig_pc, 4, 1, [9.5 1 3 8], [3 1.5]);
cellfun(@(x) set(x, 'XLim', [50 k_sep*n_shank+100], 'YLim', y_range+[-25 25], 'XColor', 'none', 'YColor', 'none'), ax_coef_2);

% pc axis
ax_hold = plt.assign_ax_to_fig(fig_pc, 4, 1, [4.75 1 2.8 8], [2.8 1.5]); cellfun(@(x) set(x, 'XLim', [-1000 1800]), ax_hold);
cellfun(@(x) yline(x, 0, ':k', LineWidth=1), ax_hold); cellfun(@(x) xline(x, 0:500:1500, ':k', LineWidth=1), ax_hold);
cellfun(@(x) ylabel(x, 'PC-1'), ax_hold);
ax_move = plt.assign_ax_to_fig(fig_pc, 4, 1, [7.75 1 0.8 8], [0.8 1.5]); cellfun(@(x) set(x, 'XLim', [0 800]), ax_move);
cellfun(@(x) yline(x, 0, ':k', LineWidth=1), ax_move);

ax_hold_2 = plt.assign_ax_to_fig(fig_pc, 4, 1, [13.25 1 2.8 8], [2.8 1.5]); cellfun(@(x) set(x, 'XLim', [-1000 1800]), ax_hold_2);
cellfun(@(x) yline(x, 0, ':k', LineWidth=1), ax_hold_2); cellfun(@(x) xline(x, 0:500:1500, ':k', LineWidth=1), ax_hold_2);
cellfun(@(x) ylabel(x, 'PC-2'), ax_hold_2);
ax_move_2 = plt.assign_ax_to_fig(fig_pc, 4, 1, [16.25 1 0.8 8], [0.8 1.5]); cellfun(@(x) set(x, 'XLim', [0 800]), ax_move_2);
cellfun(@(x) yline(x, 0, ':k', LineWidth=1), ax_move_2);

ax_var = plt.assign_ax_to_fig(fig_pc, 4, 1, [4.75 1 2.8 8], [2.8 1.5]);
ax_var_2 = plt.assign_ax_to_fig(fig_pc, 4, 1, [13.25 1 2.8 8], [2.8 1.5]);
cellfun(@(x) set(x, 'XColor', 'none', 'YColor', 'none', 'Color', 'none', 'YLim', [0 100], 'XLim', [0 100]), ax_var);
cellfun(@(x) set(x, 'XColor', 'none', 'YColor', 'none', 'Color', 'none', 'YLim', [0 100], 'XLim', [0 100]), ax_var_2);

% approach
fr_norm_app = cellfun(@(x, t_w, t_p) x(:, t_w>=t_p(1)-1000 & t_w<t_p(1)), fr_norm, t_warp, t_point, 'UniformOutput', false);

fr_c_norm_app = cell2mat(fr_norm_app(:)'); % concatenate

[coef_app, ~, ~, ~, explain_app] = pca(fr_c_norm_app');
var_app_tot = diag(coef_app' * cov_tot * coef_app);
explain_app_tot = 100 * var_app_tot ./ var_tot;

pc_c_app = (fr_c_norm'*coef_app)';
pc_app = mat2cell(pc_c_app, size(pc_c_app,1), n_t(:)');
pc_app = reshape(pc_app, 3, 2);
% fix the direction of pc axis
for i = 1:size(coef_app,2)
    d_pc = cellfun(@(x, t_w, t_p) x(i,t_w==(t_p(1)-500))-x(i,t_w==t_p(1)), pc_app, t_warp, t_point);
    if mean(d_pc, 'all')<0
        coef_app(:,i) = -coef_app(:,i);
    end
end
pc_app = cellfun(@(x) (x'*coef_app)', fr_norm, 'UniformOutput', false);


% plot unit location with projection weight
id = 1;

w_app = normalize(abs(coef_app), 1, 'range');
w_app = w_app .* sign(coef_app);
for k = 1:n_shank
    x_shift = -mean(x_k{k})+k_sep*k;
    plot(ax_coef_1{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_1{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_app(:,1))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_1{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_app(id_this, 1)), 'k');

    text(ax_coef_1{id}, mean(x_k{k})+x_shift, y_range(2)+100, sprintf('shank%d', shanks(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);

    plot(ax_coef_2{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_2{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_app(:,2))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_2{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_app(id_this, 2)), 'k');
end

% plot pc-1 and pc-2
text(ax_var{id}, 5, 100, sprintf('%.1f%%', explain_app(1)), 'FontSize', 6);
text(ax_var{id}, 5, 85, sprintf('%.1f%%', explain_app_tot(1)), 'FontSize', 6);
text(ax_var_2{id}, 5, 100, sprintf('%.1f%%', explain_app(2)), 'FontSize', 6);
text(ax_var_2{id}, 5, 85, sprintf('%.1f%%', explain_app_tot(2)), 'FontSize', 6);
for fp = 1:3
    for p = 1:2
        t_w = t_warp{fp, p};
        t_p = t_point{fp, p};
        pc_this = pc_app{fp, p};

        ax = ax_hold{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_hold_2{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move_2{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);
    end
end
y_hold = ax_hold{id}.YLim;
y_move = ax_move{id}.YLim;

y_lim = [min([y_hold, y_move]), max([y_hold, y_move])];
set(ax_hold{id}, 'YLim', y_lim);
set(ax_move{id}, 'YLim', y_lim, 'YColor', 'none');

y_hold_2 = ax_hold_2{id}.YLim;
y_move_2 = ax_move_2{id}.YLim;

y_lim_2 = [min([y_hold_2, y_move_2]), max([y_hold_2, y_move_2])];
set(ax_hold_2{id}, 'YLim', y_lim_2);
set(ax_move_2{id}, 'YLim', y_lim_2, 'YColor', 'none');

text(ax_hold{id}, -5200, mean(ax_hold{id}.YLim), 'Approach', 'Rotation', 90, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'FontSize', 8, 'HorizontalAlignment', 'center');

% hold
fr_norm_hold = cellfun(@(x, t_w, t_p) x(:, t_w>=t_p(1) & t_w<t_p(2)), fr_norm, t_warp, t_point, 'UniformOutput', false);
fr_c_norm_hold = cell2mat(fr_norm_hold(:)'); % concatenate

[coef_hold, ~, ~, ~, explain_hold] = pca(fr_c_norm_hold');
var_hold_tot = diag(coef_hold' * cov_tot * coef_hold);
explain_hold_tot = 100 * var_hold_tot ./ var_tot;

pc_c_hold = (fr_c_norm'*coef_hold)';
pc_hold = mat2cell(pc_c_hold, size(pc_c_hold,1), n_t(:)');
pc_hold = reshape(pc_hold, 3, 2);
% fix the direction of pc axis
for i = 1:size(coef_hold,2)
    d_pc = cellfun(@(x, t_w, t_p) x(i,t_w==(t_p(1)+500))-x(i,t_w==t_p(1)), pc_hold, t_warp, t_point);
    if mean(d_pc, 'all')<0
        coef_hold(:,i) = -coef_hold(:,i);
    end
end
pc_hold = cellfun(@(x) (x'*coef_hold)', fr_norm, 'UniformOutput', false);

% plot unit location with projection weight
id = 2;
w_hold = normalize(abs(coef_hold), 1, 'range');
w_hold = w_hold .* sign(coef_hold);
for k = 1:n_shank
    x_shift = -mean(x_k{k})+k_sep*k;
    plot(ax_coef_1{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_1{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_hold(:,1))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_1{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_hold(id_this, 1)), 'k');

    plot(ax_coef_2{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_2{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_hold(:,2))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_2{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_hold(id_this, 2)), 'k');
end

% plot pc-1 and pc-2
text(ax_var{id}, 5, 100, sprintf('%.1f%%', explain_hold(1)), 'FontSize', 6);
text(ax_var{id}, 5, 85, sprintf('%.1f%%', explain_hold_tot(1)), 'FontSize', 6);
text(ax_var_2{id}, 5, 100, sprintf('%.1f%%', explain_hold(2)), 'FontSize', 6);
text(ax_var_2{id}, 5, 85, sprintf('%.1f%%', explain_hold_tot(2)), 'FontSize', 6);
for fp = 1:3
    for p = 1:2
        t_w = t_warp{fp, p};
        t_p = t_point{fp, p};
        pc_this = pc_hold{fp, p};

        ax = ax_hold{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_hold_2{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move_2{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);
    end
end
y_hold = ax_hold{id}.YLim;
y_move = ax_move{id}.YLim;

y_lim = [min([y_hold, y_move]), max([y_hold, y_move])];
set(ax_hold{id}, 'YLim', y_lim);
set(ax_move{id}, 'YLim', y_lim, 'YColor', 'none');

y_hold_2 = ax_hold_2{id}.YLim;
y_move_2 = ax_move_2{id}.YLim;

y_lim_2 = [min([y_hold_2, y_move_2]), max([y_hold_2, y_move_2])];
set(ax_hold_2{id}, 'YLim', y_lim_2);
set(ax_move_2{id}, 'YLim', y_lim_2, 'YColor', 'none');

text(ax_hold{id}, -5200, mean(ax_hold{id}.YLim), 'Hold', 'Rotation', 90, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'FontSize', 8, 'HorizontalAlignment', 'center');

% response
fr_norm_resp = cellfun(@(x, t_w, t_p) x(:, t_w>=t_p(2) & t_w<t_p(3)), fr_norm, t_warp, t_point, 'UniformOutput', false);
fr_c_norm_resp = cell2mat(fr_norm_resp(:)'); % concatenate

[coef_resp, ~, ~, ~, explain_resp] = pca(fr_c_norm_resp');
var_resp_tot = diag(coef_resp' * cov_tot * coef_resp);
explain_resp_tot = 100 * var_resp_tot ./ var_tot;

pc_c_resp = (fr_c_norm'*coef_resp)';
pc_resp = mat2cell(pc_c_resp, size(pc_c_resp,1), n_t(:)');
pc_resp = reshape(pc_resp, 3, 2);
% fix the direction of pc axis
for i = 1:size(coef_resp,2)
    d_pc = cellfun(@(x, t_w, t_p) x(i,t_w==round(t_p(3)))-x(i,t_w==round(t_p(2))), pc_resp, t_warp, t_point);
    if mean(d_pc, 'all')<0
        coef_resp(:,i) = -coef_resp(:,i);
    end
end
pc_resp = cellfun(@(x) (x'*coef_resp)', fr_norm, 'UniformOutput', false);

% plot unit location with projection weight
id = 3;

w_resp = normalize(abs(coef_resp), 1, 'range');
w_resp = w_resp .* sign(coef_resp);
for k = 1:n_shank
    x_shift = -mean(x_k{k})+k_sep*k;
    plot(ax_coef_1{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_1{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_resp(:,1))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_1{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_resp(id_this, 1)), 'k');

    plot(ax_coef_2{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_2{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_resp(:,2))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_2{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_resp(id_this, 2)), 'k');
end

% plot pc-1 and pc-2
text(ax_var{id}, 5, 100, sprintf('%.1f%%', explain_resp(1)), 'FontSize', 6);
text(ax_var{id}, 5, 85, sprintf('%.1f%%', explain_resp_tot(1)), 'FontSize', 6);
text(ax_var_2{id}, 5, 100, sprintf('%.1f%%', explain_resp(2)), 'FontSize', 6);
text(ax_var_2{id}, 5, 85, sprintf('%.1f%%', explain_resp_tot(2)), 'FontSize', 6);
for fp = 1:3
    for p = 1:2
        t_w = t_warp{fp, p};
        t_p = t_point{fp, p};
        pc_this = pc_resp{fp, p};

        ax = ax_hold{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_hold_2{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move_2{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);
    end
end
y_hold = ax_hold{id}.YLim;
y_move = ax_move{id}.YLim;

y_lim = [min([y_hold, y_move]), max([y_hold, y_move])];
set(ax_hold{id}, 'YLim', y_lim);
set(ax_move{id}, 'YLim', y_lim, 'YColor', 'none');

y_hold_2 = ax_hold_2{id}.YLim;
y_move_2 = ax_move_2{id}.YLim;

y_lim_2 = [min([y_hold_2, y_move_2]), max([y_hold_2, y_move_2])];
set(ax_hold_2{id}, 'YLim', y_lim_2);
set(ax_move_2{id}, 'YLim', y_lim_2, 'YColor', 'none');

text(ax_hold{id}, -5200, mean(ax_hold{id}.YLim), 'Response', 'Rotation', 90, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'FontSize', 8, 'HorizontalAlignment', 'center');

% choice
fr_norm_move = cellfun(@(x, t_w, t_p) x(:, t_w>=t_p(3) & t_w<t_p(4)), fr_norm, t_warp, t_point, 'UniformOutput', false);
fr_c_norm_move = cell2mat(fr_norm_move(:)'); % concatenate

[coef_move, ~, ~, ~, explain_move] = pca(fr_c_norm_move');
var_move_tot = diag(coef_move' * cov_tot * coef_move);
explain_move_tot = 100 * var_move_tot ./ var_tot;

pc_c_move = (fr_c_norm'*coef_move)';
pc_move = mat2cell(pc_c_move, size(pc_c_move,1), n_t(:)');
pc_move = reshape(pc_move, 3, 2);
% fix the direction of pc axis
for i = 1:size(coef_move,2)
    d_pc = cellfun(@(x, t_w, t_p) x(i,t_w==round(t_p(3)+200))-x(i,t_w==round(t_p(3))), pc_move, t_warp, t_point);
    if mean(d_pc, 'all')<0
        coef_move(:,i) = -coef_move(:,i);
    end
end
pc_move = cellfun(@(x) (x'*coef_move)', fr_norm, 'UniformOutput', false);

% plot unit location with projection weight
id = 4;

w_move = normalize(abs(coef_move), 1, 'range');
w_move = w_move .* sign(coef_move);
for k = 1:n_shank
    x_shift = -mean(x_k{k})+k_sep*k;
    plot(ax_coef_1{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_1{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_move(:,1))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_1{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_move(id_this, 1)), 'k');

    plot(ax_coef_2{id}, (x_k{k}(1)+x_shift)*[1 1], y_k{k}, '-k');
    plot(ax_coef_2{id}, (x_k{k}(2)+x_shift)*[1 1], y_k{k}, '-k');
    id_this = abs(w_move(:,2))>coef_thres & k_unit==shanks(k);
    scatter(ax_coef_2{id}, loc_unit(1,id_this)+x_shift, loc_unit(2,id_this), 24*abs(w_move(id_this, 2)), 'k');
end

% plot pc-1 and pc-2
text(ax_var{id}, 5, 100, sprintf('%.1f%%', explain_move(1)), 'FontSize', 6);
text(ax_var{id}, 5, 85, sprintf('%.1f%%', explain_move_tot(1)), 'FontSize', 6);
text(ax_var_2{id}, 5, 100, sprintf('%.1f%%', explain_move(2)), 'FontSize', 6);
text(ax_var_2{id}, 5, 85, sprintf('%.1f%%', explain_move_tot(2)), 'FontSize', 6);
for fp = 1:3
    for p = 1:2
        t_w = t_warp{fp, p};
        t_p = t_point{fp, p};
        pc_this = pc_move{fp, p};

        ax = ax_hold{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(1,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_hold_2{id};
        ind = find(t_w<=t_p(3));
        plot(ax, t_w(ind), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);

        ax = ax_move_2{id};
        ind = find(t_w>t_p(3));
        plot(ax, t_w(ind)-t_w(ind(1)), pc_this(2,ind), 'Color', c{fp, p}, 'LineWidth', 1.5);
    end
end
y_hold = ax_hold{id}.YLim;
y_move = ax_move{id}.YLim;

y_lim = [min([y_hold, y_move]), max([y_hold, y_move])];
set(ax_hold{id}, 'YLim', y_lim);
set(ax_move{id}, 'YLim', y_lim, 'YColor', 'none');

y_hold_2 = ax_hold_2{id}.YLim;
y_move_2 = ax_move_2{id}.YLim;

y_lim_2 = [min([y_hold_2, y_move_2]), max([y_hold_2, y_move_2])];
set(ax_hold_2{id}, 'YLim', y_lim_2);
set(ax_move_2{id}, 'YLim', y_lim_2, 'YColor', 'none');

text(ax_hold{id}, -5200, mean(ax_hold{id}.YLim), 'Choice', 'Rotation', 90, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'FontSize', 8, 'HorizontalAlignment', 'center');
xlabel(ax_hold{id}, 'Time to Cent-In / ');
xlabel(ax_move{id}, 'Cent-Out (ms)');


ax_coef_all = plt.assign_ax_to_fig(fig_pc, 1, 1, [13 10.5, 3, 3], [3 3]);
ax_coef_all = ax_coef_all{1};
set(ax_coef_all, 'YLim', [.5 n_u+.5], 'XColor', 'none', 'YColor', 'none', 'YDir', 'reverse', 'CLim', [-1 1]);

% cmap_w = customcolormap(linspace(0, 1, 7), {'#d75f4e','#b5172f','#68011d', '#000000','#062e61','#2265ad','#4295c1'});
cmap_w = customcolormap_preset('pink-white-green');
colormap(ax_coef_all, cmap_w);

imagesc(ax_coef_all, 1:2, 1:n_u, w_app(ind_sort{1},1:2));
imagesc(ax_coef_all, 4:5, 1:n_u, w_hold(ind_sort{1},1:2));
imagesc(ax_coef_all, 7:8, 1:n_u, w_resp(ind_sort{1},1:2));
imagesc(ax_coef_all, 10:11, 1:n_u, w_move(ind_sort{1},1:2));

cb = colorbar(ax_coef_all, 'Units', 'centimeters', 'Position', [16.5 10.5 .25 3], 'TickDirection', 'out', 'FontSize', 7);
cb.Label.String = "Norm. projection weight";
cb.Label.FontSize = 7;

set_fig_title(fig_pc, sprintf('%s | %s', r.BehaviorClass.Subject, r.BehaviorClass.Session));

fig_name = sprintf('PCA_%s_%s', r.BehaviorClass.Subject, r.BehaviorClass.Session);
fig_path = fullfile('Figure', fig_name);
exportgraphics(fig_pc, sprintf('%s.png', fig_path), 'Resolution', 300);
exportgraphics(fig_pc, sprintf('%s.pdf', fig_path), 'ContentType', 'vector');

%% Save PCA results
Coef.Approach = coef_app;
Coef.Hold     = coef_hold;
Coef.Response = coef_resp;
Coef.Choice   = coef_move;

PCs.Approach = pc_app;
PCs.Hold     = pc_hold;
PCs.Response = pc_resp;
PCs.Choice   = pc_move;

Explain.Approach = explain_app;
Explain.Hold     = explain_hold;
Explain.Response = explain_resp;
Explain.Choice   = explain_move;
ExplainTot.Approach = explain_app_tot;
ExplainTot.Hold     = explain_hold_tot;
ExplainTot.Response = explain_resp_tot;
ExplainTot.Choice   = explain_move_tot;

PC.Coef = Coef;
PC.PCs  = PCs;
PC.Explain = Explain;
PC.ExplainTot = ExplainTot;

pcs_name = sprintf('PCA_%s_%s.mat', r.BehaviorClass.Subject, r.BehaviorClass.Session);
save(pcs_name, 'PC');

end