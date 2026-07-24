function SpikeInfo = getSpikeInfo(r, Regions, manual)
%% Calculate all features (waveform, timing-related etc.) of units.
% Cluster units by trough-peak delay and waveform PC-1 to get regular units.
% Assign units to different regions according to recording depth.
%   Input:
%       r: RT array
%       Reigions: recording regions, each column for each shank; each row for each region, for example:
%           Regions = [
%               "M2" "M2"
%               "OFC" ""
%               ];
%       manual: 1*2 logical, whether regular-cluster and recording-region should be manually set (not yet)
%   Output:
%       SpikeInfo, a table with variables:
%           Subject, Session, NP, Unit, Channel, ChUnit, SpikeTimes, Waveform,
%           Xcoords, Ycoords, Kcoords, Shank, Location, Region,
%           Amplitude, AmpTrough, tTrough, AmpPeak, tPeak, PeakTroughRatio, PeakTroughDelay, fwhh, WavePC
%           AutoCorr, ISI, Single, Regular

if nargin<3
    manual = [0 0];
end

%% Basic information
Subject = r.BehaviorClass.Subject;
Session = r.BehaviorClass.Session;

SpikeTimes = struct2table(r.Units.SpikeTimes);
SpikeNotes = r.Units.SpikeNotes;

wave_ch_all  = SpikeTimes.wave_mean;
wave_ch_max  = SpikeTimes.wave;
spike_timing = SpikeTimes.timings;
n_units = height(SpikeTimes);

Fs = 30000; % Hz
n_sample = size(wave_ch_all{1}, 2);
t_wave   = (1:n_sample) * 1000 / Fs; % in ms

if all(r.ChanMap.connected)
    NP_ver = "2.0";
    wave_scale = 5 / 4;
else
    NP_ver = "1.0";
    wave_scale = 1 / 4;
end

% table to output
SpikeInfo = table(repmat(Subject, n_units, 1), repmat(Session, n_units, 1), repmat(NP_ver, n_units, 1), (1:n_units)', SpikeNotes(:,1), SpikeNotes(:,2), spike_timing, wave_ch_all, ...
    'VariableNames', {'Subject', 'Session', 'NP', 'Unit', 'Channel', 'ChUnit', 'SpikeTimes', 'Waveform'});


%% monopolar-triangulation location
ch_location = [r.ChanMap.xcoords r.ChanMap.ycoords];
unit_location = zeros(n_units, 3);
% unit_location = getUnitLocation(r);
for i = 1:n_units
    [loc_x, loc_y, loc_z] = spikeLocation(wave_ch_all{i}, ch_location, 20, 'monopolar_triangulation');
    unit_location(i, :) = [loc_x, loc_y, loc_z];
end
SpikeInfo.Xcoords  = repmat({r.ChanMap.xcoords}, n_units, 1);
SpikeInfo.Ycoords  = repmat({r.ChanMap.ycoords}, n_units, 1);
SpikeInfo.Kcoords  = repmat({r.ChanMap.kcoords}, n_units, 1);
SpikeInfo.Shank    = r.ChanMap.kcoords(SpikeInfo.Channel);
SpikeInfo.Location = unit_location;

shanks = unique(SpikeInfo.Shank);
num_shank = length(shanks);


%% waveform features (amp_trough, t_trough, amp_peak, t_peak, peak_trough_ratio, peak_trough_delay, fwhh, pc-1)
wave_mean = zeros(length(wave_ch_max), length(t_wave));
wave_features = struct();
for i = 1:n_units
    wave_all = wave_ch_max{i};

    [feature_out, wave_m_i] = get_waveform_features(t_wave, wave_all);
    wave_mean(i, :) = wave_m_i;
    if i==1
        wave_features = feature_out;
    else
        wave_features = [wave_features; feature_out];
    end
end
wave_features = struct2table(wave_features);

% correct the difference in gain between NP1.0 and NP2.0 (which has not been considered in BuildSpikeTable.m),
wave_features.amp_trough = wave_features.amp_trough * wave_scale;
wave_features.amp_peak   = wave_features.amp_peak   * wave_scale;

% PC-1
[~, wave_pc] = pca(zscore(wave_mean, [], 2));
wave_features.pc = wave_pc(:,1);

SpikeInfo.Amplitude = wave_features.amp_peak - wave_features.amp_trough;
SpikeInfo.AmpTrough = wave_features.amp_trough;
SpikeInfo.tTrough   = wave_features.t_trough;
SpikeInfo.AmpPeak   = wave_features.amp_peak;
SpikeInfo.tPeak     = wave_features.t_peak;
SpikeInfo.PeakTroughRatio = wave_features.peak_trough_ratio;
SpikeInfo.PeakTroughDelay = wave_features.peak_trough_delay;
SpikeInfo.fwhh   = wave_features.fwhh;
SpikeInfo.WavePC = wave_features.pc;


%% spike timing features (AutoCorr, ISI, ISI violation)
[AutoCorr, tAutoCorr] = cellfun(@(x) computeAutoCorr(x, 300, 1), SpikeInfo.SpikeTimes, 'UniformOutput', false);

[ISIViolation, ISI] = cellfun(@(x) isiViolations(x), SpikeInfo.SpikeTimes, 'UniformOutput', false);
ISIViolation = cell2mat(ISIViolation);

SpikeInfo.AutoCorr  = AutoCorr;
SpikeInfo.tAutoCorr = tAutoCorr;
SpikeInfo.ISI = ISI;
SpikeInfo.ISIViolation = ISIViolation;


%% plot waveform feature
% set figure
plt = GPSPlot();
fig_unit = figure(12); clf(fig_unit);
set_fig_default(fig_unit);
set(fig_unit, 'Position', [10 10 10 8]);

% plot waveform feature
ax_wave = plt.assign_ax_to_fig(fig_unit, 1, 1, [1 4 3 3], [3 3]);
ax_wave = ax_wave{1};

scatter(ax_wave, SpikeInfo.PeakTroughDelay, SpikeInfo.WavePC, 15, 'black', 'filled', 'MarkerFaceAlpha', .6);
xline(ax_wave, 0.5, ':k');
yline(ax_wave, 0, ':k');
set(ax_wave, 'XLim', [0 1]);
xlabel(ax_wave, 'Peak-trough delay (ms)');
ylabel(ax_wave, 'Waveform PC-1');

%% Cluster units by waveform features (peak_trough_delay vs. pc-1)
Single  = SpikeNotes(:,3)==1;
Regular = false(n_units, 1);
if ~manual(1)
    % regular units should have longer peak_trough_delay and normal PC
    cluster_X = [wave_features.peak_trough_delay, wave_features.pc];
    % at first, take units with longer peak-trough-dealy into consideration,
    ind_regular = cluster_X(:,1)>0.5;
%     % find units with mahalanobis-distance < 3
%     mal_dist = sqrt(mahal(cluster_X, cluster_X(ind_regular, :)));
    % find units within 3 standard deviance in both dimensions
    mal_dist_1 = sqrt(mahal(cluster_X(:,1), cluster_X(ind_regular,1)));
    mal_dist_2 = sqrt(mahal(cluster_X(:,2), cluster_X(ind_regular,2)));
    mal_dist = max([mal_dist_1 mal_dist_2], [], 2);
    ind_regular = mal_dist<3;
    % then, re-calculate mahalanobis-distance based the units found before,
%     % cluster units with mahalanobis-distance < 3 as "regular"
%     mal_dist = sqrt(mahal(cluster_X, cluster_X(ind_regular, :)));
%     ind_regular = mal_dist<3;
    % cluster units within 3 standard deviance in both dimensions as "regular"
    mal_dist_1 = sqrt(mahal(cluster_X(:,1), cluster_X(ind_regular,1)));
    mal_dist_2 = sqrt(mahal(cluster_X(:,2), cluster_X(ind_regular,2)));
    mal_dist = max([mal_dist_1 mal_dist_2], [], 2);
    ind_regular = mal_dist<3;

    Regular(ind_regular) = true;
end

SpikeInfo.Single  = Single;
SpikeInfo.Regular = Regular;


%% Update plot
c_r  = [0, 135, 100]/255;
c_ir = [210, 40, 35]/255;

% update waveform feature scatter
cla(ax_wave);
scatter(ax_wave, SpikeInfo.PeakTroughDelay(Regular & Single), SpikeInfo.WavePC(Regular & Single), 15, c_r, 'filled', 'MarkerFaceAlpha', .6);
scatter(ax_wave, SpikeInfo.PeakTroughDelay(Regular & ~Single), SpikeInfo.WavePC(Regular & ~Single), 10, c_r, 'LineWidth', .75, 'MarkerEdgeAlpha', .7);
scatter(ax_wave, SpikeInfo.PeakTroughDelay(~Regular), SpikeInfo.WavePC(~Regular), 10, c_ir, 'LineWidth', .75, 'MarkerEdgeAlpha', .7);
xline(ax_wave, 0.5, ':k');
yline(ax_wave, 0, ':k');

% show example waveform
n_show = 25;
ax_exam = plt.assign_ax_to_fig(fig_unit, 1, 2, [1 1 3 1.5], [1.35 1.35]);
cellfun(@(x) set(x, 'XLim', [0 max(t_wave)], 'YLim', [-500 300], 'XColor', 'none', 'YColor', 'none'), ax_exam);
% regular-single
regular_show = find(Regular & Single & SpikeInfo.AmpTrough>-500 & SpikeInfo.AmpPeak<200);
if length(regular_show)>n_show
    regular_show = regular_show(randperm(length(regular_show), n_show));
end
ax = ax_exam{1};
for i = 1:length(regular_show)
    id_i = regular_show(i);
    wave_i = SpikeInfo.Waveform{id_i}(SpikeInfo.Channel(id_i), :) * wave_scale;
    plot(ax, t_wave, wave_i, 'Color', c_r);
end
title(ax, 'regular-single');
% scalar
plot(ax, [.1 .6], [-450 -450], '-k', 'LineWidth', 1);
plot(ax, [.1 .1], [-450 -250], '-k', 'LineWidth', 1);
text(ax, .1, -450, '500 μs', 'FontSize', 6, 'VerticalAlignment', 'top');
text(ax, .1, -350, '200 μV ', 'FontSize', 6, 'HorizontalAlignment', 'right');

% irregular-narrow
irregular_show = find(~Regular & SpikeInfo.PeakTroughDelay<0.45 & SpikeInfo.AmpTrough>-500 & SpikeInfo.AmpPeak<200);
if length(irregular_show)>n_show
    irregular_show = irregular_show(randperm(length(irregular_show), n_show));
end
ax = ax_exam{2};
for i = 1:length(irregular_show)
    id_i = irregular_show(i);
    wave_i = SpikeInfo.Waveform{id_i}(SpikeInfo.Channel(id_i), :) * wave_scale;
    plot(ax, t_wave, wave_i, 'Color', c_ir);
end
title(ax, 'irregular');

%% plot unit location on probe
if num_shank==1
    ax_loc_width = 1.5;
else
    ax_loc_width = 3;
end
ax_loc = plt.assign_ax_to_fig(fig_unit, 1, 1, [5.5 1 ax_loc_width 6], [ax_loc_width 6]);
ax_loc = ax_loc{1};

dist_horizon = max(r.ChanMap.xcoords) - min(r.ChanMap.xcoords);
dist_horizon = max([.2*dist_horizon, 50]);
set(ax_loc, 'XLim', [min(r.ChanMap.xcoords)-dist_horizon max(r.ChanMap.xcoords)+dist_horizon], 'YLim', [min(r.ChanMap.ycoords)-50 max(r.ChanMap.ycoords)+50], ...
    'XTick', [min(r.ChanMap.xcoords) max(r.ChanMap.xcoords)], 'YTick', [min(r.ChanMap.ycoords) max(r.ChanMap.ycoords)]);
xlabel(ax_loc, 'Horizontal location (μm)');
ylabel(ax_loc, 'Distance from tip (μm)');

%% Assign units to different regions according to recording depth (LocY)
SpikeInfo.Region = strings(n_units, 1);
if size(Regions, 2)==1
    Regions = repmat(Regions, 1, num_shank);
end
for k = 1:num_shank % for each shank
    ind_k = find(SpikeInfo.Shank==shanks(k));
    scatter(ax_loc, SpikeInfo.Location(ind_k,1), SpikeInfo.Location(ind_k,2), 15, 'black', 'filled', 'LineWidth', 1, 'MarkerFaceAlpha', .6);
    scatter(ax_loc, ch_location(ind_k,1), ch_location(ind_k,2), 2, 'black', '.', 'LineWidth', 1, 'MarkerEdgeAlpha', .6);

    % assign regions
    regions_k = Regions(:,k);
    regions_k = regions_k(regions_k~="");
    n_regions = length(regions_k);
    if n_regions > 1
        depth_k = SpikeInfo.Location(ind_k, 2);
        if ~manual(2)
            % fit a gaussian-mixture model to cluster units to different region
            gm = fitgmdist(depth_k, n_regions, 'RegularizationValue', 0.01);
            [~, sort_idx] = sortrows(gm.mu, 1, "descend");
            gm_sorted = gmdistribution(gm.mu(sort_idx,:), gm.Sigma(:,:,sort_idx), gm.ComponentProportion(sort_idx));

            region_idx = cluster(gm_sorted, depth_k);

            for m = 1:n_regions
                region_m = find(region_idx==m);
                [~, ind_rm] = rmoutliers(depth_k(region_m), 'mean');
                region_m = region_m(~ind_rm);
                SpikeInfo.Region(ind_k(region_m)) = regions_k(m)';
            end
        else
            % separated by hand
            disp('Manually choose locations to separate regions');

            sep_k = zeros(n_regions-1, 1);
            for m = 1:n_regions-1
                [~, y_m]  = getpts(ax_loc);
                sep_k(m) = mean(y_m(1:end-1));
            end
            sep_k = sort(sep_k, 'descend');
            SpikeInfo.Region(ind_k) = regions_k(1);
            for m = 1:n_regions-1
                region_m = depth_k<sep_k(m);
                SpikeInfo.Region(ind_k(region_m)) = regions_k(m+1);
            end
        end
    else
        SpikeInfo.Region(ind_k) = regions_k;
    end
end

%% update unit location plot
cla(ax_loc);
ind_region = SpikeInfo.Region~="";
scatter(ax_loc, SpikeInfo.Location(Regular & ~Single & ind_region, 1), SpikeInfo.Location(Regular & ~Single & ind_region, 2), 10, c_r, 'LineWidth', .75, 'MarkerEdgeAlpha', .7);
scatter(ax_loc, SpikeInfo.Location(~Regular & ind_region, 1), SpikeInfo.Location(~Regular & ind_region, 2), 10, c_ir, 'LineWidth', .75, 'MarkerEdgeAlpha', .7);
scatter(ax_loc, SpikeInfo.Location(Regular & Single & ind_region, 1), SpikeInfo.Location(Regular & Single & ind_region, 2), 15, c_r, 'filled', 'MarkerFaceAlpha', .6);

scatter(ax_loc, ch_location(:,1), ch_location(:,2), 2, 'black', '.', 'LineWidth', 1, 'MarkerEdgeAlpha', .6);

% show regions
region_all = unique(Regions);
region_all = region_all(region_all~="");
n_regions  = length(region_all);
if n_regions > 1
    y_max = zeros(n_regions, 1);
    y_min = zeros(n_regions, 1);
    y_med = zeros(n_regions, 1);
    for i = 1:n_regions
        y_i = SpikeInfo.Location(SpikeInfo.Region==region_all(i), 2);
        y_max(i) = max(y_i);
        y_min(i) = min(y_i);
        y_med(i) = median(y_i);
    end
    [~, id_sort] = sort(y_med);
    region_all = region_all(id_sort);
    y_max = y_max(id_sort);
    y_min = y_min(id_sort);
    y_med = y_med(id_sort);

    y_sep = zeros(n_regions-1, 1);
    for i = 1:n_regions-1
        y_sep(i) = median([y_min(i+1) y_max(i)]);
        yline(ax_loc, y_sep(i), '--k', 'LineWidth', 1, 'Alpha', 1);
    end
else
    y_med = median(SpikeInfo.Location(:,2));
end
for i = 1:n_regions
    text(ax_loc, ax_loc.XLim(2), y_med(i), region_all(i), 'FontSize', 7, 'FontWeight', 'bold');
end

%% save figure
save_name = sprintf('SpikeInfo_%s_%s', Subject, Session);
exportgraphics(fig_unit, save_name+".png", 'Resolution', 600);

save(save_name+".mat", 'SpikeInfo');

end % getSpikeInfo
