function WarpOut = PSTHLiteWarpedSimple(r, id)
% 3.9.2024 plot a simple version of PSTH
% 4/23/2024 Plot warped PSTH
% 5/15/2025 

%%
% figure out which unit it is
spk_note = r.Units.SpikeNotes;
ind_unit = find(spk_note(:, 1)==id(1) & spk_note(:, 2)==id(2));
spkwave  = r.Units.SpikeTimes(ind_unit).wave;
spktimes = r.Units.SpikeTimes(ind_unit).timings; % in ms
spktimes_ = zeros(1, ceil(max(spktimes)));
spktimes_(round(spktimes)) = 1;
[ar, lags]  = xcorr(spktimes_, 25);
ar(lags==0) = 0;

% get base information
FPs   = r.PopPSTH.FPs;
Ports = r.PopPSTH.Ports;

nFPs   = length(FPs);
nPorts = length(Ports);

% Set colors
c = GPSColor();
p_c = .9;
c_port = {p_c*c.Contra+(1-p_c)*[1 1 1], p_c*c.Ipsi+(1-p_c)*[1 1 1]};

% For PSTH and raster plots
c_cent_in  = [5 191 219] / 255;
c_trigger  = [247 182 45] / 255;
c_cent_out = [238 5 219] / 255;
c_reward   = [164 208 164] / 255;
c_precor   = [0 0 0];
c_preerr   = [160 82 45] / 255;

press_col   = [5 5 5]/255;
trigger_col = [242 182 250]/255;
release_col = [87, 108, 188]/255;
poke_col    = [164, 208, 164]/255;

if nFPs == 2
    FP_cols = [192, 127, 0; 76, 61, 61]/255;
else
    FP_cols = [255, 217, 90; 192, 127, 0; 76, 61, 61]/255;
end

%%
%% Plot press
sdf_centin_mean = cell(nFPs, nPorts);
sdf_centin_ci = cell(nFPs, nPorts);

% the whole duration of time warp
pre_  = 2.5; % 2.5 sec before cent-in, about the lower boundary of shuttle time
post_ = 2.5; % 2.5 sec after choice
post_keep = 2;
% 
% hf=73;
% figure(hf); clf
% set(gcf, 'unit', 'centimeters', 'position', [2 2 17 16], 'paperpositionmode', 'auto',...
%     'renderer','Painters', 'Visible', 'on', 'color', 'w')
% x_now = 2;
% x_now_org = x_now;
% y_now_org = 13;
% x_size = 5;
% y_size = 2;
% y_size2 = 2;
% nplot = 10;

% time range aligned to different events
CentInRange  = [-2500 2500];
CentOutRange = [-1000 3000];
ChoiceRange  = [-1500 1000];
WarpedRange  = [-2000 5000];
% 
% hspacing = .6;
% x_now2 = x_now + x_size + hspacing;
% x_size_centout = x_size*(diff(CentOutRange)) / diff(CentInRange);
% x_now3 = x_now2 +x_size_centout + hspacing;
% 
% x_size_choice = x_size*(diff(ChoiceRange)) / diff(CentInRange);
% x_size_warped = x_size*(diff(WarpedRange)) / diff(CentInRange);
% 
% yrange = [0 40];
% vspacing = .5;
% vspacing_raster = 0.25;
% rng(0)
% FRmax = 20;
% y_now = y_now_org;

warped_time_points = cell(nFPs, nPorts);

spikes_trials = cell(1000, 4); % contains most data, 1:[cent-in cent-out choice]; 2:spike times in this duration; 3:iFP; 4:iPort
n_count = 0;
latency = 2.5; % max cent-out to choice latency

for ind_FP = 1:nFPs
    iFP = FPs(ind_FP);
    for ind_Port = 1:nPorts
        iPort = Ports(ind_Port);

        ind_this = find(r.PSTH.Events.CentIn.FP{1}==iFP & r.PSTH.Events.CentIn.Port{1}==ind_Port);

        centin_times  = sort(r.PSTH.Events.CentIn.Time{ind_this});
        centout_times = sort(r.PSTH.Events.CentOut.Time{ind_this});
        choice_times  = sort(r.PSTH.Events.Pokes.RewardPoke.Time{ind_FP, ind_Port});

        % cent-in times and cent-out times should be matched in num. this is what
        % we stretch
        hold_pairs = [(centin_times), (centout_times)];
        ioc_seq = []; % cent(I)n - cent(O)ut - (C)hoice sequence (each row represents a trial)
        ioc_spktimes = {};
        ioc_sdfs = {};
    
        % check if a choice-poke follows cent-out time, if not, we won't include it
        for j = 1:size(hold_pairs, 1)
            j_centout = hold_pairs(j, 2);
            % where is the choice poke
            j_choice = choice_times(find(choice_times>j_centout, 1, 'first'));

            if ~isempty(j_choice) && isempty(find(centin_times>j_centout & centin_times<j_choice, 1)) && (j_choice-j_centout<latency*1000)
                % if there is a choice poke after cent-out && there is no
                % cent-in between cent-out and choice poke && there is a
                % choice poke within a certain latency after cent-out
                ioc_seq   = [ioc_seq; hold_pairs(j, :) j_choice];
                total_dur = round(j_choice-hold_pairs(j, 1)) + pre_*1000 + post_*1000; % the unit is ms, [pre, from cent-in to choice, post]
                this_spk_train = spktimes(spktimes>=hold_pairs(j, 1)-pre_*1000 & spktimes<=j_choice+post_*1000);
                ioc_spktimes = [ioc_spktimes {this_spk_train}];

                % save this sequence
                n_count = n_count + 1;
                spikes_trials(n_count, :) = {[hold_pairs(j,:) j_choice], this_spk_train, iFP, iPort};
    
                % map spike times to spike matrix, then calculate spike density function
                tspk = (0:total_dur-1) - pre_*1000; % in ms
                spkmat = zeros(1, total_dur); % each entry represents for 1 ms
                if ~isempty(this_spk_train)
                    % align to cent-in
                    this_spk_train = this_spk_train - (hold_pairs(j, 1));
                    % convert the spike train to sdf
                    [~, ind_spikes] = intersect(round(tspk), round(this_spk_train));
                    spkmat(ind_spikes) = 1;
                end
                spkout   = sdf(tspk/1000, spkmat, 20); % spkout = sdf(tspk, spkin, kernel_width)
                ioc_sdfs = [ioc_sdfs [tspk; spkout']]; % ioc_sdfs is a cell
            end
        end
    
        % let's warp this thing
        median_hold_duration = median(ioc_seq(:,2)-ioc_seq(:,1)); % median, in ms
        median_movement_time = median(ioc_seq(:,3)-ioc_seq(:,2)); % median, in ms
    
        % jt_template has the following structures: end of FP, median hold
        % duration, median movement time
        jt_template = [0 iFP*1000 median_hold_duration median_movement_time+median_hold_duration]; % [cent-in, trigger, cent-out, choice] in ms
        warped_time_points{ind_FP, ind_Port} = jt_template;
        dt = 1; % 1 ms
        jt_target_time = (0:dt:median_movement_time+median_hold_duration);
        sprintf('Critical time points are %2.2f\n', jt_template(1), jt_template(2), jt_template(3), jt_template(4));
        t_warptarget_first  = jt_target_time(jt_target_time>=iFP*1000 & jt_target_time<jt_template(3)); % from trigger to cent-out
        t_warptarget_second = jt_target_time(jt_target_time>=jt_template(3) & jt_target_time<jt_template(4)); % from cent-out to choice

        ioc_spk_warped{ind_FP, ind_Port} = [];
        for j = 1:size(ioc_seq) % for each trial
            % jt is [cent-in, cent-out, choice]
            jt = ioc_seq(j, :) - ioc_seq(j, 1); % normalize so that the first point is time 0
            jsdf = ioc_sdfs{j};
            tsdf = jsdf(1, :); % this is the time of sdf, defined previously as: tspk = (0:total_dur-1)-pre_*1000; % in ms
            jsdf = jsdf(2, :); % this is the sdf
            not_warped = jsdf(tsdf<iFP*1000);
            %        sprintf('pre-press duration is %2.2f ms', length(not_warped))
            towarp_first   = jsdf(tsdf>=iFP*1000 & tsdf<jt(2)); % from trigger to cent-out
            t_towarp_first = tsdf(tsdf>=iFP*1000 & tsdf<jt(2)); % from trigger to cent-out
            towarp_second   = jsdf(tsdf>=jt(2) & tsdf<jt(3)); % from cent-out to choice
            t_towarp_second = tsdf(tsdf>=jt(2) & tsdf<jt(3)); % from cent-out to choice
            not_warped2 = jsdf(tsdf>jt(3));
            %        sprintf('post-poke duration is %2.2f ms', length(not_warped2))
            not_warped2 = not_warped2(1:post_keep*1000); % max 2 sec after choice poke
    
            dt = 1; % 1 ms
            if ~isempty(t_towarp_first) && ~isempty(t_towarp_second)
                sdf_warped_first  = SpikesGPS.SRT.warp_sdf(t_towarp_first, towarp_first, t_warptarget_first); % input is V, X, and duration to warp
                sdf_warped_second = SpikesGPS.SRT.warp_sdf(t_towarp_second, towarp_second, t_warptarget_second);
                new_sdf = [not_warped, sdf_warped_first, sdf_warped_second, not_warped2];
                ioc_spk_warped{ind_FP, ind_Port} = [ioc_spk_warped{ind_FP, ind_Port}; new_sdf];
            end
        end
    
        t_warped{ind_FP, ind_Port} = (-pre_*1000:dt:length(new_sdf)-pre_*1000-dt);
%     
%         % Note that centin_times, spkmat, etc. have been sorted by reaction time
%         spkmat         = r.PSTH.PSTHs(ind_unit).CentIn{ind_FP, ind_Port}{3};  
%         t_spkmat       = r.PSTH.PSTHs(ind_unit).CentIn{ind_FP, ind_Port}{4};
%         reaction_times = r.PSTH.PSTHs(ind_unit).CentIn{ind_FP, ind_Port}{6};
%         centin_times   = r.PSTH.Events.CentIn.Time{ind_this};
%         centout_times  = r.PSTH.Events.CentOut.Time{ind_this};
%     
%         centin_raster{ind_FP, ind_Port} = {t_spkmat, spkmat};
%     
%         % plot raster
%         nplot = size(spkmat, 2); % number of trials
%         ha_CentIn_cell(ind_FP, ind_Port) = axes('unit', 'centimeters', ...
%             'position', [x_now y_now x_size y_size2],...
%             'nextplot', 'add', ...
%             'xlim', CentInRange, 'ylim', [0 nplot+1], 'ydir', 'reverse', ...
%             'ytick', [0 20], 'xtick', -2000:1000:2000,...
%             'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%             'XTickLabelRotation', 0, 'color', 'none', ...
%             'ticklength',[.015 .1]);
%         axis off
%         line([0 0], [0 nplot+1], 'color', 'k', 'linewidth', 1);
%         line([iFP iFP]*1000, [0 nplot+1], 'color', c_trigger, 'linewidth', 1);
%     
%         y_now = y_now - y_size2 - vspacing_raster;
%         ind_plot = (1:nplot);
%     
%         spkmat_plot  = spkmat(:, ind_plot);
%         centin_times = centin_times(ind_plot);
%         reaction_times = reaction_times(ind_plot);
%     
%         for k = 1:nplot
%             xx = t_spkmat((spkmat_plot(:, k)==1));
%             yy = [0 .8] + k;
%             if ~isempty(xx)
%                 line([xx; xx], yy, 'color', c_port{ind_Port});
%             end
%             % Plot reaction time
%             xx = centout_times(k) - centin_times(k);
%             yy = [0 1] + k;
%             if ~isempty(xx)
%                 line([xx; xx], yy, 'color', c_cent_out, 'linewidth', 1)
%             end
%             % find poke time
%             k_poketime = find(choice_times>centin_times(k), 1, 'first');
%             xx = choice_times(k_poketime) - centin_times(k);
%             yy = .5 + k;
%             if ~isempty(xx)
%                 plot(xx, yy, 'o', 'color', c_reward, 'linewidth', 1, 'markerfacecolor', c_reward, 'markersize', 4, 'markeredgecolor', 'w')
%             end
%         end
%     
%         spkout   = sdf(t_spkmat/1000, spkmat, 20);
%         sdf_mean = mean(spkout, 2);
%         sdf_ci   = transpose(bootci(1000, @mean, spkout'));
%         sdf_centin_mean{ind_FP, ind_Port} = sdf_mean;
%         sdf_centin_ci{ind_FP, ind_Port} = sdf_ci;
    end
end
% 
% y_now2 = y_now-vspacing;
% 
% %% plot PSTH
% ha_Press_PSTH =  axes('unit', 'centimeters', ...
%     'position', [x_now  y_now2 x_size y_size],...
%     'nextplot', 'add', ...
%     'xlim', [CentInRange], 'ylim', yrange, ...
%     'ytick', [0:20:100], 'xtick', [-2000:1000:2000],...
%     'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%     'XTickLabelRotation', 0, 'color', 'none', ...
%     'ticklength',[.025 .1]);
% line_press = line([0 0], yrange, 'color', 'k', 'linewidth', 1);
% press_psth = {t_spkmat, sdf_centin_mean, sdf_centin_ci};
% psth_info_explained = 'time, sdf_mean, sdf_ci';
% 
% for i =1:length(FPs)
%     plotshaded(t_spkmat, sdf_centin_ci{i}, [.6 .6 .6])
%     plot(t_spkmat, sdf_centin_mean{i}, 'linewidth', 2, 'color', FP_cols(i,:))
%     if max(sdf_centin_mean{i})>FRmax
%         FRmax =max(sdf_centin_mean{i});
%     end
%     if i == 2
%         xlabel('Time from press (ms)')
%         ylabel('Spike rate (Hz)')
%     end
% end
% 
% %% Plot release
% y_now =y_now_org;
% sdf_release_mean = cell(1,nFPs);
% sdf_release_ci = cell(1, nFPs);
% 
% % all presses
% tPresses = sort(r.PSTH.Events.Presses.Time{end});
% 
% for ind_FP = 1:nFPs
% 
%     tRelease_indFP = r.PSTH.PSTHs(ind_unit).Releases{ind_FP}{end-1};
%     spkmat = r.PSTH.PSTHs(ind_unit).Releases{ind_FP}{3};
%     nplot = size(spkmat, 2);
%     ind_plot = (1:nplot);
%     t_spkmat = r.PSTH.PSTHs(ind_unit).Releases{ind_FP}{4};
%     reaction_times = r.PSTH.PSTHs(ind_unit).Releases{ind_FP}{6};
%     release_raster{ind_FP} = {t_spkmat, spkmat};
%     % plot raster
%     ha_Release_cell =  axes('unit', 'centimeters', ...
%         'position', [x_now2  y_now x_size_centout y_size2],...
%         'nextplot', 'add', ...
%         'xlim', [CentOutRange], 'ylim', [0 nplot+1], 'ydir', 'reverse', ...
%         'ytick', [0 20], 'xtick', [-2000:1000:2000],...
%         'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%         'XTickLabelRotation', 0,  'color', 'none', ...
%         'ticklength',[.025 .1]);
%     axis off
%     y_now = y_now-y_size2 -vspacing_raster;
%     line([0 0], [0 nplot+1], 'color', 'k', 'linewidth', 1)
% 
%     spkmat_plot = spkmat(:, ind_plot);
%     reaction_times = reaction_times(ind_plot);
% 
%     for k =1:nplot
%         xx = t_spkmat((spkmat_plot(:, k)==1));
%         yy = [0 .8]+k;
%         if ~isempty(xx)
%             line([xx; xx], yy, 'color', FP_cols(ind_FP, :))
%         end
% 
%         tk_Release = tRelease_indFP(k);
%         tk_Press = tPresses(find(tPresses<tk_Release, 1, 'last'));
% 
%         % Plot reaction time
%         xx = -reaction_times(k);
%         xx = -(tk_Release-tk_Press-FPs(ind_FP));
%         yy = [0 1]+k;
%         if ~isempty(xx)
%             line([xx; xx], yy, 'color', trigger_col, 'linewidth', 1)
%         end
%     end
% 
%     spkout=sdf(t_spkmat/1000, spkmat, 20);
%     sdf_mean = mean(spkout, 2);
%     sdf_ci   = transpose(bootci(1000, @mean, spkout'));
%     sdf_release_mean{ind_FP} = sdf_mean;
%     sdf_release_ci{ind_FP} = sdf_ci;
% 
% end
% y_now2 = y_now-vspacing;
% 
% % plot PSTH
% ha_Release_PSTH =  axes('unit', 'centimeters', ...
%     'position', [x_now2  y_now2 x_size_centout y_size],...
%     'nextplot', 'add', ...
%     'xlim', [CentOutRange], 'ylim', yrange, ...
%     'ytick', [0:20:100], 'xtick', [-2000:1000:2000],...
%     'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%     'yticklabel', [],...
%     'XTickLabelRotation', 0, 'color', 'none', ...
%     'ticklength',[.025 .1]);
% line_release = line([0 0], yrange, 'color', 'k', 'linewidth', 1);
% 
% if i == 2
%     xlabel('Release (ms)')
% end
% 
% for i =1:length(FPs)
%     plotshaded(t_spkmat, sdf_release_ci{i}, [.6 .6 .6])
%     plot(t_spkmat, sdf_release_mean{i}, 'linewidth', 2, 'color', FP_cols(i, :))
%     if i == 2
%         xlabel('Release (ms)')
%     end
%     if max(sdf_release_mean{i})>FRmax
%         FRmax =max(sdf_release_mean{i});
%     end
% end
% release_psth = {t_spkmat, sdf_release_mean, sdf_release_ci};
% 
% %% Plot poke
% y_now =y_now_org;
% sdf_poke_mean = cell(1,nFPs);
% sdf_poke_ci = cell(1, nFPs);
% 
% for ind_FP = 1:nFPs
%     spkmat = r.PSTH.PSTHs(ind_unit).RewardPokes{ind_FP}{3};
%     nplot = size(spkmat, 2);
%     ind_plot = (1:nplot);
%     t_spkmat = r.PSTH.PSTHs(ind_unit).RewardPokes{ind_FP}{4};
%     move_times = r.PSTH.PSTHs(ind_unit).RewardPokes{ind_FP}{6};
%     poke_raster{ind_FP} = {t_spkmat, spkmat};
% 
%     % plot raster
%     ha_Poke_cell =  axes('unit', 'centimeters', ...
%         'position', [x_now3  y_now x_size_choice y_size2],...
%         'nextplot', 'add', ...
%         'xlim', [ChoiceRange], 'ylim', [0 nplot+1], 'ydir', 'reverse',...
%         'ytick', [0 20], 'xtick', [-2000:1000:2000],...
%         'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%         'XTickLabelRotation', 0,  'color', 'none', ...
%         'ticklength',[.025 .1]);
%     axis off
%     y_now = y_now-y_size2 -vspacing_raster;
%     line([0 0], [0 nplot+1], 'color', 'k', 'linewidth', 1)
%     spkmat_plot = spkmat(:, ind_plot);
%     move_times = move_times(ind_plot);
% 
%     for k =1:nplot
%         xx = t_spkmat((spkmat_plot(:, k)==1));
%         yy = [0 .8]+k;
%         if ~isempty(xx)
%             line([xx; xx], yy, 'color', FP_cols(ind_FP, :))
%         end
%         % Plot reaction time
%         xx = -move_times(k);
%         yy = [0 1]+k;
%         if ~isempty(xx)
%             line([xx; xx], yy, 'color', release_col, 'linewidth', 1)
%         end
%     end
% 
%     spkout=sdf(t_spkmat/1000, spkmat, 20);
%     sdf_mean = mean(spkout, 2);
%     sdf_ci   = transpose(bootci(1000, @mean, spkout'));
%     sdf_poke_mean{ind_FP} = sdf_mean;
%     sdf_poke_ci{ind_FP} = sdf_ci;
% 
% end
% y_now2 = y_now-vspacing;
% 
% % plot PSTH
% ha_Poke_PSTH =  axes('unit', 'centimeters', ...
%     'position', [x_now3  y_now2 x_size_choice y_size],...
%     'nextplot', 'add', ...
%     'xlim', [ChoiceRange], 'ylim', yrange, ...
%     'ytick', [0:20:100], 'xtick', [-2000:1000:2000],...
%     'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%     'yticklabel', [],...
%     'XTickLabelRotation', 0, 'color', 'none', ...
%     'ticklength',[.025 .1]);
% line_poke = line([0 0], yrange, 'color', 'k', 'linewidth', 1);
% 
% for i =1:length(FPs)
%     plotshaded(t_spkmat, sdf_poke_ci{i}, [.6 .6 .6])
%     plot(t_spkmat, sdf_poke_mean{i}, 'linewidth', 2, 'color', FP_cols(i,:))
%     if i == 2
%         xlabel('Poke (ms)')
%     end
%     if max(sdf_poke_mean{i})>FRmax
%         FRmax =max(sdf_poke_mean{i});
%     end
% end
% poke_psth = {t_spkmat, sdf_poke_mean, sdf_poke_ci};
% 
% % update a few plots
% set(ha_Press_PSTH, 'ylim', [0 FRmax]*1.2);
% set(ha_Release_PSTH, 'ylim', [0 FRmax]*1.2);
% set(ha_Poke_PSTH, 'ylim', [0 FRmax]*1.2);
% set(line_press, 'ydata', [0 FRmax]*1.2);
% set(line_release, 'ydata', [0 FRmax]*1.2);
% set(line_poke, 'ydata', [0 FRmax]*1.2);
% 
% %% Plot spike waveforms
% y_now = 1.5;
% x_now = x_now3 + x_size_choice+1.5;
% height0 = 1.5;

spks = spkwave;
tspk = (1:size(spks, 2))/30; % in ms
spk_mean = mean(spks, 1);
spk_std = std(spks, 0, 1);
% yrange = [min(spk_mean)*1.2 max(spk_mean)*1.2];
% 
% axes('unit', 'centimeters', ...
%     'position', [x_now  y_now height0  height0],...
%     'nextplot', 'add', ...
%     'xlim', [0 tspk(end)], 'ylim', yrange,...
%     'ytick', [0 20], 'xtick', [0:5],...
%     'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%     'XTickLabelRotation', 40, 'yticklabel',[], 'xticklabel', [], 'color', 'none');
% line([1 2], [yrange(1) yrange(1)], 'color', 'k')
% 
% plotshaded(tspk, [spk_mean-spk_std; spk_mean+spk_std], [0.5 0.5 0.5]);
% plot(tspk, spk_mean, 'k', 'linewidth', 1)
% 
% axis off
% 
% %% Plot auto-correlation
% y_now = 1.5+height0+2;
% 
% axes('unit', 'centimeters', ...
%     'position', [x_now  y_now 2  height0],...
%     'nextplot', 'add', ...
%     'xlim', [-20 20], 'ylim', [0 100],...
%     'xtick', [-20:10:20],...
%     'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%     'XTickLabelRotation', 40, 'color', 'none', 'ticklength', [0.05 0.1]);
% bar(lags, ar)
% xlabel('Lag (ms)')
% ylabel('Frequency')
% axis 'auto y'
% 
% %% Add warped PSTH
% x_now = x_now_org;
% y_now = 1.5;
% x_size_warped = 10;
% ha_Warped_PSTH =  axes('unit', 'centimeters', ...
%     'position', [x_now  y_now x_size_warped y_size*1.5],...
%     'nextplot', 'add', ...
%     'xlim', [WarpedRange], 'ylim', [0 FRmax]*1.2, ...
%     'ytick', [0:20:100], 'xtick', [-2000:1000:2000],...
%     'xscale', 'linear', 'yscale', 'linear', 'ticklength', [0.02, 1], ...
%     'XTickLabelRotation', 0, 'color', 'none', ...
%     'ticklength',[.025 .1]);
% line_press = line([0 0], yrange, 'color', 'k', 'linewidth', 1);
% warped_time_points{ind_FP}

sdf_warped_ci = cell(nFPs, nPorts);
sdf_warped_mean = cell(nFPs, nPorts);
for i = 1:nFPs
    for j = 1:nPorts
        if size(ioc_spk_warped{i}, 1)>10
            sdf_warped_ci{i,j} = bootci(1000, @mean, ioc_spk_warped{i,j});
            % plotshaded(t_warped{i}, sdf_warped_ci{i}, [.6 .6 .6])
        end

        sdf_warped_mean{i,j} = mean(ioc_spk_warped{i,j}, 1);
        % plot(t_warped{i}, sdf_warped_mean{i}, 'linewidth', 2, 'color', FP_cols(i, :))
        % line([FPs(i) FPs(i)], [0 FRmax]*1.2, 'color', trigger_col, 'linestyle', '--', 'linewidth', 1)
    end
end
% xlabel('Time from press (ms)')
% ylabel('Spike rate (Hz)')

%%
% 
% %% added 4/28/2025
% % Trigger-Release-Poke sequence.
% % 500 ms before trigger, then trigger, Trigger-Release (RT, need warping to
% % median RT), Release to Poke (also needs to be warped)
% % code initially written in cal_spk_out.m
% 
% % warp psth in this way: approach to press: no change; press to
% % press+750 ms, no change. trigger-250 ms to release: warp; release to
% % poke: warp;
% % function warp_out=warp_srt_sdfs_pooled(prp_sequence, prp_sdfs, toplot)
% 
% all_pokes                   =  sort(cell2mat(r.PSTH.Events.Pokes.RewardPoke.Time));
% pre_    =   2.75; % before press
% post_   =   1; % post press
% pre__   =   0.75; % pre trigger
% post__  =   1; % post poke
% 
% press_spktrains_all = cell(1, length(FPs));
% press_sdfs_all      = cell(1, length(FPs));
% release_sdfs_all      = cell(1, length(FPs));
% press_release_seqs_all = cell(1, length(FPs));
% press_release_seqs_org = cell(1, length(FPs));
% 
% release_event_seqs_all = [];
% release_spktrains_all = cell(1, length(FPs));
% 
% for ind_FP =1:length(FPs)
%     % find a few good sequence
%     centin_times              =  sort(r.PSTH.Events.Presses.Time{ind_FP});
%     trigger_times            =  centin_times + FPs(ind_FP);
%     centout_times            =  sort(r.PSTH.Events.Releases.Time{ind_FP});
%     trigger_times_   = [];
%     release_times_    = [];
%     poke_times_     = [];
%     for m =1:length(centout_times)
%         if ~isempty(find(all_pokes>centout_times(m), 1, 'first'))
%             poke_times_     = [poke_times_; all_pokes((find(all_pokes>centout_times(m), 1, 'first')))];
%             release_times_  = [release_times_; centout_times(m)];
%             trigger_times_  = [trigger_times_; trigger_times(m)];
%         end
%     end
% 
%     reaction_times_        =      release_times_  -   trigger_times_;
%     retrieval_durs           =      poke_times_     -   release_times_;
%     ind_included            =      ~isoutlier(retrieval_durs, 'ThresholdFactor', 10) & reaction_times_>100; % make sure only responses with RT over 100 ms are counted
% 
%     release_times_=release_times_(ind_included);
%     trigger_times_=trigger_times_(ind_included);
%     poke_times_=poke_times_(ind_included);
%     % also get rid of press times that are not coupled to a sequence
%     press_times_ = zeros(size(trigger_times_));
%     for m =1:size(trigger_times_)
%         press_times_(m) = centin_times(find(centin_times<trigger_times_(m), 1, 'last'));
%     end
%     % here, number of presses is likely larger than the number of
%     % trigger-release-poke sequence
%     press_release_seqs = [press_times_ trigger_times_ release_times_ poke_times_]; % warping required.
%     spikes_trials = cell(1000, 2); % contains most data, a holder
% 
%     % Go through these events and make sdf out of it.
%     sigma_kernel = 25; % gaussian kernel to make spike density function
%     dt = 1; % time bins, 1 ms
%     press_sdfs =  {};
%     press_spktrains = {};
%     press_release_seqs_= [];
%     % Go through each trigger-release-poke sequence
%     release_spktrains = {};
%     release_sdfs = {};
% 
%     k_ = 0;
%     for k =1:size(press_release_seqs, 1)
%         k_press     =       press_release_seqs(k, 1);
%         k_trigger   =       press_release_seqs(k, 2);
%         k_release   =       press_release_seqs(k, 3);
%         k_poke      =       press_release_seqs(k, 4);
%         if k_press-pre_*1000>0 && k_press+post_*1000<max(spktimes) && k_trigger-pre__*1000>0 && k_poke+post__*1000<max(spktimes)
%             k_                                      =       k_+1;
%             press_release_seqs_(k_, :)   =       press_release_seqs(k, :);
%             t_range                             =       [k_press-pre_*1000 k_press+post_*1000];
%             k_spktimes                       =       spktimes(spktimes>=t_range(1) & spktimes<=t_range(2));
%             k_spktimes                       =       k_spktimes-k_press; %  time aligned to press time
%             t_range                             =       t_range-k_press;
%             press_spktrains{k_}          =       k_spktimes; % spk train saved to a cell array.
%             [spkout, tspk]                  =       sdf25(k_spktimes, t_range, sigma_kernel, dt);  %  spkout=sdf(tspk, spkin, kernel_width
%             press_sdfs{k_}                  =       [tspk' spkout'];
% 
%             t_range                         =   [k_trigger-pre__*1000 k_poke+post__*1000];
%             k_spktimes                  =   spktimes(spktimes>=t_range(1) & spktimes<=t_range(2));
%             k_spktimes                  =   k_spktimes-k_trigger; %  time aligned to trigger time
%             release_spktrains{k_}   =   k_spktimes; % spk train saved to a cell array.
%             % document the relative timings
%             rel_timing                   = [0 k_release-k_trigger k_poke-k_release];
%             release_event_seqs_all = [release_event_seqs_all; rel_timing];
%             [spkout, tspk]              =   sdf25(k_spktimes, t_range-k_trigger, sigma_kernel, dt);  %  spkout=sdf(tspk, spkin, kernel_width)
%             release_sdfs                 =   [release_sdfs, [tspk' spkout']];
%         end
%     end
% 
%     press_spktrains_all{ind_FP}                    =          press_spktrains;
%     press_sdfs_all{ind_FP}                          =            press_sdfs;
%     press_release_seqs_all{ind_FP}              =           press_release_seqs_;
%     release_spktrains_all{ind_FP}                =           release_spktrains;
%     release_sdfs_all{ind_FP}                        =           release_sdfs;
%     
% end
% 
% warp_out.event_sequence_label         =          {'press', 'trigger', 'release', 'poke'};
% warp_out.event_sequence                  =           press_release_seqs_all;
% press_sdf_FPs                                      =           cell(1, length(press_sdfs_all));
% press_sdf_FPs_mean_ci                       =           cell(1, length(press_sdfs_all));
% press_sdf_pooled = [];
% 
% for ii =1:length(press_sdfs_all)
%     for jj = 1:length(press_sdfs_all{ii})
%         press_sdf_pooled = [press_sdf_pooled press_sdfs_all{ii}{jj}(:, 2)];
%         press_sdf_FPs{ii} = [press_sdf_FPs{ii} press_sdfs_all{ii}{jj}(:, 2)];
%      end
%     press_sdf_FPs_mean_ci{ii} = [press_sdfs_all{ii}{1}(:, 1) mean(press_sdf_FPs{ii}, 2) transpose(bootci(1000, @mean, press_sdf_FPs{ii}'))];
%  end
% 
% warp_out.press.spk_train_FPs_explained          =      'These are press-related spike trains and sdfs, each cell is a FP';
% warp_out.press.spk_train_FPs                           =       press_spktrains_all;
% warp_out.press.time                                         =       press_sdfs_all{1}{1}(:, 1);
% warp_out.press.sdf_FPs                                    =       press_sdf_FPs;
% warp_out.press.sdf_FPs_mean_ci                      =       press_sdf_FPs_mean_ci;  
% 
% warp_out.press.sdf_pooled.time                      =       press_sdfs_all{1}{1}(:, 1);
% warp_out.press.sdf_pooled.sdf_mean                  =       mean(press_sdf_pooled, 2);
% warp_out.press.sdf_pooled.sdf_ci                    =       transpose(bootci(1000, @mean, press_sdf_pooled'));
% warp_out.press.sdf_pooled.sdf_trials               =       press_sdf_pooled;
% 
% warp_out.release.spk_train_FPs_explained          =       'These release events are actually trigger-related spike trains and sdfs, each cell is a FP';
% warp_out.release.spk_train_FPs                    =        release_spktrains_all; % note that we actually align the spikes with trigger, not release. 
% warp_out.release.sdf_FPs                          =        release_sdfs_all; 
%  
% % merge different FPs
% release_spktrains_all = horzcat(release_spktrains_all{:});
% release_sdfs_all = horzcat(release_sdfs_all{:});
% % Calculate the median value of trigger-to-release latency and
% % release-to-poke latency
% rt_median                               =       round(median(release_event_seqs_all(:, 2)));
% retrieval_median                    =       round(median(release_event_seqs_all(:, 3)));
% target_time_trigger                 =       (-pre__*1000:-1);
% target_time_rt                          =       (0:rt_median);
% target_time_retreival               =       (rt_median+1:retrieval_median+rt_median);
% target_time_postpoke            =       (retrieval_median+rt_median+1:retrieval_median+rt_median+post__*1000);
% target_time                             =       [target_time_trigger target_time_rt target_time_retreival target_time_postpoke];
% release_sdfs_warped              =        [];
% 
% for ii =1:length(release_sdfs_all)
%     t_ii                =   release_sdfs_all{ii}(:, 1);
%     sdf_ii              =   release_sdfs_all{ii}(:, 2);
%     ii_rt                = release_event_seqs_all(ii, 2);
%     ii_retrieval        = release_event_seqs_all(ii, 3);
% 
%     ind_preTrigger      = find(t_ii<=target_time_trigger(end));
%     sdf_ii_preTrigger   = sdf_ii(ind_preTrigger); % won't warp
%     ind_Release         = find(t_ii>=0 & t_ii<=ii_rt);
%     t_ii_Release        = t_ii(ind_Release);
%     sdf_ii_Release      = sdf_ii(ind_Release); % warp to target_time_rt
%     % function vw = warp_sdf(t1, v1, tw)
%     sdf_ii_Release_warped = Spikes.SRT.warp_sdf(t_ii_Release, sdf_ii_Release, target_time_rt);
%     ind_Poke            = find(t_ii>=ii_rt+1 & t_ii<=ii_rt+ii_retrieval);
%     t_ii_Poke           = t_ii(ind_Poke);
%     sdf_ii_Poke         = sdf_ii(ind_Poke); % warp to target_time_rt
%     % function vw = warp_sdf(t1, v1, tw)
%     sdf_ii_Poke_warped = Spikes.SRT.warp_sdf(t_ii_Poke, sdf_ii_Poke, target_time_retreival);
%     % post poke activity, not warped.
%     ind_PostPoke            = find(t_ii>(ii_retrieval+ii_rt));
%     ind_PostPoke            = ind_PostPoke(1:length(target_time_postpoke));
%     sdf_ii_PostPoke         = sdf_ii(ind_PostPoke); % warp to target_time_rt
%     % put them together
%     sdf_ii_warped = [sdf_ii_preTrigger; sdf_ii_Release_warped'; sdf_ii_Poke_warped'; sdf_ii_PostPoke];
%     release_sdfs_warped = [release_sdfs_warped sdf_ii_warped];
% end
% 
% warp_out.release.sdf_warped.time = target_time';
% warp_out.release.sdf_warped.sdf_mean = mean(release_sdfs_warped, 2)';
% warp_out.release.sdf_warped.sdf_ci = transpose(bootci(1000, @mean, release_sdfs_warped'));
% warp_out.release.sdf_warped.sdf_trials = release_sdfs_warped;
% warp_out.release.sdf_warped.rt_retrieval_time = [rt_median, retrieval_median];

%%
this_unit = r.BehaviorClass.Subject + '|Session' + r.BehaviorClass.Session + '|Ch' + num2str(id(1)) + '|Unit' + num2str(id(2));
% uicontrol('Style', 'text', 'unit', 'normalized', 'Position', [.01 .9 .8 .08], 'String', this_unit, 'Fontname', 'dejavu sans',...
%     'fontsize', 12, 'fontweight', 'bold', 'backgroundcolor', 'w')

%  Pack output
WarpOut.meta.unit_name        = this_unit;
WarpOut.meta.subject          = r.BehaviorClass.Subject;
WarpOut.meta.session          = r.BehaviorClass.Session;
WarpOut.meta.Channel          = num2str(id(1));
WarpOut.meta.Unit             = num2str(id(2));
% press_psth = {t_spkmat, sdf_press_mean, sdf_press_ci};
% psth_info_explained = 'time, sdf_mean, sdf_ci';
% WarpOut.psth_explained = psth_info_explained;
% WarpOut.press_psth = press_psth;
% WarpOut.press_raster = centin_raster;
% WarpOut.release_psth = release_psth;
% WarpOut.release_raster = release_raster;
% WarpOut.poke_psth = poke_psth;
% WarpOut.poke_raster = poke_raster;
WarpOut.spike_waveform        = {tspk, tspk, spk_mean, spk_std};
WarpOut.spikes_trials         = spikes_trials;
WarpOut.spike_train_explained = {'[cent-in, cent-out, choice] times', 'spk times', 'FPs', 'Ports'};
WarpOut.twarp                 = t_warped;
WarpOut.sdf_all               = ioc_spk_warped;
WarpOut.sdf_avg               = sdf_warped_mean;
WarpOut.sdf_ci                = sdf_warped_ci;
WarpOut.median_time_points    = warped_time_points;
WarpOut.pre_post_durations    = [pre_, post_, post_keep];
% WarpOut.warp2 = warp_out;
% 
% styleAllAxesInFigureGeneral(gcf)
% 
% if toplot
%     % save this figure
%     anm_name             =     r.BehaviorClass.Subject;
%     session              =     r.BehaviorClass.Date;
%     thisFolder = fullfile(pwd, 'Figures_WarpedPSTHs');
%     if ~exist(thisFolder, 'dir')
%         mkdir(thisFolder)
%     end
%     tosavename2= fullfile(thisFolder, [anm_name '_' session '_Ch'  num2str(id(1)) '_Unit' num2str(id(2)) '_Lite']);
%     print (gcf,'-dpng', tosavename2)
% end
% 
% close(hf)
end