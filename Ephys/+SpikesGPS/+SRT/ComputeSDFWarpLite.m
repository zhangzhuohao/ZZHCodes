function SDFAll = ComputeSDFWarpLite(spktimes, beh_info, cal_ci, sigma, dt, type)

if nargin < 3
    cal_ci = 1;
    sigma = 20;
    dt = 10;
    type = 'gaussian';
elseif nargin<4
    sigma = 20;
    dt = 10;
    type = 'gaussian';
elseif  nargin<5
    dt = 10;
    type = 'gaussian';
elseif nargin<6
    type = 'gaussian';
end

%% Task information
FPs = unique(beh_info.FP(beh_info.Stage==1));
FPs = FPs(~isinf(FPs));
NumFPs = length(FPs);

Ports = unique(beh_info.PortCorrect);
NumPorts = length(Ports);

%% Gather all
SDFAll = SpikesGPS.ComputeSDFAllLite(spktimes, beh_info, sigma, dt, type);
% for i = 1:length(SDFAll.sdf)
%     [SDFAll.sdf{i}, SDFAll.t_sdf{i}] = downsample_sdf(SDFAll.sdf{i}, SDFAll.t_sdf{i}, dt);
% end
SDFTable = struct2table(SDFAll);

%% Warp sdf for each condition
% warp correct trials only
SDFCorrect = SDFTable(SDFAll.Outcome=="Correct", :);

% the whole duration of time warp
pre_  = 2.5; % 2.5 sec before cent-in, about the lower boundary of shuttle time
post_ = 3; % take 3 sec after choice at first
post_keep = 2.5; % should be shorter than post_, to avoid length error

latency = 2.5; % max cent-out to choice latency (max movement time)

% trim whole trial sdf to time warp duration
ind_valid = find(SDFCorrect.MT<latency*1000 & SDFCorrect.RT>dt);
SDFCorrect = SDFCorrect(ind_valid, :);
trial_dur = num2cell(SDFCorrect.t_choice - SDFCorrect.t_centin);
SDFCorrect.sdf = cellfun(@(x,t,t_p) x(t>=-pre_*1000 & t<=t_p+post_*1000), SDFCorrect.sdf, SDFCorrect.t_sdf, trial_dur, 'UniformOutput', false);
SDFCorrect.t_sdf = cellfun(@(t,t_p) t(t>=-pre_*1000 & t<=t_p+post_*1000), SDFCorrect.t_sdf, trial_dur, 'UniformOutput', false);

%% Warp for each condition
SDFCorrect.t_point  = cell(length(ind_valid), 1);
SDFCorrect.t_warp   = cell(length(ind_valid), 1);
SDFCorrect.sdf_warp = cell(length(ind_valid), 1);

t_points = cell(NumFPs, NumPorts);
sdf_warp = cell(NumFPs, NumPorts);
t_warped = cell(NumFPs, NumPorts);
for i = 1:NumFPs
    i_fp = FPs(i);
    for j = 1:NumPorts
        j_port = Ports(j);

        ind_ij = find(SDFCorrect.FP==i_fp & SDFCorrect.Port==j_port);
        if isempty(ind_ij)
            continue
        end

        SDF_ij = SDFCorrect(ind_ij, :);
        median_hold_duration = median(SDF_ij.HD, 'omitnan'); % median, in ms
        median_movement_time = median(SDF_ij.MT, 'omitnan'); % median, in ms
        % cent-in, trigger, cent-out, choice (relative to cent-in moment)
        t_points{i,j} = [0 1000*i_fp median_hold_duration median_hold_duration+median_movement_time];
        t_temp = [-pre_*1000 t_points{i,j} t_points{i,j}(end)+post_keep*1000];

        t_target = t_temp(1)+dt/2:dt:t_temp(end)-dt/2;
        sdf_warp{i,j} = nan(length(ind_ij), length(t_target));
        for k = 1:length(ind_ij) % for each trial
            % time points in this trial normalize so that the first point is time 0
            k_t_point = [SDF_ij.t_centin(k) SDF_ij.t_trigger(k) SDF_ij.t_centout(k) SDF_ij.t_choice(k)] - SDF_ij.t_centin(k);
            k_t_point = [-pre_*1000 k_t_point k_t_point(end)+post_keep*1000];

            k_t_sdf = SDF_ij.t_sdf{k}; % in ms
            k_sdf   = SDF_ij.sdf{k}; % this is the sdf

            k_sdf_warp = cell(1, length(t_temp)-1);
            for m = 1:length(t_temp)-1
                m_t_temp = t_temp([m m+1]);
                m_t_warp = t_target(t_target>=m_t_temp(1) & t_target<m_t_temp(2));

                m_t_point = k_t_point([m m+1]);
                m_sdf     = k_sdf(k_t_sdf>=m_t_point(1) & k_t_sdf<m_t_point(2));
                m_t_sdf   = k_t_sdf(k_t_sdf>=m_t_point(1) & k_t_sdf<m_t_point(2));

                m_sdf_warp = warp_sdf(m_t_sdf, m_sdf, m_t_warp);

                k_sdf_warp{m} = m_sdf_warp;
            end
            sdf_warp{i,j}(k,:) = cat(2, k_sdf_warp{:});

            SDFCorrect.t_point{ind_ij(k)}  = t_points{i,j};
            SDFCorrect.t_warp{ind_ij(k)}   = t_target;
            SDFCorrect.sdf_warp{ind_ij(k)} = cat(2, k_sdf_warp{:});
        end
        t_warped{i,j} = t_target;
    end
end

% Calculate mean and CI of warped sdf
sdf_warped_ci = cell(NumFPs, NumPorts);
sdf_warped_m  = cell(NumFPs, NumPorts);
for i = 1:NumFPs
    for j = 1:NumPorts
        sdf_warped_m{i,j} = mean(sdf_warp{i,j}, 1);
        if cal_ci
            if size(sdf_warp{i,j}, 1)>=5
                sdf_warped_ci{i,j} = bootci(1000, @mean, sdf_warp{i,j});
            else
                sdf_warped_ci{i,j} = nan(2, length(sdf_warped_m{i,j}));
            end
        end
    end
end

SDFWarp.t_points = t_points;
SDFWarp.t_points_description = ["CentIn", "Trigger", "CentOut", "Choice"];
SDFWarp.pre_post_durations = [pre_, post_, post_keep];

SDFWarp.t_warp   = t_warped;
SDFWarp.sdf_warp = sdf_warp;
SDFWarp.sdf_mean = sdf_warped_m;
SDFWarp.sdf_ci   = sdf_warped_ci;

%% Another warp method
% Pool all FPs for each Port
% Approach to Cent-In: no change
% Cent-In to Cent-In + 250 ms: no change
% Trigger - 250 ms to Trigger: no change
% Trigger to Cent-Out: warp
% Cent-Out to Choice: warp

% around cent-in, not warp
pre_centin  = 2.5; % 2.5 sec before cent-in, about the lower boundary of shuttle time
post_centin = 0.25; % 0.5 sec after cent-in

SDFCorrect.t_pool_centin   = cellfun(@(t) t(t>=-pre_centin*1000 & t<post_centin*1000), SDFCorrect.t_warp, 'UniformOutput', false);
SDFCorrect.sdf_pool_centin = cellfun(@(x,t) x(t>=-pre_centin*1000 & t<post_centin*1000), SDFCorrect.sdf_warp, SDFCorrect.t_warp, 'UniformOutput', false);

pool_t_points.centin = cell(1, NumPorts);
pool_sdf_warp.centin = cell(1, NumPorts);
pool_t_warped.centin = cell(1, NumPorts);
for j = 1:NumPorts
    j_port = Ports(j);
    ind_j = find(SDFCorrect.Port==j_port);
    SDF_j = SDFCorrect(ind_j, :);

    pool_t_points.centin{j} = 0;
    pool_sdf_warp.centin{j} = cat(1, SDFCorrect.sdf_pool_centin{ind_j});
    pool_t_warped.centin{j} = SDFCorrect.t_pool_centin{ind_j(1)};
end

% around trigger, warp
pre_trigger = 0.25;
post_choice = 2.5; % take 3 sec after choice at first

t_point_trig = cell(1,NumPorts);
pool_t_points.trigger = cell(1, NumPorts);
pool_sdf_warp.trigger = cell(1, NumPorts);
pool_t_warped.trigger = cell(1, NumPorts);
for j = 1:NumPorts
    j_port = Ports(j);
    ind_j = find(SDFCorrect.Port==j_port);
    SDF_j = SDFCorrect(ind_j, :);

    median_reaction_time = median(SDF_j.RT, 'omitnan'); % median, in ms
    median_movement_time = median(SDF_j.MT, 'omitnan'); % median, in ms
    % trigger, cent-out, choice (relative to trigger moment)
    t_point_trig{j} = [0 median_reaction_time median_reaction_time+median_movement_time];
    t_temp = [-pre_trigger*1000 t_point_trig{j} t_point_trig{j}(end)+post_choice*1000];

    t_target = t_temp(1)+dt/2:dt:t_temp(end)-dt/2;
    for k = 1:length(ind_j) % for each trial
        % time points in this trial normalize so that the first point is time 0
        k_t_point = [SDF_j.t_trigger(k) SDF_j.t_centout(k) SDF_j.t_choice(k)] - SDF_j.t_trigger(k);
        k_t_point = [-pre_trigger*1000 k_t_point k_t_point(end)+post_choice*1000];

        k_t_sdf = SDF_j.t_sdf{k} - (SDF_j.t_trigger(k) - SDF_j.t_centin(k)); % in ms
        k_sdf   = SDF_j.sdf{k}; % this is the sdf

        k_sdf_warp = cell(1, length(t_temp)-1);
        for m = 1:length(t_temp)-1
            m_t_temp = t_temp([m m+1]);
            m_t_warp = t_target(t_target>=m_t_temp(1) & t_target<m_t_temp(2));

            m_t_point = k_t_point([m m+1]);
            m_sdf     = k_sdf(k_t_sdf>=m_t_point(1) & k_t_sdf<m_t_point(2));
            m_t_sdf   = k_t_sdf(k_t_sdf>=m_t_point(1) & k_t_sdf<m_t_point(2));

            m_sdf_warp = warp_sdf(m_t_sdf, m_sdf, m_t_warp);
            k_sdf_warp{m} = m_sdf_warp;
        end
        SDFCorrect.t_point_trigger{ind_j(k)}  = t_point_trig{j};
        SDFCorrect.t_pool_trigger{ind_j(k)}   = t_target;
        SDFCorrect.sdf_pool_trigger{ind_j(k)} = cat(2, k_sdf_warp{:});
    end
    pool_t_points.trigger{j} = t_point_trig{j};
    pool_sdf_warp.trigger{j} = cat(1, SDFCorrect.sdf_pool_trigger{ind_j});
    pool_t_warped.trigger{j} = SDFCorrect.t_pool_trigger{ind_j(1)};
end

% Calculate mean and CI of warped sdf
pool_sdf_warped_mean.centin  = cell(1, NumPorts);
pool_sdf_warped_ci.centin    = cell(1, NumPorts);
pool_sdf_warped_mean.trigger = cell(1, NumPorts);
pool_sdf_warped_ci.trigger   = cell(1, NumPorts);
for j = 1:NumPorts
    % cent-in around
    pool_sdf_warped_mean.centin{j} = mean(pool_sdf_warp.centin{j}, 1);
    if cal_ci
        if size(pool_sdf_warp.centin{j}, 1)>=5
            pool_sdf_warped_ci.centin{j} = bootci(1000, @mean, pool_sdf_warp.centin{j});
        else
            pool_sdf_warped_ci.centin{j} = nan(2, length(pool_sdf_warped_mean.centin{j}));
        end
    end

    % trigger around
    pool_sdf_warped_mean.trigger{j} = mean(pool_sdf_warp.trigger{j}, 1);
    if cal_ci
        if size(pool_sdf_warp.trigger{j}, 1)>=5
            pool_sdf_warped_ci.trigger{j} = bootci(1000, @mean, pool_sdf_warp.trigger{j});
        else
            pool_sdf_warped_ci.trigger{j} = nan(2, length(pool_sdf_warped_mean.trigger{j}));
        end
    end
end

SDFWarpPool.t_points = pool_t_points;
SDFWarpPool.t_points_description.trigger = ["Trigger", "CentOut", "Choice"];
SDFWarpPool.t_points_description.centin  = "CentIn";
SDFWarpPool.pre_post_durations.trigger = [pre_trigger, post_choice, post_keep];
SDFWarpPool.pre_post_durations.centin  = [pre_centin, post_centin];

SDFWarpPool.t_warp   = pool_t_warped;
SDFWarpPool.sdf_mean = pool_sdf_warped_mean;
SDFWarpPool.sdf_ci   = pool_sdf_warped_ci;

%% Save sdf of correct trials
SDFAll.SDFCorrect = SDFCorrect;
SDFAll.SDFWarp = SDFWarp;
SDFAll.SDFWarpPool = SDFWarpPool;

end % ComputeSDFWarpLite
