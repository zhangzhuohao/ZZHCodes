function DurOut = calDur(Dur, FP, varargin)
% RTOut = calRT(hold_dur, forepriod, 'Remove100ms', 1, 'RemoveOutliers', 1, 'ToPlot', 0, 'CalSE', 0)
% RTOut = calRT(reaction_time, [], 'Remove100ms', 1, 'RemoveOutliers', 1, 'ToPlot', 0, 'CalSE', 0)
% JY 8.19.2022
% revised by ZZH 4.06.2023, if set CalSE to 1, add two addtional field to the output
% ZZH 3.25.2024, remove 50ms as very fast response (used to be 100ms)

if isempty(FP)
    RelTime = Dur;
else
    RelTime = Dur-FP;
end

RelTime = RelTime(~isnan(RelTime));

% computing_method = 'cdf'; % could also be mean, ecdf
remove_50ms = 1; % remove very fast responses, which are deemed anticipatory responses
remove_outliers = 1;
toplot = 1;
calse = 1;

if nargin>2
    for i=1:2:size(varargin,2)
        switch varargin{i}
            case {'Remove100ms'}
                remove_50ms = varargin{i+1}; %
            case {'RemoveOutliers'}
                remove_outliers = varargin{i+1}; %
            case {'ToPlot'}
                toplot = varargin{i+1}; %
            case {'Calse', 'CalSE'}
                calse = varargin{i+1}; %
            otherwise
                errordlg('unknown argument')
        end
    end
end

switch remove_50ms
    case 1
        RelTime = RelTime(RelTime>=0.05);
    case 0
        RelTime = RelTime(RelTime>0.0);
end

% remove outliers
if remove_outliers
    RelTimeOrg = RelTime;
    %     thisMAD = mad(RelTime);
    [RelTime, indremoved] = rmoutliers_custome(RelTime);
    if toplot
        figure(10); clf
        set(10, 'Visible', 'on')
        plot(RelTimeOrg, 'ko');
        hold on
        plot(indremoved, RelTimeOrg(indremoved), 'rx', 'linewidth',2,'markersize', 8);
        set(gca, 'yscale','log')
    end
end


DurOut.N = length(RelTime);
if isempty(RelTime) || DurOut.N<5
    DurOut.median = nan;
    DurOut.median_ksdensity = nan;
    DurOut.q1 = nan;
    DurOut.q3 = nan;
    DurOut.mean = nan;
    DurOut.std = nan;
    DurOut.sem = nan;
else
    DurOut.median = median(RelTime, 'omitnan');
    DurOut.median_ksdensity = kscdf_med(RelTime);
    DurOut.q1 = prctile(RelTime, 25);
    DurOut.q3 = prctile(RelTime, 75);
    DurOut.mean = mean(RelTime, 'omitnan');
    DurOut.std = std(RelTime, 'omitnan');
    DurOut.sem = DurOut.std / sqrt(DurOut.N);
    if calse
        DurOut.median_std = std(bootstrp(1000, @(x) median(x, 'omitnan'), RelTime));
        DurOut.median_ksdensity_std = std(bootstrp(1000, @kscdf_med, RelTime));
    end
end
end

function rtout = kscdf_med(rt)
% use a kernel density method to estimate the median
rt_space = linspace(min(rt), max(rt), 100);
f = ksdensity(rt, rt_space, 'function', 'cdf');
% find f near 0.5
ind_med = find(f>=0.5, 1, 'first');
rtout = rt_space(ind_med);

end

