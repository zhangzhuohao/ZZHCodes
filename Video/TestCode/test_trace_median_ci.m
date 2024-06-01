obj = TrajProgressClass;
%%
info    = obj.TrialInfo;
info    = info(info.Label=="Control", :);
ind     = info.Stage==1 & info.Outcome=="Correct";
t       = obj.TimeMatIn;
mat_ang = obj.TraceMatrix.Control.In.AngleHead;

ang_sorted = obj.sort_data(mat_ang(ind,:), {info.FP(ind), info.PortCorrect(ind)}, {obj.TargetFP, obj.LeftRight});
m = cellfun(@(x, t, t_bound) obj.cal_trace_median(x, t, t_bound, 0.05), ang_sorted, repmat({t}, 3, 2), {[0 500], [0 500]; [0 1000] [0 1000]; [0 1500] [0 1500]}, 'UniformOutput', false);

figure(); hold on;
fill([m{2,1}.time flip(m{2,1}.time)], [m{2,1}.ci(1,:) flip(m{2,1}.ci(2,:))], 'r', 'FaceColor', GPSColor.PortL, 'FaceAlpha', .2, 'EdgeColor', 'none')
fill([m{2,1}.time flip(m{2,1}.time)], [m{2,2}.ci(1,:) flip(m{2,2}.ci(2,:))], 'r', 'FaceColor', GPSColor.PortR, 'FaceAlpha', .2, 'EdgeColor', 'none')

plot(m{2,1}.time, m{2,1}.trace, 'Color', GPSColor.PortL, 'LineWidth', 2);
plot(m{2,2}.time, m{2,2}.trace, 'Color', GPSColor.PortR, 'LineWidth', 2);

%%
info    = obj.TrialInfo;
info    = info(info.Label=="Chemo", :);
ind     = info.Stage==1 & info.Outcome=="Correct";
t       = obj.TimeMatIn;
mat_ang = obj.TraceMatrix.Chemo.In.AngleHead;

ang_sorted = obj.sort_data(mat_ang(ind,:), {info.FP(ind), info.PortCorrect(ind)}, {obj.TargetFP, obj.LeftRight});
m = cellfun(@(x, t, t_bound) obj.cal_trace_median(x, t, t_bound, 0.05), ang_sorted, repmat({t}, 3, 2), {[0 500], [0 500]; [0 1000] [0 1000]; [0 1500] [0 1500]}, 'UniformOutput', false);

figure(); hold on;
fill([m{2,1}.time flip(m{2,1}.time)], [m{2,1}.ci(1,:) flip(m{2,1}.ci(2,:))], 'r', 'FaceColor', GPSColor.PortL, 'FaceAlpha', .2, 'EdgeColor', 'none')
fill([m{2,1}.time flip(m{2,1}.time)], [m{2,2}.ci(1,:) flip(m{2,2}.ci(2,:))], 'r', 'FaceColor', GPSColor.PortR, 'FaceAlpha', .2, 'EdgeColor', 'none')

plot(m{2,1}.time, m{2,1}.trace, 'Color', GPSColor.PortL, 'LineWidth', 2);
plot(m{2,2}.time, m{2,2}.trace, 'Color', GPSColor.PortR, 'LineWidth', 2);