function result = testTrace(time_points, M1, M2, shuffle_iters, alpha, label)

N1 = mean(cellfun(@(x) height(x), M1));
N2 = mean(cellfun(@(x) height(x), M2));

if min([N1 N2]) < 10
    result = [];
    fprintf("Pass.\n")
    return
end
%%

switch label
    case {'Chosen'}
        label = "port_chosen";
    case {'Correct'}
        label = "port_correct";
end

auc = zeros(length(time_points), 1);
sig = zeros(length(time_points), 1);
auc_s_mean = zeros(length(time_points), 1);
auc_s_std  = zeros(length(time_points), 1);
auc_s_sem  = zeros(length(time_points), 1);

for i = 1:length(time_points)
    roc_this = rocmetrics([M1{i}.(label); M2{i}.(label)], [M1{i}.ang; M2{i}.ang], "L");
    auc_this = 0.5 + abs(roc_this.AUC - 0.5);
    auc(i) = auc_this;

    auc_shuffle = zeros(1, shuffle_iters);
    for s = 1:shuffle_iters
        labels = [M1{i}.(label); M2{i}.(label)];
        labels_shuffled = labels(randperm(length(labels)));
        roc_shuffle = rocmetrics(labels_shuffled, [M1{i}.ang; M2{i}.ang], "L");
        auc_shuffle(s) = 0.5 + abs(roc_shuffle.AUC - 0.5);
    end
    auc_s_mean(i) = mean(auc_shuffle);
    auc_s_std(i)  = std(auc_shuffle);
    auc_s_sem(i)  = auc_s_std(i) / sqrt(shuffle_iters);

    sig_this = auc_this >= prctile(auc_shuffle, 100*(1-alpha));
    sig(i) = sig_this;

end
fprintf("Done.\n")

time_points = time_points';
result = table(time_points, auc, sig, auc_s_mean, auc_s_std, auc_s_sem);

end

