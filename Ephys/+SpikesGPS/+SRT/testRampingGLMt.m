function stat_out = testRampingGLMt(t_spike, ephys_table)

% Use those two fields in r. For faster calculation.

% time bin
bin_width = 0.1;
time_bin = 0:bin_width:1;
time_pos = movmean(time_bin, 2, 'Endpoints', 'discard');
num_bin  = length(time_pos);

% use 0.5 s to separate holding period
time_sep = 0.5;

% initialize spike count matrix
num_trial = height(ephys_table);
spk_count = zeros(num_trial, num_bin);

% get spike count in each bin
for i = 1:num_bin
%     spk_out = getSpikeCount(rs, 'tCentIn', time_bin(i+[0 1])*1000, ku);
    spk_out = getSpikeTimingsWithin(t_spike, ephys_table.tCentIn, time_bin(i+[0 1])*1000);
    spk_count(:,i) = cellfun(@length, spk_out);
end

% spike count threshold
threshold = 20;

% initialize output
stat_out = struct();
% stat_out.unit_id = repmat(ku, 2, 1);
stat_out.port = nan(2,1);
stat_out.status = strings(2,1);
stat_out.mdl  = cell(2,1);
% pre
stat_out.pval_pre = nan(2,1);
stat_out.beta_pre = nan(2,1);
stat_out.se_pre   = nan(2,1);
stat_out.beta_pre_ci_l = nan(2,1);
stat_out.beta_pre_ci_u = nan(2,1);
% post
stat_out.pval_post = nan(2,1);
stat_out.beta_post = nan(2,1);
stat_out.se_post   = nan(2,1);
stat_out.beta_post_ci_l = nan(2,1);
stat_out.beta_post_ci_u = nan(2,1);
% interact
stat_out.pval_int = nan(2,1);
stat_out.beta_int = nan(2,1);
stat_out.se_int   = nan(2,1);
stat_out.beta_int_ci_l = nan(2,1);
stat_out.beta_int_ci_u = nan(2,1);

% use trials whose FP >= 1 sec and HD >= 1 sec (including Correct, Late, Probe <and Wrong> trials)
id_trial = ephys_table.FP>=1 & ephys_table.HD>=1000;

for p = 1:2
    % test for each port
    id_trial_p = id_trial & ephys_table.PortCorrect==p;
    stat_out.port(p) = p;

    % make spike table for GLM
    num_trial_this = sum(id_trial_p);
    spk_count_this = spk_count(id_trial_p, :);
    spk_time_this  = repmat(time_pos, num_trial_this, 1);
    trial_id_this  = repmat(ephys_table.Trials(id_trial_p), 1, num_bin);

    spk_tbl = table(spk_count_this(:), spk_time_this(:), categorical(trial_id_this(:)), 'VariableNames', {'SpikeCount', 'Time', 'Trial'});
    spk_tbl.Period = categorical(spk_tbl.Time < time_sep, [true false], {'Pre', 'Post'});

    % check unit activity
    if sum(spk_tbl.SpikeCount) < threshold
        stat_out.status(p) = "low_activity";
        continue
    end

    % fit GLM with poisson
    try
        lastwarn(''); % clear warnings
        opts = statset('MaxIter', 1000, 'TolFun', 1e-8); % Stricter tolerance
        mdl = fitglm(spk_tbl, 'SpikeCount ~ 1 + Time*Period', ...
            'Distribution', 'poisson', 'Link', 'log', 'Offset', log(bin_width), ...
            'Options', opts);

        % check for warnings indicating convergence issues
        [wmsg, wid] = lastwarn;
        if ~isempty(wid) && (contains(wmsg, 'converge') || contains(wmsg, 'ill-conditioned') || contains(wmsg, 'Iteration limit'))
            stat_out.status(p) = "warning";
            continue
        end

        stat_out.status(p) = "glm_fitted";
        stat_out.mdl{p} = mdl;
    catch ME
        disp(['GLM failed for ' event ': ' ME.message '. Skipping.']);
        stat_out.status(p) = "fit_failed";
    end

    % disp(mdl.Coefficients);

    z_crit = norminv(0.975);
    % post-hoc
    coef_name = mdl.CoefficientNames';
    beta = mdl.Coefficients.Estimate;
    cov_mat = mdl.CoefficientCovariance;
    % find index
    idx_time = find(strcmp(coef_name, 'Time'));
    idx_int  = find(contains(coef_name, 'Time:Period_Post'));

    % pre
    contrast_pre = zeros(1, length(beta));
    contrast_pre(idx_time) = 1;
    % p value
    stat_out.pval_pre(p) = coefTest(mdl, contrast_pre);
    % estimation (log scale)
    est = contrast_pre * beta;
    se  = sqrt(contrast_pre * cov_mat * contrast_pre');
    ci_l = est - z_crit*se;
    ci_u = est + z_crit*se;
    % transform to original scale
    stat_out.beta_pre(p) = est;
    stat_out.se_pre(p)   = se;
    stat_out.beta_pre_ci_l(p) = ci_l;
    stat_out.beta_pre_ci_u(p) = ci_u;

    % post
    contrast_post = zeros(1, length(beta));
    contrast_post([idx_time idx_int]) = 1;
    % p value
    stat_out.pval_post(p) = coefTest(mdl, contrast_post);
    % estimation (log scale)
    est = contrast_post * beta;
    se  = sqrt(contrast_post * cov_mat * contrast_post');
    ci_l = est - z_crit*se;
    ci_u = est + z_crit*se;
    % transform to original scale
    stat_out.beta_post(p) = est;
    stat_out.se_post(p)   = se;
    stat_out.beta_post_ci_l(p) = ci_l;
    stat_out.beta_post_ci_u(p) = ci_u;

    % interaction
    contrast_int = zeros(1, length(beta));
    contrast_int(idx_int) = 1;
    % p value
    stat_out.pval_int(p) = coefTest(mdl, contrast_int);
    % estimation (log scale)
    est = contrast_int * beta;
    se  = sqrt(contrast_int * cov_mat * contrast_int');
    ci_l = est - z_crit*se;
    ci_u = est + z_crit*se;
    % transform to original scale
    stat_out.beta_int(p) = est;
    stat_out.se_int(p)   = se;
    stat_out.beta_int_ci_l(p) = ci_l;
    stat_out.beta_int_ci_u(p) = ci_u;

end

stat_out = struct2table(stat_out);

end % testRampingGLM