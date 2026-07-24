function [aucs, se, ci, pval] = cal_auc_batch(labels, scores_matrix, label_pos, varargin)
    % auc_batch - Fast vectorized computation of AUC for multiple score columns sharing identical labels.
    %
    % Inputs:
    %   labels        : N x 1 vector, binary label
    %   scores_matrix : N x P matrix (where N is samples and P is timepoints/features/variables/channels)
    %   label_pos     : scalar
    %
    % Name-Value Pair Arguments:
    %   'Alpha'       - Significance level for Confidence Intervals (default: 0.05 for 95% CI).
    %   'Tail'        - Hypothesis test tail against chance level (0.5):
    %                   'both'  (default) H0: AUC = 0.5 vs H1: AUC ~= 0.5
    %                   'right'           H0: AUC <= 0.5 vs H1: AUC > 0.5
    %                   'left'            H0: AUC >= 0.5 vs H1: AUC < 0.5
    %
    % Outputs:
    %   aucs          : 1 x P vector containing the calculated AUC value for each column.
    %   ci            - Matrix (2 x P) of Confidence Intervals [lower_bound; upper_bound].
    %   se            - Row vector (1 x P) of DeLong Standard Errors.
    %   p_values      - Row vector (1 x P) of p-values testing against chance level (0.5).

    % Parse optional Name-Value arguments using inputParser
    p = inputParser;
    addRequired(p, 'labels');
    addRequired(p, 'scores_matrix');
    addRequired(p, 'label_pos');
    addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'Tail', 'both', @(x) ismember(lower(x), {'both', 'right', 'left'}));

    parse(p, labels, scores_matrix, label_pos, varargin{:});
    
    alpha = p.Results.Alpha;
    tail = lower(p.Results.Tail);

    if ~ismember(label_pos, labels)
        error('Positive class is not found in the input data.');
    end
    pos_mask = labels==label_pos;
    neg_mask = ~pos_mask;

    n_pos = sum(pos_mask);
    n_neg = sum(neg_mask);

    if n_pos == 0 || n_neg == 0
        aucs = nan(1, size(scores_matrix,2));
        se = nan(1, size(scores_matrix,2));
        ci = nan(2, size(scores_matrix,2));
        pval = nan(1, size(scores_matrix,2));

        return
    end

    % Compute ranks column-wise in parallel (handles ties automatically)
    R_all = tiedrank(scores_matrix); % N x P matrix

    % Extract ranks for positive samples and compute AUC
    R_pos = R_all(pos_mask, :);
    sum_ranks_pos = sum(R_pos, 1); % 1 x P
    aucs = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg);

    % Compute DeLong Standard Errors (SE) and CIs
    if nargout>1
        R_pos_internal = tiedrank(scores_matrix(pos_mask, :)); % n_pos x P
        R_neg_internal = tiedrank(scores_matrix(neg_mask, :)); % n_neg x P

        V10 = (R_pos - R_pos_internal) / n_neg;                   % n_pos x P
        V01 = 1 - (R_all(neg_mask, :) - R_neg_internal) / n_pos;  % n_neg x P

        S10 = var(V10, 0, 1); % 1 x P
        S01 = var(V01, 0, 1); % 1 x P

        se = sqrt(S10 / n_pos + S01 / n_neg);

        % CIs via normal approximation
        z_alpha = norminv(1 - alpha / 2);
        ci = [max(0, aucs - z_alpha .* se); min(1, aucs + z_alpha .* se)];
    end

    % Hypothesis Testing against H0: AUC = 0.5
    if nargout > 3
        se_clean = se;
        se_clean(se_clean == 0) = eps; % Avoid division by zero for constant features

        z_scores = (aucs - 0.5) ./ se_clean;

        switch tail
            case 'both'
                pval = 2 * normcdf(-abs(z_scores));
            case 'right'
                pval = normcdf(-z_scores);
            case 'left'
                pval = normcdf(z_scores);
        end
    end

end