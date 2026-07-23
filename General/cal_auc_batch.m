function aucs = cal_auc_batch(labels, scores_matrix, label_pos)
    % auc_batch - Fast vectorized computation of AUC for multiple score columns sharing identical labels.
    %
    % Inputs:
    %   labels        : N x 1 vector, binary label
    %   scores_matrix : N x P matrix (where N is samples and P is timepoints/features/variables/channels)
    %   label_pos     : scalar
    %
    % Outputs:
    %   aucs          : 1 x P vector containing the calculated AUC value for each column.

    if ~ismember(label_pos, labels)
        error('Positive class is not found in the input data.');
    end
    pos_mask = labels==label_pos;
    n_pos = sum(pos_mask);
    n_neg = length(labels) - n_pos;

    % 1. Compute ranks column-wise in parallel (handles ties automatically)
    R = tiedrank(scores_matrix); % N x P 矩阵

    % 2. Extract ranks for positive samples and sum along each column
    sum_ranks_pos = sum(R(pos_mask, :), 1); % 1 x P

    % 3. Compute AUC directly using the Mann-Whitney U statistic formula
    aucs = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg);
end