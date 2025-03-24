function [B, TFrm, TFoutlier, L, U, C] = rmoutliers_tail(data)

[~, ~, ~, L, U0, C] = rmoutliers(data, 'median');

Gap = U0 - C;
U = C + 2*Gap;

TFrm = data>=L & data<=U;
TFoutlier = ~TFrm;
B = data(TFrm);

end