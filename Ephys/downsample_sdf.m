function [sdf_out, t_sdf_out] = downsample_sdf(sdf_in, t_sdf, ds_ratio)

n_bins = floor(length(t_sdf) / ds_ratio);

t_sdf = t_sdf(1:n_bins*ds_ratio);
sdf_in   = sdf_in(1:n_bins*ds_ratio);

t_sdf_mat = reshape(t_sdf, ds_ratio, n_bins);
t_sdf_out = median(t_sdf_mat, 1);

sdf_mat = reshape(sdf_in, ds_ratio, n_bins);
sdf_out = mean(sdf_mat, 1);

end

