function [sdf_out, t_sdf] = sdf26(tspk, time_range, sigma, dt)
% spike density function, 2025. JY
% modified by (ZZH, 2026) to calculate sdf for multi trials respectively
%
% Inputs:
%   tspk       - Vector of spike times in milliseconds, or cell for spike times in each trial
%   time_range - [start end] in ms, defines the time window for analysis
%   sigma      - Standard deviation of Gaussian kernel (in ms), controls width
%   dt         - Time step for output sdf_out (in ms), e.g., 1 ms
%
% Output:
%   sdf_out    - Spike Density Function (firing rate in spikes/sec)
%   t          - Time points of sdf
%
% Example:
%   tspk = [100, 150, 200, 250]; % Spike times in ms
%   sdf_out = sdf_out25(tspk, [0 500], 20, 1);

% Input validation
if nargin < 4
    error('All inputs (tspk, time_range, sigma, dt) are required.');
end

if length(time_range) ~= 2 || time_range(1) >= time_range(2)
    error('time_range must be a 2-element vector [start end] with start < end.');
end
if sigma <= 0 || dt <= 0
    error('sigma and dt must be positive.');
end

if iscell(tspk)
    num_trials = length(tspk);
elseif isvector(tspk)
    tspk = {tspk};
    num_trials = 1;
else
    error('tspk must be a vector or a cell array');
end

% Define time vector based on time_range and dt
t_start = floor(time_range(1));
t_end = floor(time_range(2));
t_sdf = t_start:dt:t_end;  % Time points for sdf_out output
bin_edges = t_start-dt/2:dt:t_end+dt/2;

% Initialize sdf_out array
sdf_out = zeros(num_trials, length(t_sdf));

% Gaussian kernel parameters
% Kernel size should cover ±3σ for 99.7% of the Gaussian
kernel_half_width = ceil(4 * sigma);  % In ms
kernel_t = -kernel_half_width:dt:kernel_half_width;
gaussian_kernel = (1 / (sigma * sqrt(2 * pi))) * exp(-kernel_t.^2 / (2 * sigma^2));
% Normalize kernel to ensure area = 1 when scaled to spikes/sec later
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);

for i = 1:num_trials
    tspk_i = tspk{i};
    if isempty(tspk_i)
        continue
    end
    % Convert tspk to spike train
    spk_train = histcounts(tspk_i, bin_edges);
    % Perform convolution
    sdf_out(i,:) = conv(spk_train, gaussian_kernel, 'same');
    % Convert sdf_out to firing rate (spikes/sec)
    sdf_out(i,:) = sdf_out(i,:) * (1000 / dt);  % Scale from per-ms to per-second
end
end