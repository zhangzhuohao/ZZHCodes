function f = exGaussCdf(x,parms)

% Retrieve parameters
tau = parms(1);
mu  = parms(2);
sig = parms(3);

% Calculate cumulative distribution of ex-Gaussian
f = -exp(-x/tau+mu/tau+sig^2/2/tau^2) .* ...
    normcdf((x-mu-sig^2/tau)/sig) + ...
    normcdf((x-mu)/sig);

% REFERENCES 
% Heathcote et al (1991) Quantitative Methods In Psychology, 
% vol. 109 no. 2, pp 340-347.