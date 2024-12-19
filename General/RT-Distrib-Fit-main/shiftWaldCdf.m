function f = shiftWaldCdf(x,parms)

% Retrieve parameters
alpha = parms(1); 
theta = parms(2); 
gamma = parms(3); 

% Calculate cumulative distribution of shifted-Wald
x = x-theta; x(x<0) = 0; 
f = normcdf((gamma.*x-alpha)./sqrt(x))+exp(2.*alpha.*gamma).* ...
    normcdf(-((gamma.*x+alpha)./sqrt(x))); 

% REFERENCES 
% Heathcote (2004) Behaviour Research Methods, vol. 36 no. 4, pp
% 678-694.
% Note: Heathcote denotes alpha, a, gamma, m, theta, s, and x, w. 