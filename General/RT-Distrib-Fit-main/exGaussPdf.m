function f = exGaussPdf(parms,x)

% Retrieve parameters
tau = parms(1);
mu  = parms(2);
sig = parms(3);

% Calculate ex-Gaussian PDF
f = (1/tau)*exp(((mu-x)/tau)+((sig^2)/(2*tau.^2))) .* ...
    .5.*(1+erf((((x-mu)/sig)-(sig/tau))/sqrt(2)));

% REFERENCES 
% Lacouture & Cousineau (2008) Tutorials in Quantitative Methods for
% Psychology, vol. 4 no. 1, pp 35-45.