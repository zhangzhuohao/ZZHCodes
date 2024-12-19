function f = shiftWaldPdf(parms,x)

alpha = parms(1);
theta = parms(2);
gamma = parms(3);

f = (alpha./(sqrt(2.*pi.*(x(x>=theta)-theta).^3))).* ...
    exp(-((alpha-gamma.*(x(x>=theta)-theta)).^2./(2.*(x(x>=theta)-theta))));
prfx = nan(1,sum(x < theta));
if (iscolumn(f) ~= 1), f = f'; end
if (iscolumn(prfx) ~= 1), prfx = prfx'; end
f = [prfx; f];

% REFERENCES 
% Heathcote (2004) Behaviour Research Methods, vol. 36 no. 4, pp
% 678-694.
% Note: Heathcote denotes alpha, a, gamma, m, theta, s, and x, w.