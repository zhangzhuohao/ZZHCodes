function [startVec1,startVec2,startVec3,lB,uB] = ...
    genStartingParameters(data,chooseDistrib)

% Generate starting paramater values for appropriate distribution
if (chooseDistrib == 0)
        
    % Generate sensible starting parameter values for ex-Gaussian
    tau = std(data).*0.8;
    mu  = mean(data)-tau;
    sig = sqrt(var(data)-(tau^2));

    % For each parameter, create a vector of starting parameter values 
    % containing the above value, plus a lower bound that is 50% less than that
    % value, and an upper bound that is 50% higher than that value
    startVec1 = [.5*tau tau tau+tau*.5];
    startVec2 = [.5*mu mu mu+mu*.5];
    startVec3 = [.5*sig sig sig+sig*.5]; 
      
else
   
    % Generate sensible starting parameter values for Shifted Wald
    alpha = 70;
    theta = min(data)-10; % This must be less than the quickest RT
    gamma = 0.2;
    
    % For each parameter, create a vector of starting parameter values in
    % same way as for ex-Gaussian (see above). However, the starting points
    % for theta must all be less than the fastest RT
    startVec1 = [.5*alpha alpha alpha+alpha*.5];
    startVec2 = [theta*.9 theta*.5 theta];
    startVec3 = [.5*gamma gamma gamma+gamma*.5];
    
end

% Generate lower and upper reflection boundaries on parameters
lB = [0 0 0];
uB = [max(data) range(data) range(data)];

% REFERENCES 
% Lacouture & Cousineau (2008) Tutorials in Quantitative Methods for
% Psychology, vol. 4 no. 1, pp 35-45.
% Heathcote (2004) Behaviour Research Methods, vol. 36 no. 4, pp
% 678-694.