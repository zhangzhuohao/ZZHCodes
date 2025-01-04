function Stage = getStage(obj)
%GETSTAGE Summary of this function goes here
%   Detailed explanation goes here
Stage = zeros(length(obj.FP), 1);
Stage(obj.FP==1.5) = 1;
end
