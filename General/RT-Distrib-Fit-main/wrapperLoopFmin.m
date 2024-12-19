function [bestX,bestFval,bestDummy,bestOutput] = ...
    wrapperLoopFmin(parms,conData,startVec1,startVec2,startVec3,...
    lB,uB,chooseDistrib)

bestFval   = realmax;
bestX      = parms.*0;
bestDummy  = 0;
bestOutput = 0;

for p1 = startVec1         % tau | alpha
    for p2 = startVec2     %  mu | theta
        for p3 = startVec3 % sig | gamma
            [x,fval,dummy,output] = fminSearchBnd(@logMaxLikelihood,...
                [p1 p2 p3],lB,uB);
            if fval < bestFval
                bestFval = fval;
                bestX = x;
                bestDummy = dummy;
                bestOutput = output;
            end
        end
    end
end

% Objective function for calculating log-maximum likelihood
    function lnL = logMaxLikelihood(parms)
        % Retrieve appropriate PDF
        if (chooseDistrib == 0)
            pdf = exGaussPdf(parms,conData);
        else
            pdf = shiftWaldPdf(parms,conData);
        end
        % The PDF for the shifted Wald occassionally contains missing
        % values - make sure they are removed
        pdf(isnan(pdf)) = 0;

        % Make sure we don't have any zeros in PDF
        pdf(pdf<eps) = eps;

        % Log-Likelihood estimate
        lnL  = -sum(log(pdf));

    end
end