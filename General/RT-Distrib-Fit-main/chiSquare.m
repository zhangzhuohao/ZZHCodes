function chi2 = chiSquare(parms,x,edges,chooseDistrib)

% Probability masses for observed response times
histCount = histc(x,[-Inf,edges,Inf]);
histCount = histCount(1:end-1);
pmObs     = histCount./sum(histCount);
pmObs     = pmObs(:);

% Probability masses for predicted response times
if (chooseDistrib == 0)
    pmPrd = diff([0,exGaussCdf(edges,parms),1]);
else
    pmPrd = diff([0,shiftWaldCdf(edges,parms),1]);
end
pmPrd = pmPrd(:);

% Compute chi square statistic
nObs = numel(x);
chi2 = nObs.*sum(((pmObs-pmPrd).^2)./pmPrd);

% References
% Ratcliff & Tuerlinckx (2002) Psychonomic Bulletin & Review, vol. 9 no. 3,
% pp 438-481. 