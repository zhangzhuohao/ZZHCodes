
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                               %
% SCRIPT FOR FITTING THE EX-GAUSSIAN &          %
% SHIFTED-WALD TO RT DATA                       % 
%                                               %
% Fits the data of individual participants      %
% as well as the group-level data               %
%                                               %
% Mark Hurlstone                                %
% Lancaster University                          %
% m.hurlstone@lancaster.ac.uk                   %
% 21-05-21                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET NUMBER OF PARTICIPANTS, CONDITIONS, AND OBSERVATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exitFlag = 0; count = 0;
while (exitFlag < 1)
    count = count + 1;
    par = ['Participant_' num2str(count) '.txt'];
    if (~exist(par,'file'))
        nPars = count-1;
        dataSmpl = dlmread(['Participant_' num2str(nPars) '.txt']);
        nObs = size(dataSmpl,1);
        nCons = size(dataSmpl,2);
        exitFlag = 1;
    end
end     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTION PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What distribution are we fitting?
chooseDistrib = 0; % 0 = ex-Gaussian; 1 = Shifted-Wald

% Distribution parameters
if (chooseDistrib == 0)
    % ex-Gaussian parameters
    p.tau = 250;
    p.mu  = 500;
    p.sig = 100;
else
    % Shifted-Wald parameters
    p.alpha = 70;
    p.theta = 625;
    p.gamma = 0.2;
end
    
% Convert parameters for fminsearch
tempParms = struct2cell(p); 
for i = 1:length(tempParms)
	parms(i) = tempParms{i};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET-UP FMINSEARCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fminsearch options
x          = parms; 
nfuncevals = 500;
tolerance  = .1;
defopts    = optimset ('fminsearch');
options    = optimset (defopts,'Display','iter','TolFun',...
    tolerance,'MaxFunEvals', nfuncevals);

% Set-up an array for storing best-fitting parameters 
nParms    = length(parms);
bestParms = zeros(nPars+1,nParms*nCons);

% Initialise an aggregate data file 
allData  = zeros(nObs,nCons,nPars);

% zTransorm of response times for creating group distribution
zTrans = zeros(nObs,nCons,nPars);

% Create an array for storing log-likelihoods and chi-square GOFs
lnLs = zeros(nPars+1,nCons);
chi2 = zeros(nPars+1,nCons);

% Create an array for storing results of Kolmogorov-Smirnoff tests on data
ks = zeros(nPars+1,3*nCons); 

% Determine quantile averages for defining observed RT bins
cumProb = .05:.05:.95;                    
edges   = zeros(length(cumProb),nCons);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT CHOSEN DISTRIBUTION TO EACH PARTICIPANT AND CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pars = 1:nPars+1 
    
    if (pars < nPars+1)
        % Retrieve current participant data
        parData = dlmread(['Participant_' num2str(pars) '.txt']);
    else
        % Now all the individual participant data has been fitted, shift 
        % the z-Transformed data of all participants by the average mean.
        % This is used to represent and fit the group data.
        meanMean = zeros(nObs,nCons,nPars);
        meanStd = zeros(nObs,nCons,nPars);
        for pars2 = 1:nPars
            meanMean(:,:,pars2) = repmat(mean(mean(allData),3),nObs,1);
            meanStd(:,:,pars2)  = repmat(mean(std(allData),3),nObs,1);
        end
        groupData = meanStd .* zTrans + meanMean;
    end
       
    % Initialise indexing variables for storing best-fitting parameters
    a = 1; b = nParms;
    
    % Fit the data for each condition
    for con = 1:nCons
        
        if (pars < nPars+1)
            % Retrieve current condition data for *CURRENT* participant
            conData = parData(:,con); 
        else
            % Retrieve current condition data for *ALL* participants
            conData = reshape(squeeze(groupData(:,con,:)), ...
                size(groupData,1)*nPars,1); 
        end
        
        % Get rid of zeros and NaNs     
        noZeros = sort(nonzeros(conData));    
        conData = noZeros(~isnan(noZeros(:)));
        
        % Run Kolmogorov-Smirnoff test
        [h,p,ksstat] = kstest(conData); 
        ks(pars,a:b) = [h,p,ksstat];     

        % Get RT bin edges  
        edges(:,con) = edges(:,con) + quantile(conData,cumProb)';
                
        % Build an aggregate data file as we go along
        if (pars < nPars+1),allData(:,con,pars) = conData; end
    
        % z-Transform of all response times
        if (pars < nPars+1), zTrans(:,con,pars)=(conData-mean(conData))...
                /std(conData); end
        
        % Generate starting parameters based on to-be-fitted data
        [startVec1,startVec2,startVec3,lB,uB] = ...
            genStartingParameters(conData,chooseDistrib);
        
        % Initiate distribution fitting
        [x,fval,dummy,output] = wrapperLoopFmin(parms,conData,...
            startVec1,startVec2,startVec3,lB,uB,chooseDistrib);
        
        % Store best-fitting parameters
        bestParms(pars,a:b) = x;
        
        % Get log-likelihoods
        if (chooseDistrib == 0)
            pdf = exGaussPdf(x,conData); 
        else
            pdf = shiftWaldPdf(x,conData);
        end
        % Log-Likelihood estimate
        lnLs(pars,con) = -sum(log(pdf)); 
        
        % Get chi-square goodness-of-fit statistic
        chi2(pars,con) = chiSquare(x,conData,quantile(conData,...
            cumProb),chooseDistrib);
        
        % Update indexing variables
        a = a + nParms; b = b + nParms;
        
    end
    
end
% Quantile averaging
edges = edges./nPars; % *Not currently used*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE RESULTS TO FILE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rows 1 to N (where N is the number of participants) represent estimates 
% obtained for individual participants, whereas row N+1 provides estimates 
% for the group as a whole.
dlmwrite('Kolmogorov_Smirnov_Test.txt',ks,'delimiter','\t');
dlmwrite('Fitted_Participant_Parameters.txt',bestParms,'delimiter','\t');
dlmwrite('Log_Maximum_Likelihoods.txt',lnLs,'delimiter','\t');
dlmwrite('ChiSquare_Goodness_of_Fit.txt',chi2,'delimiter','\t');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DATA AND BEST-FITTING EX-GAUSSIAN PDFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The histogram plots of the data are produced by creating a group 
% distribution through z-Transformation of the raw response times and
% re-scaling back to the original scale. However, we need to make a
% decision regarding whether to plot the best-fitting parameter averaged 
% ex-Gaussian PDF for the fits to individual participants, or the 
% best-fitting ex-Gaussian PDF for the fit to the group distribution.

% Get best-fitting parameters by condition
pdfPlot = 1; % 0 = mean individual-level parms; 1 = group-level parms
if (pdfPlot == 0)
    bestParmsPlot = mean(bestParms(1:nPars,:),1);
else
    bestParmsPlot = bestParms(nPars+1,:);
end

% Determine number of rows (assumes 3 plots per column)
nRows = ceil(nCons/3); xy = zeros(1,nCons);
 
% Initialise indexing variables for retrieving best-fitting parameters
a = 1; b = nParms;
 
% Number of histogram bins
nBins = 20;

% Get maximum response time
maxRT = max(max(max((allData))));

% Generate sub-plots for each condition
for con = 1:nCons

    xy(con) = subplot(nRows,3,con,'fontSize',12); 
    cla; hold on;
    h     = histogram(groupData(:,con,:),nBins,'FaceColor','w',...
        'Normalization','pdf');
    p1    = bestParmsPlot(a);   
    p2    = bestParmsPlot(a+1); 
    p3    = bestParmsPlot(b);   
    parms = [p1 p2 p3];       
    x     = linspace(0,max(max(max(sort(allData(:,con,:))))),2000'); 
    % Retrieve the appropriate PDF
    if (chooseDistrib == 0)
        f = exGaussPdf(parms,x);
        title(sprintf('Tau = %.2f, Mu = %.2f, Sigma = %.2f',parms),...
            'fontsize',8);
    else
        f = shiftWaldPdf(parms,x);
        title(sprintf('Alpha = %.2f, Theta = %.2f, Gamma = %.2f',parms),...
            'fontsize',8);
    end   
    fNorm = numel(x)*f*h.BinWidth; % Not currently used
    plot(x,f,'r-','LineWidth',2.5);
    xlabel('Response Time (ms)','fontsize',8);
    ylabel('Density','fontsize',8)
    xlim([0 ceil(maxRT/1000)*1000])
    ylim([0 max(f)*1.5])
    a = a + nParms; b = b + nParms;

end
linkaxes(xy,'xy');
