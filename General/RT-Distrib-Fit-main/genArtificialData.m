
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
% SCRIPT FOR PRODUCING ARTIFICAL DATA FOR FITTING BY GENERATING   %
% RESPONSE TIMES USING THE EX-GAUSSIAN OR SHIFTED-WALD            %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars

% Choose distribution from which to generate response times
chooseDistrib = 0; % 0 = ex-Gaussian; 1 = Shifted-Wald

if (chooseDistrib == 0)

    % Assume three conditions with the following expected distribution 
    % parameter values: 
    %   Condition 1: tau = 250; mu = 500; sig = 100;
    %   Condition 2: tau = 300; mu = 600; sig = 100;
    %   Condition 3: tau = 350; mu = 700; sig = 100;
    parms = [250 500 100; 
             300 600 100;
             350 700 100];
        
else
    % Assume three conditions with the following expected distribution 
    % parameter values: 
    %   Condition 1: alpha = 70; theta = 625; gamma = 0.18;
    %   Condition 2: alpha = 45; theta = 625; gamma = 0.16;
    %   Condition 3: alpha = 30; theta = 625; gamma = 0.14;    
    parms = [70 625 0.18; 
             44 625 0.16;
             28 625 0.14];
               
end

% Lower and upper bounds on parameters
lB = parms*.5;
uB = parms*1.5; 
            
% Number of parameters
nParms = size(parms,2);

% Number of participants to generate data for
nPars = 15;

% Number of conditions
nCons = size(parms,1);

% Number of observations to generate per condition
nObs = 150;

% Store randomly generated parameter values
recParms = zeros(nPars,size(parms,2)*nCons);

for i = 1:nPars
    
    data = zeros(nObs,nCons);
    a = 1; b = nParms;
    
    for j = 1:nCons    
        
        % Generate uniformly distributed parameter values
        p1 = (uB(j,1)-lB(j,1)).*rand(1,1)+lB(j,1); % tau | alpha
        p2 = (uB(j,2)-lB(j,2)).*rand(1,1)+lB(j,2); %  mu | theta
        p3 = (uB(j,3)-lB(j,3)).*rand(1,1)+lB(j,3); % sig | gamma
        
        % Generate ex-Gaussian | shifted Wald random numbers according to 
        % distribution parameters        
        if (chooseDistrib == 0)
            data(:,j) = exGaussRanNum(p1,p2,p3,nObs,1);
        else
            data(:,j) = shiftWaldRanNum(p1,p2,p3,nObs,1);
        end
               
        % Store parameter values
        recParms(i,a:b) = [p1 p2 p3];
        a = a + nParms; b = b + nParms;
                    
    end
    % Write the data to file
    dlmwrite(['Participant_' num2str(i) '.txt'],...
        data,'delimiter','\t');
end
% Write the distribution parameters for each participant and condition to
% a common file
dlmwrite('Artificial_Participant_Parameters.txt',recParms,'delimiter','\t');
