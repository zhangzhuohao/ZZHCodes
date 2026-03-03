function PopOut = PopulationActivity(r, varargin)
% 5/13/2023 Jianing Yu
% Build upon SRTSpikesPopulation
% Since we add PSTH to r, we will just use previously computed results to
% build Population activity. 
% we avoid re-compute the same stuff. 

% % if nargin<1
% %     load(Spikes.r_name)
% % end

Units = r.Units.SpikeNotes;

TrialTypes = ["Cued"; "Uncued"; "Premature"];
% Extract these event-related activity
nType = length(TrialTypes);
[~, nPort] = size(r.PSTH.PSTHs(1).CentIn);
% nFP = length(r.BehaviorClass.TargetFP);
% nPort = length(r.BehaviorClass.Ports);

%
PSTH_CentIn     = cell(nType, nPort); % premature, uncued and cued
PSTH_CentInZ    = cell(nType, nPort); % premature, uncued and cued
PSTH_CentInStat = cell(nType, nPort); % this gives the statistics of cent-in

PSTH_CentInAll     = []; % merge premature, uncued and cued
PSTH_CentInAllZ    = []; % merge premature, uncued and cued
PSTH_CentInAllStat = [];

PSTH_CentOut     = cell(nType, nPort);
PSTH_CentOutZ    = cell(nType, nPort);
PSTH_CentOutStat = cell(nType, nPort); % this gives the statistics of cent-out

PSTH_CentOutAll     = []; % merge premature, uncued and cued
PSTH_CentOutAllZ    = []; % merge premature, uncued and cued
PSTH_CentOutAllStat = [];

PSTH_Trigger     = cell(nType, nPort);
PSTH_TriggerZ    = cell(nType, nPort);
PSTH_TriggerStat = cell(nType, nPort); % this gives the statistics of trigger

PSTH_Reward     = cell(nType, nPort);
PSTH_RewardZ    = cell(nType, nPort);
PSTH_RewardStat = cell(nType, nPort); % this gives the statistics of release

if ~isfield(r, 'PSTH')
    error('Compute PSTH first. Run " >>SRTSpikes(r, []);" ')
end

n_unit = length(r.PSTH.PSTHs);

for i = 1:n_unit
    PSTH_baseline = [];
    if i==1 % first row for time info
        PSTH_CentInAll(1, :)  = r.PSTH.PSTHs(i).CentInAll{2};
        PSTH_CentInAllZ(1, :) = r.PSTH.PSTHs(i).CentInAll{2};
    end
    PSTH_CentInAll = [PSTH_CentInAll; r.PSTH.PSTHs(i).CentInAll{1}];
    tspkmat_centinall = r.PSTH.PSTHs(i).CentInAll{4};
    trialspxmat_centinall = r.PSTH.PSTHs(i).CentInAll{3};

    StatOut = ExamineTaskResponsive(tspkmat_centinall, trialspxmat_centinall);
    StatOut.CellIndx = r.Units.SpikeNotes(i, :);
    PSTH_CentInAllStat.StatOut(i) = StatOut;
    PSTH_Trials = zeros(nType, nPort);

    for jtype = 1:nType
        for kport = 1:nPort
            if jtype==3 % premature
                psth_centin_this  = r.PSTH.PSTHs(i).PrematureCentIn{1, kport};
                psth_centout_this = r.PSTH.PSTHs(i).PrematureCentOut{1, kport};
                psth_reward_this  = cell(1,6); % set it empty
            else
                psth_centin_this  = r.PSTH.PSTHs(i).CentIn{jtype, kport};
                psth_centout_this = r.PSTH.PSTHs(i).CentOut{jtype, kport};
                psth_reward_this  = r.PSTH.PSTHs(i).RewardChoice{jtype, kport};
            end
            if jtype==1 % cued
                psth_trig_this = r.PSTH.PSTHs(i).Triggers{jtype, kport};
            else
                psth_trig_this = cell(1,6);
            end

            if i==1 % first row for time info
                PSTH_CentIn{jtype, kport}(1,:)   = psth_centin_this{2};
                PSTH_CentInZ{jtype, kport}(1,:)  = psth_centin_this{2};
                PSTH_CentOut{jtype, kport}(1,:)  = psth_centout_this{2};
                PSTH_CentOutZ{jtype, kport}(1,:) = psth_centout_this{2};
                if jtype~=3 % correct trials
                    PSTH_Reward{jtype, kport}(1,:)  = psth_reward_this{2};
                    PSTH_RewardZ{jtype, kport}(1,:) = psth_reward_this{2};
                end
                if jtype==1 % correct cued trials
                    PSTH_Trigger{jtype, kport}(1,:)  = psth_trig_this{2};
                    PSTH_TriggerZ{jtype, kport}(1,:) = psth_trig_this{2};
                end
            end

            % CentIn
            PSTH_CentIn{jtype, kport} = [PSTH_CentIn{jtype, kport}; psth_centin_this{1}];
            tspkmat     = psth_centin_this{4};
            trialspxmat = psth_centin_this{3};

            % restricting activity to a relatively narrow window.
            tCentInWin = [-2500 1500];
            indWin = find(tspkmat>=tCentInWin(1) & tspkmat<= tCentInWin(2));
            tspkmat = tspkmat(indWin);
            trialspxmat = trialspxmat(indWin, :);

            StatOut = ExamineTaskResponsive(tspkmat, trialspxmat);
            StatOut.CellIndx = r.Units.SpikeNotes(i, :);
            PSTH_CentInStat{jtype, kport}.StatOut(i) = StatOut;
            PSTH_Trials(jtype, kport) = size(trialspxmat, 2);
            PSTH_baseline = [PSTH_baseline psth_centin_this{1}];

            % CentOut
            PSTH_CentOut{jtype, kport} = [PSTH_CentOut{jtype, kport}; psth_centout_this{1}];
            tspkmat     = psth_centout_this{4};
            trialspxmat = psth_centout_this{3};

            % restricting activity to a relatively narrow window.
            tCentOutWin = [-500 1000];
            indWin = find(tspkmat>=tCentOutWin(1) & tspkmat<= tCentOutWin(2));
            tspkmat = tspkmat(indWin);
            trialspxmat = trialspxmat(indWin, :);

            StatOut = ExamineTaskResponsive(tspkmat, trialspxmat);
            StatOut.CellIndx =  r.Units.SpikeNotes(i, :);
            PSTH_CentOutStat{jtype, kport}.StatOut(i) = StatOut;
            PSTH_baseline = [PSTH_baseline psth_centout_this{1}];

            % Trigger
            PSTH_Trigger{jtype, kport} = [PSTH_Trigger{jtype, kport}; psth_trig_this{1}];
            tspkmat     = psth_trig_this{4};
            trialspxmat = psth_trig_this{3};

            StatOut = ExamineTaskResponsive(tspkmat, trialspxmat);
            StatOut.CellIndx = r.Units.SpikeNotes(i, :);
            PSTH_TriggerStat{jtype, kport}.StatOut(i) =  StatOut;

            % Reward
            PSTH_Reward{jtype, kport} = [PSTH_Reward{jtype, kport}; psth_reward_this{1}];
            tspkmat     = psth_reward_this{4};
            trialspxmat = psth_reward_this{3};

            % restricting activity to a relatively narrow window.
            tRewardWin = [-1000 1000];
            indWin = find(tspkmat>=tRewardWin(1) & tspkmat<=tRewardWin(2));
            tspkmat = tspkmat(indWin);
            trialspxmat = trialspxmat(indWin, :);

            StatOut = ExamineTaskResponsive(tspkmat, trialspxmat);
            StatOut.CellIndx = r.Units.SpikeNotes(i, :);
            PSTH_RewardStat{jtype, kport}.StatOut(i) = StatOut;
            PSTH_baseline = [PSTH_baseline psth_reward_this{1}];

            % compute z
            mean_baseline = mean(PSTH_baseline, 'omitnan');
            sd_baseline   = std(PSTH_baseline, 'omitnan');
            PSTH_CentInZ{jtype, kport}  = [PSTH_CentInZ{jtype, kport}; (psth_centin_this{1}-mean_baseline)/sd_baseline];
            PSTH_CentOutZ{jtype, kport} = [PSTH_CentOutZ{jtype, kport}; (psth_centout_this{1}-mean_baseline)/sd_baseline];
            PSTH_RewardZ{jtype, kport}  = [PSTH_RewardZ{jtype, kport}; (psth_reward_this{1}-mean_baseline)/sd_baseline];

            PSTH_CentInZ{jtype, kport}(isnan(PSTH_CentInZ{jtype, kport}) | isinf(PSTH_CentInZ{jtype, kport}))    = 0;
            PSTH_CentOutZ{jtype, kport}(isnan(PSTH_CentOutZ{jtype, kport}) | isinf(PSTH_CentOutZ{jtype, kport})) = 0;
            PSTH_RewardZ{jtype, kport}(isnan(PSTH_RewardZ{jtype, kport}) | isinf(PSTH_RewardZ{jtype, kport}))    = 0;
        end
    end
end


PopOut.Name       = r.BehaviorClass.Subject;
PopOut.ImplantLateral = r.ImplantLateral;
PopOut.FP         = r.BehaviorClass.TargetFP;
PopOut.CueUncue   = r.BehaviorClass.CueUncue;
PopOut.TrialTypes = TrialTypes;
PopOut.Ports      = r.BehaviorClass.Ports;
PopOut.Trials     = PSTH_Trials;
PopOut.Session    = r.BehaviorClass.Session;
PopOut.Date       = strrep(r.Meta(1).DateTime(1:11), '-','_');
PopOut.Units      = Units;
PopOut.CentIn     = PSTH_CentIn;
PopOut.CentInZ    = PSTH_CentInZ;
PopOut.CentInStat = PSTH_CentInStat;

PopOut.CentInAll     = PSTH_CentInAll;
PopOut.CentInAllZ    = PSTH_CentInAllZ;
PopOut.CentInAllStat = PSTH_CentInAllStat;

PopOut.CentOut     = PSTH_CentOut;
PopOut.CentOutZ    = PSTH_CentOutZ;
PopOut.CentOutStat = PSTH_CentOutStat;

PopOut.CentOutAll     = PSTH_CentOutAll;
PopOut.CentOutAllZ    = PSTH_CentOutAllZ;
PopOut.CentOutAllStat = PSTH_CentOutAllStat;

PopOut.Reward     = PSTH_Reward;
PopOut.RewardZ    = PSTH_RewardZ;
PopOut.RewardStat = PSTH_RewardStat;

PopOut.Trigger     = PSTH_Trigger;
PopOut.TriggerZ    = PSTH_TriggerZ;
PopOut.TriggerStat = PSTH_TriggerStat;

[PopOut.IndSort, PopOut.IndUnmodulated] = SpikesGPS.Timing.RankActivityFlat(PopOut);


PopOut.IndSort = SpikesGPS.Timing.VisualizePopPSTH(PopOut, 'Each');
SpikesGPS.Timing.VisualizePopPSTH(PopOut, 'Ipsi');
SpikesGPS.Timing.VisualizePopPSTH(PopOut, 'Contra');

r.PopPSTH = PopOut;
r_name = 'RTarray_'+r.PSTH.PSTHs(1).ANM_Session{1}+'_'+strrep(r.PSTH.PSTHs(1).ANM_Session{2}, '_', '')+'.mat';
save(r_name, 'r');
% Save a copy of PSTHOut to a collector folder
tosavename = 'PopOut_'+r.PSTH.PSTHs(1).ANM_Session{1}+'_'+strrep(r.PSTH.PSTHs(1).ANM_Session{2}, '_', '')+'.mat';
save(tosavename, 'PopOut');
