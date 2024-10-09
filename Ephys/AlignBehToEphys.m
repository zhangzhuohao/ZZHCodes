function EventOut = AlignBehToEphys(EventOut, BehClass)

% 3/3/2021 Jianing Yu - AlignMED2BR
% EventOut comes from BlackRock's digital input
% bMED is the b array coming from MED data
% Time of some critical behavioral events (e.g., Trigger stimulus) needs to be mapped to EventOut
% Alignment is performed using press onset data
% Alignment of each trigger stimulus needs to be adjusted to the preceding press
% onset

% 4/2/2023 - AlignBehaviorClassToBR - revised from AlignMED2BR but use behavior class instead
% the goal of this function is to find the behavioral meaning of each blackrock's press 
% Lever presses and releases recorded in blackrock

% 10/9/2024 Zhuohao Zhang - AlignBehToEphys - revised from AlignBehaviorClassToBR for GPS tasks
% the goal of this function is to find the behavioral meaning of each neuropixels' PokeCentIn 
% PokeChoiceIn and PokeInitIn. And complete behavior events in ephys
% recording (Trigger, PokeCentOut)

% these are times for poke-in center, choice and init recorded in ephys.
PokeCentInEphys   = EventOut.Onset{strcmp(EventOut.EventsLabels, 'PokeCentIn')};
% PokeChoiceInEphys = EventOut.Onset{strcmp(EventOut.EventsLabels, 'PokeChoiceIn')};
% PokeInitInEphys   = EventOut.Onset{strcmp(EventOut.EventsLabels, 'PokeInitIn')};

NumTrialsEphys    = length(PokeCentInEphys);

% these are times for poke-in center, choice and init recorded in bpod
Beh = BehClass.BehavTable;

PokeCentInBeh   = Beh.TrialCentInTime*1000;
PokeCentOutBeh  = (Beh.TrialStartTime+Beh.CentOutTime)*1000;
PokeChoiceInBeh = (Beh.TrialStartTime+Beh.ChoicePokeTime)*1000;
PokeInitInBeh   = (Beh.TrialStartTime+Beh.InitInTime)*1000;
TriggerBeh      = (Beh.TrialStartTime+Beh.TriggerCueTime)*1000;

PerformanceBeh  = Beh.Outcome;
TrialIndexBeh   = Beh.Trials;
PortCorrectBeh  = Beh.PortCorrect;
FPBeh           = Beh.FP;
CueBeh         = Beh.Cued;
StageBeh        = Beh.Stage;

%
% Start to map
% PokeCentInBeh is the time recorded in Bpod
% PokeCentInEphys is the time recorded in Ephys system
% IndMatched gives the matching index. 

IndMatched = findseqmatchrev(PokeCentInBeh, PokeCentInEphys, 0, 1); % note that IndMatched has the same length as PokeCentInEphys
IndNaN     = find(isnan(IndMatched));
IndValid   = find(~isnan(IndMatched));

%
TrialIndexEphys                = 1:NumTrialsEphys; % this is 1 to the total number of poke-cent recorded in ephys
TrialIndexEphys(IndNaN)        = NaN; % some ephys-poke-cent cannot be mapped
TrialIndexEphys(IndValid)      = TrialIndexBeh(IndMatched(IndValid)); % only for these presses we can find matching ones in behavior

PerformanceEphys               = strings(1, NumTrialsEphys); % each trial corresponds to a behavioral outcome
PerformanceEphys(IndNaN)       = "";
PerformanceEphys(IndValid)     = PerformanceBeh(IndMatched(IndValid));

CueEphys                       = ones(1, NumTrialsEphys); % each press corresponds to a behavioral outcome
CueEphys(IndNaN)               = NaN;
CueEphys(IndValid)             = CueBeh(IndMatched(IndValid));

FPEphys                        = zeros(1, NumTrialsEphys); % each press corresponds to a behavioral outcome
FPEphys(IndNaN)                = NaN;
FPEphys(IndValid)              = FPBeh(IndMatched(IndValid));

PortCorrectEphys               = zeros(1, NumTrialsEphys); % each press corresponds to a behavioral outcome
PortCorrectEphys(IndNaN)       = NaN;
PortCorrectEphys(IndValid)     = PortCorrectBeh(IndMatched(IndValid));

StageEphys                     = zeros(1, NumTrialsEphys); % each press corresponds to a behavioral outcome
StageEphys(IndNaN)             = NaN;
StageEphys(IndValid)           = StageBeh(IndMatched(IndValid));

% add information to EventOut
EventOut.TrialIndexEphys       = TrialIndexEphys;
EventOut.OutcomeEphys          = PerformanceEphys;
EventOut.CueEphys              = CueEphys;
EventOut.FPEphys               = FPEphys;
EventOut.PortCorrectEphys      = PortCorrectEphys;
EventOut.StageEphys            = StageEphys;

% there are no marks for trigger and poke-cent-out during GPS neuropixels
% recording, we can fill them in through poke-cent-in time
TriggerMapped = nan(1, NumTrialsEphys);
for i = 1:NumTrialsEphys
    iPokeCentInEphys = PokeCentInEphys(i);
    if ~isnan(TrialIndexEphys(i))
        i_trial = TrialIndexEphys(i);
        iPokeCentInBeh = PokeCentInBeh(i_trial);
        iTriggerBeh = TriggerBeh(i_trial);
        if i <= NumTrialsEphys
            if TriggerBeh(i_trial) > 0
                TriggerMapped(i) = iTriggerBeh-iPokeCentInBeh+iPokeCentInEphys;
            end
        end
    end
end
TriggerMapped = TriggerMapped(~isnan(TriggerMapped));

PokeCentOutMapped = nan(1, NumTrialsEphys);
for i = 1:NumTrialsEphys
    iPokeCentInEphys = PokeCentInEphys(i);
    if ~isnan(TrialIndexEphys(i))
        i_trial = TrialIndexEphys(i);
        iPokeCentInBeh  = PokeCentInBeh(i_trial);
        iPokeCentOutBeh = PokeCentOutBeh(i_trial);
        if ~isnan(iPokeCentInBeh)
            PokeCentOutMapped(i) = iPokeCentOutBeh-iPokeCentInBeh+iPokeCentInEphys;
        end
    end
end

% find out the min distance between triggertimemapped and the ones recorded
% in blackrock
figure;
subplot(2, 1, 1)
plot(PokeCentInEphys, 5, 'ko');
hold on
line([PokeCentOutMapped' PokeCentOutMapped'], [4 6]', 'color', 'm')
set(gca, 'ylim', [4 6])
ylabel('Poke-out center')

subplot(2, 1, 2)
plot(PokeCentInEphys, 5, 'ko');
hold on
line([TriggerMapped' TriggerMapped'], [4 6]', 'color', 'm')
set(gca, 'ylim', [4 6])
ylabel('Trigger')

% add event time info to EventOut
if isempty(find(strcmp(EventOut.EventsLabels, 'Trigger'), 1))
    EventOut.EventsLabels{end+1}='Trigger';
end
TriggerDur = 250; % 250 ms trigger events
EventOut.Onset{strcmp(EventOut.EventsLabels, 'Trigger')}  = TriggerMapped;
EventOut.Offset{strcmp(EventOut.EventsLabels, 'Trigger')} = TriggerMapped + TriggerDur;

if isempty(find(strcmp(EventOut.EventsLabels, 'PokeCentOut'), 1))
    EventOut.EventsLabels{end+1}='PokeCentOut';
end
EventOut.Onset{strcmp(EventOut.EventsLabels, 'PokeCentOut')}  = PokeCentOutMapped;
EventOut.Offset{strcmp(EventOut.EventsLabels, 'PokeCentOut')} = PokeCentOutMapped + 200;

end % AlignBehToEphys