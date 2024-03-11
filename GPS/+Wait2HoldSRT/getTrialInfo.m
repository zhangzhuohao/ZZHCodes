function obj = getTrialInfo(obj, SessionData)
% Get trial information of events' time points

obj.FP = zeros(obj.NumTrials, 1);
for i = 1:SessionData.nTrials
    obj.FP(i) = SessionData.RawEvents.Trial{1, i}.Events.BNC1High(1) - SessionData.RawEvents.Trial{1, i}.States.FP(1);
end
obj.FP = round(10 * obj.FP) * 0.1;
obj.MixedFP = unique(obj.FP);

obj.RW = SessionData.Custom.RW';

obj.InitPokeInTime = cell(obj.NumTrials, 1);
obj.InitPokeOutTime = cell(obj.NumTrials, 1);

obj.CentPokeInTime = cell(obj.NumTrials, 1);
obj.CentPokeOutTime = cell(obj.NumTrials, 1);

obj.ChoiceCueTime = zeros(obj.NumTrials, 2);
obj.TriggerCueTime = zeros(obj.NumTrials, 1);
obj.ChoicePokeTime = zeros(obj.NumTrials, 1);

obj.PortCorrect = zeros(obj.NumTrials, 1);
obj.PortChosen = zeros(obj.NumTrials, 1);
obj.Outcome = cell(obj.NumTrials, 1);

for i = 1:obj.NumTrials
    iStates = SessionData.RawEvents.Trial{i}.States;
    iEvents = SessionData.RawEvents.Trial{i}.Events;

    obj.InitPokeInTime{i} = iEvents.Port3In(iEvents.Port3In>=iStates.Wait4Sample(1) & iEvents.Port3In<=iStates.Wait4Center(2));
    obj.InitPokeOutTime{i} = iEvents.Port3Out(iEvents.Port3Out>=iStates.Wait4Sample(1) & iEvents.Port3Out<=iStates.Wait4Center(2));

    if isfield(iStates, 'Wait4Out')
        obj.CentPokeInTime{i} = [iStates.FP(:, 1); iStates.Wait4Out(2:end, 1); iStates.Late(2:end, 1)];
    else
        obj.CentPokeInTime{i} = iStates.FP(:, 1);
    end

    obj.PortCorrect(i) = SessionData.Custom.CorPort(i);

    if ~isnan(iStates.Premature(1)) % Premature
        obj.Outcome{i} = 'Premature';
        obj.CentPokeOutTime{i} = iStates.FP(:, 2);
        obj.ChoiceCueTime(i, :) = [nan, nan]; % when premature, the choice light would not lit up
        obj.TriggerCueTime(i) = nan;
        obj.ChoicePokeTime(i) = nan;
        obj.PortChosen(i) = nan;

        if iStates.FP(end, 2) - iStates.FP(1, 1) > obj.FP(i)
            obj.Outcome{i} = 'Bug';
        end
    elseif ~isnan(iStates.Late(1)) % Late
        if isfield(iStates, 'Wait4Out')
            obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Out(1:end-1, 2); iStates.Late(:, 2)];
        else
            obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Late(2)];
        end
        obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Late(2)];
        obj.TriggerCueTime(i) = iStates.ChoiceCue(1);

        if ~any(isfield(iStates, ["LateWrong", "LateCorrect"]))
            obj.Outcome{i} = 'Late';
            obj.ChoicePokeTime(i) = nan;
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            obj.PortChosen(i) = nan;
        elseif ~isnan(iStates.LateWrong(1))
            obj.Outcome{i} = 'LateWrong';
            obj.ChoicePokeTime(i) = iStates.LateWrong(1);
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            if isfield(iEvents, 'Port2In') && sum(ismember(iEvents.Port2In, iStates.LateWrong(1)))
                obj.PortChosen(i) = 2;
            elseif  isfield(iEvents, 'Port1In') && sum(ismember(iEvents.Port1In, iStates.LateWrong(1)))
                obj.PortChosen(i) = 1;
            else % sometimes, they poke during the choice light presentation
                obj.PortChosen(i) = NaN;
            end
        elseif ~isnan(iStates.LateCorrect(1))
            obj.Outcome{i} = 'LateCorrect';
            obj.ChoicePokeTime(i) = iStates.LateCorrect(1);
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            if isfield(iEvents, 'Port2In') && sum(ismember(iEvents.Port2In, iStates.LateCorrect(1)))
                obj.PortChosen(i) = 2;
            elseif  isfield(iEvents, 'Port1In') && sum(ismember(iEvents.Port1In, iStates.LateCorrect(1)))
                obj.PortChosen(i) = 1;
            else % sometimes, they poke during the choice light presentation
                obj.PortChosen(i) = NaN;
            end
        else
            obj.Outcome{i} = 'LateMiss';
            obj.ChoicePokeTime(i) = nan;
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            obj.PortChosen(i) = nan;
        end
    elseif ~isnan(iStates.WrongPort(1)) % selected the wrong port
        obj.Outcome{i} = 'Wrong';
        if isfield(iStates, 'Wait4Out')
            if ~isnan(iStates.Wait4Out(2))
                obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Out(:, 2)];
            else
                obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.ChoiceCue(2)];
            end
        else
            obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Choice(1)];
        end
        obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Wait4Choice(1)]; % duration that the choice light lit up.
        obj.TriggerCueTime(i) = iStates.ChoiceCue(1);
        obj.ChoicePokeTime(i) = iStates.Wait4Choice(2);
        % figure out the port situation: WrongPort(1)
        % should match a port entry time
        if isfield(iEvents, 'Port2In') && sum(ismember(iEvents.Port2In, iStates.WrongPort(1)))
            obj.PortChosen(i) = 2;
        elseif  isfield(iEvents, 'Port1In') && sum(ismember(iEvents.Port1In, iStates.WrongPort(1)))
            obj.PortChosen(i) = 1;
        else % sometimes, they poke during the choice light presentation
            obj.PortChosen(i) = NaN;
        end
    elseif ~isnan(iStates.Wait4Reward(1))
        obj.Outcome{i} = 'Correct';
        if isfield(iStates, 'Wait4Out')
            if ~isnan(iStates.Wait4Out(2))
                obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Out(:, 2)];
            else
                obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.ChoiceCue(2)];
            end
        else
            obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Choice(1)];
        end
        obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Wait4Choice(1)]; % duration that the choice light lit up.
        obj.TriggerCueTime(i) = iStates.ChoiceCue(1);
        obj.ChoicePokeTime(i) = iStates.Wait4Choice(2);
        % find out the port
        if isfield(iEvents, 'Port2In') && sum(ismember(iEvents.Port2In, iStates.Wait4Reward(1)))
            obj.PortChosen(i) = 2;
        elseif  isfield(iEvents, 'Port1In') && sum(ismember(iEvents.Port1In, iStates.Wait4Reward(1)))
            obj.PortChosen(i) = 1;
        else % sometimes, they poke during the choice light presentation
            obj.PortChosen(i) = NaN;
        end
    else
        obj.Outcome{i} = 'Bug';
    end
end

obj.FP = roundn(obj.FP, -1);

%% remove bug trials
ind_bug = strcmp(obj.Outcome, 'Bug') | obj.FP <= 0;

obj.NumTrials = obj.NumTrials - sum(ind_bug);

obj.Trials(ind_bug)             = [];
obj.TrialStartTime(ind_bug)     = [];

obj.FP(ind_bug)                 = [];
obj.RW(ind_bug)                 = [];

obj.InitPokeInTime(ind_bug)     = [];
obj.InitPokeOutTime(ind_bug)    = [];

obj.CentPokeInTime(ind_bug)     = [];
obj.CentPokeOutTime(ind_bug)    = [];

obj.ChoiceCueTime(ind_bug, :)   = [];
obj.TriggerCueTime(ind_bug)     = [];
obj.ChoicePokeTime(ind_bug)     = [];

obj.PortCorrect(ind_bug)        = [];
obj.PortChosen(ind_bug)         = [];
obj.Outcome(ind_bug)            = [];

end