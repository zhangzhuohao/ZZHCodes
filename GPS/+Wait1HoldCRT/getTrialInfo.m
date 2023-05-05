function obj = getTrialInfo(obj, SessionData)

obj.FP = zeros(obj.NumTrials, 1);
for i = 1:SessionData.nTrials
    obj.FP(i) = SessionData.RawEvents.Trial{1, i}.Events.BNC1High(1) - SessionData.RawEvents.Trial{1, i}.States.FP(1);
end
obj.FP = round(10 * obj.FP) / 10;
obj.MixedFP = unique(obj.FP);

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

    obj.CentPokeInTime{i} = iStates.FP(:, 1);

    obj.PortCorrect(i) = SessionData.Custom.CorPort(i);

    if ~isnan(iStates.Premature(1)) % Premature
        obj.Outcome{i} = 'Premature';
        obj.CentPokeOutTime{i} = iStates.FP(:, 2);
        obj.ChoiceCueTime(i, :) = [nan, nan]; % when premature, the choice light would not lit up
        obj.TriggerCueTime(i) = nan;
        obj.ChoicePokeTime(i) = nan;
        obj.PortChosen(i) = nan;
    elseif ~isnan(iStates.WrongPort(1)) % selected the wrong port
        obj.Outcome{i} = 'Wrong';
        obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Choice(1)];
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
        obj.CentPokeOutTime{i} = [iStates.FP(1:end-1, 2); iStates.Wait4Choice(1)];
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
    end
end

obj.Stage = obj.FP==1.5;

end