 function obj = getTrialInfo(obj, SessionData)
% WaitHold task
% Get trial information of events' time points

if contains(obj.Protocol, 'SRT')
    paradigm = "SRT";
elseif contains(obj.Protocol, 'CRT')
    paradigm = "CRT";
elseif contains(obj.Protocol, 'Hold')
    paradigm = "Hold";
elseif contains(obj.Protocol, 'Free')
    paradigm = "Free";
end

%
obj.NumTrials = SessionData.nTrials;

obj.Trials = (1:obj.NumTrials)';
obj.TrialStartTime = SessionData.TrialStartTimestamp(obj.Trials)';

obj.InitPokeInTime = cell(obj.NumTrials, 1);
obj.InitPokeOutTime = cell(obj.NumTrials, 1);

obj.CentPokeInTime = cell(obj.NumTrials, 1);
obj.CentPokeOutTime = cell(obj.NumTrials, 1);

obj.ChoiceCueTime = zeros(obj.NumTrials, 2);
obj.TriggerCueTime = zeros(obj.NumTrials, 1);
obj.ChoicePokeTime = zeros(obj.NumTrials, 1);

obj.PortCorrect = zeros(obj.NumTrials, 1);
obj.PortChosen = zeros(obj.NumTrials, 1);
obj.Outcome = strings(obj.NumTrials, 1);
obj.LateChoice = strings(obj.NumTrials, 1);

if ~isfield(SessionData.Custom, "Version")
    obj.FP = SessionData.Custom.FP';
else
    switch SessionData.Custom.Version
        case {'20230312', '20230618'}
            obj.FP = zeros(obj.NumTrials, 1);
            for i = 1:obj.NumTrials
                if ~isfield(SessionData.RawEvents.Trial{1, i}.Events, 'BNC1High')
                    obj.Outcome(i) = "Bug";
                else
                    obj.FP(i) = SessionData.RawEvents.Trial{1, i}.Events.BNC1High(1) - SessionData.RawEvents.Trial{1, i}.States.FP(1);
                end
            end
    end
end
obj.FP = roundn(obj.FP, -1);
obj.TargetFP = roundn((1:15) / 10, -1);

if ~isfield(SessionData.Custom, 'RW')
    obj.RW = nan(obj.NumTrials, 1);
else
    obj.RW = SessionData.Custom.RW';
end
obj.Cued = ones(obj.NumTrials, 1);

for i = 1:obj.NumTrials

    if obj.Outcome(i)=="Bug"
        continue;
    end

    iStates = SessionData.RawEvents.Trial{i}.States;
    iEvents = SessionData.RawEvents.Trial{i}.Events;

    obj.InitPokeInTime{i} = iEvents.Port3In(iEvents.Port3In>=iStates.Wait4Sample(1) & iEvents.Port3In<=iStates.Wait4Center(2));
    obj.InitPokeOutTime{i} = iEvents.Port3Out(iEvents.Port3Out>=iStates.Wait4Sample(1) & iEvents.Port3Out<=iStates.Wait4Center(2));

    obj.PortCorrect(i) = SessionData.Custom.CorPort(i);

    if ~isnan(iStates.PrematureCorrect(1)) || ~isnan(iStates.PrematureWrong(1)) % Premature
        obj.Outcome(i) = "Premature";
        obj.CentPokeInTime{i} = [iStates.FP(1) iEvents.Port4In(iEvents.Port4In>=iStates.FP(1) & iEvents.Port4In<=iStates.FP(2))];
        obj.CentPokeOutTime{i} = iEvents.Port4Out(iEvents.Port4Out>=iStates.FP(1) & iEvents.Port4Out<=iStates.FP(2));
        switch paradigm
            case {'CRT', 'Hold'}
                obj.ChoiceCueTime(i, :) = [nan, nan]; % when premature, the choice light would not lit up
            case {'SRT'}
                obj.ChoiceCueTime(i, :) = [iStates.FP(1) iStates.FP(2)];
        end
        obj.TriggerCueTime(i) = nan;
        if ~isnan(iStates.PrematureCorrect(1))
            obj.ChoicePokeTime(i) = iStates.PrematureCorrect(1);
            obj.PortChosen(i) = obj.PortCorrect(i);
        else
            obj.ChoicePokeTime(i) = iStates.PrematureWrong(1);
            obj.PortChosen(i) = 3 - obj.PortCorrect(i);
        end
        if iStates.FP(end, 2) - iStates.FP(1, 1) > obj.FP(i)
            obj.Outcome(i) = "Bug";
        end
    elseif isfield(iStates, 'Late') && ~isnan(iStates.Late(1)) % Late
        obj.Outcome(i) = "Late";
        obj.CentPokeInTime{i} = [iStates.FP(1) iEvents.Port4In(iEvents.Port4In>=iStates.FP(1) & iEvents.Port4In<=iStates.Late(end, 2))];
        obj.CentPokeOutTime{i} = iEvents.Port4Out(iEvents.Port4Out>=iStates.FP(1) & iEvents.Port4Out<=iStates.Late(end, 2));
        switch paradigm
            case {'CRT', 'Hold'}
                obj.ChoiceCueTime(i, :) = [iStates.FP(2), iStates.Late(end, 2)];
            case {'SRT'}
                obj.ChoiceCueTime(i, :) = [iStates.FP(1), iStates.Late(end, 2)];
        end
        obj.TriggerCueTime(i) = iStates.ChoiceCue(1);

        if ~any(isfield(iStates, ["LateWrong", "LateCorrect"]))
            obj.ChoicePokeTime(i) = nan;
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            obj.PortChosen(i) = nan;
        elseif ~isnan(iStates.LateWrong(1))
            obj.LateChoice(i) = "Wrong";
            obj.ChoicePokeTime(i) = iStates.LateWrong(1);
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            if isfield(iEvents, 'Port2In') && sum(ismember(iEvents.Port2In, iStates.LateWrong(1)))
                obj.PortChosen(i) = 2;
            elseif  isfield(iEvents, 'Port1In') && sum(ismember(iEvents.Port1In, iStates.LateWrong(1)))
                obj.PortChosen(i) = 1;
            else % sometimes, they poke during the choice light presentation
                obj.PortChosen(i) = nan;
            end
        elseif ~isnan(iStates.LateCorrect(1))
            obj.LateChoice(i) = "Correct";
            obj.ChoicePokeTime(i) = iStates.LateCorrect(1);
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            if isfield(iEvents, 'Port2In') && sum(ismember(iEvents.Port2In, iStates.LateCorrect(1)))
                obj.PortChosen(i) = 2;
            elseif  isfield(iEvents, 'Port1In') && sum(ismember(iEvents.Port1In, iStates.LateCorrect(1)))
                obj.PortChosen(i) = nan;
            else % sometimes, they poke during the choice light presentation
                obj.PortChosen(i) = nan;
            end
        else
            obj.LateChoice(i) = "Miss";
            obj.ChoicePokeTime(i) = nan;
            % figure out the port situation: WrongPort(1)
            % should match a port entry time
            obj.PortChosen(i) = nan;
        end
    elseif ~isnan(iStates.WrongPort(1)) % selected the wrong port
        obj.Outcome(i) = "Wrong";
        obj.CentPokeInTime{i} = [iStates.FP(1) iEvents.Port4In(iEvents.Port4In>=iStates.FP(1) & iEvents.Port4In<=iStates.WrongPort(1))];
        obj.CentPokeOutTime{i} = iEvents.Port4Out(iEvents.Port4Out>=iStates.FP(1) & iEvents.Port4Out<=iStates.WrongPort(1));
        switch paradigm
            case {'CRT'}
                obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Wait4Choice(1)];
            case {'Hold', 'Free'}
                obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Wait4Choice(2)]; % duration that the choice light lit up.
            case {'SRT'}
                obj.ChoiceCueTime(i, :) = [iStates.FP(1) iStates.Wait4Choice(2)];
        end
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
        obj.Outcome(i) = 'Correct';
        obj.CentPokeInTime{i} = [iStates.FP(1) iEvents.Port4In(iEvents.Port4In>iStates.FP(1) & iEvents.Port4In<=iStates.Wait4Reward(1))];
        obj.CentPokeOutTime{i} = iEvents.Port4Out(iEvents.Port4Out>=iStates.FP(1) & iEvents.Port4Out<=iStates.Wait4Reward(1));
        switch paradigm
            case {'CRT'}
                obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Wait4Choice(1)];
            case {'Hold', 'Free'}
                obj.ChoiceCueTime(i, :) = [iStates.ChoiceCue(1) iStates.Wait4Choice(2)]; % duration that the choice light lit up.
            case {'SRT'}
                obj.ChoiceCueTime(i, :) = [iStates.FP(1) iStates.Wait4Choice(2)];
        end        
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
        obj.Outcome(i) = "Bug";
    end
end

obj.FP = roundn(obj.FP, -1);

%% remove bug trials
ind_bug = strcmp(obj.Outcome, "Bug") | obj.FP <= 0;

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
obj.LateChoice(ind_bug)         = [];

obj.Cued(ind_bug)               = [];
end