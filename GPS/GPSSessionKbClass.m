classdef GPSSessionKbClass
    % Turn SessionData file into a class
    %   Detailed explanation goes here
    %

    properties
        % Experiment information
        BpodFile
        BpodFilePath
        BpodFileName
        Subject % Name of the animal, extrated from data file name

        ANMInfoFile
        
        Task % Name of task, extrated from data file name and set as {'Autoshaping', 'Wait1Hold', 'Wait1HoldCRT', 'Wait2HoldCRT', 'ThreeFPHoldCRT', 'ThreeFPHoldSRT', 'Kornblum'}
        Session % Date of this session, extrated from data file name

        % Session information
        SessionStartTime % Start time of this session, in the form of hh:mm:ss, extraced from SessionData.Info
        NumTrials % Total trial number, extraced from SessionData.Info
        Trials % Trial index, [NumTrials x 1]
        TrialStartTime % Time when trial started (last trial ended) relative to SessionStartTime, [NumTrials x 1]

        % Trial information
        InitPokeInTime % Time when the animal first poked-in port_init, [NumTrials x 1]
        InitPokeOutTime % Time when the animal poked-out port_init, and went to port_cent, [NumTrials x 1]
        CentPokeInTime % Time when the animal poked-in port_cent, [NumTrials x 1]
        CentPokeOutTime % Time when the animal poked-out port_cent, and went to make choice, [NumTrials x 1]
        ChoicePokeTime % Time when the animal poked-in one of the port_choice (port_left/port_right), [NumTrials x 1]
        ChoiceCueTime % Duration that choice light lit up, [NumTrials x 2]
        TriggerCueTime

        PortCorrect % The correct for each trial (1 for port_left, 2 for port_right), [NumTrials x 1]
        PortChosen % The port chosen by the animal in each trial, [NumTrials x 1]
        Outcome % Outcome of each trial, (1 for correct; 0 for choice-error; -1 for premature; -2 for late), [NumTrials x 1]

        MixedFP % FP configuration, for ThreeFP paradigm only, [NumTrials x N (number of FPs)]
        FP % foreperiod in each trial, for all paradigms except Autoshaping, [NumTrials x 1]
        RW

        Cued
    end

    properties (Constant)
        PerformanceType = ["Correct", "Premature", "Late", "Wrong"];
        PerformanceCode = [1, -1, -2, 0];
        Definition = {
            'ShuttleTime: Time of moving from port_init to port_cent';
            'ReactionTime: For hold paradigm, it is time from tone to port_cent poke-out; for free paradigm, it is the time from tone to port_choice poke-in';
            'MovementTime: For hold paradigm, it is the time from port_cent poke-out to port_choice poke-in; for free paradigm, leave it empy; for autoshaping paradigm, it is the time from port_cent poke-in to port_choice poke-in';
            'HoldDuration: Time from port_cent poke out to port_choice poke in (for hold paradigm only)';
            };
        Ports = ["L", "R"];
        CueUncue = [1 0];
    end

    properties (Dependent)
        Strain % Strain of the animal, set manually
        Gender
        Experimenter % Name (initials) of the Experimenter, set manually
        Treatment % Set manually
        Dose % Set manually
        Label
        RTSortedLabels

        Stage % warm-up or 3FP
        Ind
        Bins

        ShuttleTime % check
        ShuttleTimeDistribution % check
        ShuttleTimeStat

        MovementTime
        MovementTimeSorted
        MovementTimePDF
        MovementTimeCDF
        MovementTimeStat

        ChoiceTime
        ChoiceTimeSorted
        ChoiceTimePDF
        ChoiceTimeCDF
        ChoiceTimeStat

        RT % check
        RTSorted
        RTPDF
        RTCDF
        RTStat % check

        HoldDuration %
        HoldDurationSorted
        HoldDurationPDF
        HoldDurationCDF
        HoldDurationStat

        Interruption

        Performance % Performance table, grouped by PortCorrect (and FP)
        PerformanceTrackCue
        PerformanceTrackUncue

        BehavTable
    end

    methods
        %% Initiate
        function obj = GPSSessionKbClass(BpodFile, AnmInfoFile)
            % Process SessionData to get properties
            load(BpodFile, 'SessionData');
            obj.BpodFile = BpodFile;
            [obj.BpodFilePath, obj.BpodFileName] = fileparts(obj.BpodFile);

            obj.ANMInfoFile = AnmInfoFile;

            % Rat name
            obj.Subject = extractBefore(obj.BpodFileName, '_');
            Protocol = extractAfter(obj.BpodFileName, [obj.Subject '_']);
            switch Protocol(1:end-16) % e.g. 'NEW_03_Wait3FPFlash'
                case {'GPS_08_KornblumHold'}
                    obj.Task = 'KornblumSRT';
                case {'GPS_09_KornblumHoldEmpty'}
                    obj.Task = 'KornblumSRTEmp';
            end

            % Session meta-information
            obj.Session = Protocol(end-14:end-7);
            obj.SessionStartTime = SessionData.Info.SessionStartTime_UTC;
            obj.NumTrials = SessionData.nTrials;
            if isfield(SessionData.RawEvents.Trial{end}.Events, 'WavePlayer1_2') % Sometime the protocol might shut down
                obj.NumTrials = obj.NumTrials - 1;
            end
            obj.Trials = (1:obj.NumTrials)';
            obj.TrialStartTime = SessionData.TrialStartTimestamp(obj.Trials)';

            % Get paradigm specific trial information
            obj = feval([obj.Task '.getTrialInfo'], obj, SessionData);
        end

        %% Manually set
        function obj = set.Experimenter(obj, person_name)
            % Manually set the name (initials) of the experimenter, and
            % turn it to string
            obj.Experimenter = string(person_name);
        end

        function obj = set.Treatment(obj,treatment)
            % Manually set the treatment of this session, NaN/Saline/DCZ
            if ismember(treatment, {'None', 'Saline', 'DCZ'})
                obj.Treatment = string(treatment);
            else
                error('Treatment can only be: None, Saline, DCZ')
            end
        end

        function obj = set.Dose(obj,dose)
            % Manually set the dose of injection
            if isnumeric(dose)
                obj.Dose = dose;
            else
                error('"dose" must be a scalar')
            end
        end

        function obj = set.Label(obj,label)
            % Manually set the treatment of this session, NaN/Saline/DCZ
            if ismember(label, {'None', 'Control', 'Chemo'})
                obj.Treatment = string(label);
            else
                error('Treatment can only be: None, Control, Chemo')
            end
        end

        function obj = set.MixedFP(obj,  mFP)
            if isnumeric(mFP)
                if mean(mFP)>100
                    error('Please use seconds')
                else
                    obj.MixedFP = mFP;
                end
            end
        end

        function obj = set.Strain(obj,strain)
            if ismember(strain, {'BN', 'Wistar', 'LE', 'SD', 'Hybrid'})
                obj.Strain = string(strain);
            else
                error("Strain can only be: 'BN', 'Wistar', 'LE', 'SD', 'Hybrid'")
            end
        end

        function obj = set.Gender(obj,strain)
            if ismember(strain, {'F', 'M'})
                obj.Strain = string(strain);
            else
                error("Strain can only be: 'F', 'M'")
            end
        end

        %% Get animal information
        function value = get.Strain(obj)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", "ANM", "TextType", "string");
            strain = anm_info.Strain(strcmp(anm_info.Name, obj.Subject));

            value = strain;
        end

        function value = get.Gender(obj)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", "ANM", "TextType", "string");
            gender = anm_info.Gender(strcmp(anm_info.Name, obj.Subject));

            value = gender;
        end

        function value = get.Treatment(obj)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");
            treatment = anm_info.Treatment(strcmp(string(anm_info.Session), obj.Session));

            value = treatment;
        end

        function value = get.Dose(obj)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");
            dose = anm_info.Dose(strcmp(string(anm_info.Session), obj.Session));

            value = dose;
        end

        function value = get.Label(obj)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");
            label = anm_info.Label(strcmp(string(anm_info.Session), obj.Session));

            value = label;
        end

        function value = get.Experimenter(obj)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");
            experimenter = anm_info.Experimenter(strcmp(string(anm_info.Session), obj.Session));

            value = experimenter;
        end

        function value = get.RTSortedLabels(obj)
            switch obj.Task
                case {'KornblumSRT'}
                    value = 'Cued x Ports';
                otherwise
                    value = 'FPs x Ports';
            end
        end

        %%
        function value = get.Stage(obj)
            stage = feval([obj.Task '.getStage'], obj);
            value = stage;
        end

        function value = get.Ind(obj)
            ind = feval([obj.Task '.getInd'], obj);
            value = ind;
        end

        function value = get.Bins(obj)
            bins = feval([obj.Task '.getBins']);
            value = bins;
        end

        %% Time interval info (same way to get for each paradigm)
        %% Shuttle time
        function value = get.ShuttleTime(obj)
            % find the last poke out time before center poke
            shuttle_time = cellfun(@(x, y) y(1)-x(end), obj.InitPokeOutTime, obj.CentPokeInTime);
            value = shuttle_time;
        end

        % Shuttle time statistics
        function value = get.ShuttleTimeStat(obj)
            [dataout, interq] = rmoutliers_custome(obj.ShuttleTime);

            ST_Mean = mean(dataout);
            ST_Median = median(dataout);
            ST_IQR = interq;
            ST_STD = std(dataout);

            Subjects = [string(obj.Subject)];
            Sessions = [string(obj.Session)];

            StatVars = {'Subjects', 'Sessions', 'Mean (s)',  'Median (s)', 'IQR (s)', 'STD (s)'};
            thisTable = table(Subjects, Sessions, ST_Mean, ST_Median, ST_IQR, ST_STD, 'VariableNames', StatVars);
            value = thisTable;
        end

        function value = get.ShuttleTimeDistribution(obj)
            datain = obj.ShuttleTime;
            dataout = rmoutliers_custome(datain);

            % make it logarithmic
            dataout_log = log(dataout);
            binEdges = obj.Bins.ShuttleTimeLog;
            binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;
            Ncounts = histcounts(dataout_log, binEdges);
            Counts = Ncounts';
            ST_Centers = exp(binCenters)';
            value = table(ST_Centers, Counts);
        end

        %% Reaction time
        function value = get.RT(obj)
            % find the last poke out time before center poke
            reaction_time = cellfun(@(x) x(end), obj.CentPokeOutTime) - obj.TriggerCueTime;
            value = reaction_time;
        end

        function value = get.RTSorted(obj)
            value = obj.sortData("RT", "correct");
        end

        function value = get.RTPDF(obj)
            value = obj.getPDF("RT", "correct");
        end

        function value = get.RTCDF(obj)
            value = obj.getCDF("RT", "correct");
        end

        function value = get.RTStat(obj)
            value = obj.getStat("RT", "correct");
        end

        %% Hold duration
        function value = get.HoldDuration(obj)
            % find the last poke out time before center poke
            hold_duration = cellfun(@(x, y) y(end) - x(1), obj.CentPokeInTime, obj.CentPokeOutTime);
            value = hold_duration;
        end

        function value = get.HoldDurationSorted(obj)
            value = obj.sortData("HoldDuration", "port");
        end

        function value = get.HoldDurationPDF(obj)
            value = obj.getPDF("HoldDuration", "port");
        end

        function value = get.HoldDurationCDF(obj)
            value = obj.getCDF("HoldDuration", "port");
        end

        function value = get.HoldDurationStat(obj)
            value = obj.getStat("HoldDuration", "port");
        end

        %% Movement time
        function value = get.MovementTime(obj)
            movement_time = obj.ChoicePokeTime - cellfun(@(x) x(end), obj.CentPokeOutTime);
            value = movement_time;
        end

        function value = get.MovementTimeSorted(obj)
            value = obj.sortData("MovementTime", "correct");
        end

        function value = get.MovementTimePDF(obj)
            value = obj.getPDF("MovementTime", "correct");
        end

        function value = get.MovementTimeCDF(obj)
            value = obj.getCDF("MovementTime", "correct");
        end

        function value = get.MovementTimeStat(obj)
            value = obj.getStat("MovementTime", "correct");
        end

        %% Choice time
        function value = get.ChoiceTime(obj)
            choice_time = obj.RT + obj.MovementTime;
            value = choice_time;
        end

        function value = get.ChoiceTimeSorted(obj)
            value = obj.sortData("ChoiceTime", "correct");
        end

        function value = get.ChoiceTimePDF(obj)
            value = obj.getPDF("ChoiceTime", "correct");
        end

        function value = get.ChoiceTimeCDF(obj)
            value = obj.getCDF("ChoiceTime", "correct");
        end

        function value = get.ChoiceTimeStat(obj)
            value = obj.getStat("ChoiceTime", "correct");
        end

        %%
        function value = get.Interruption(obj)

            inter_on  = cellfun(@(x, y) x(1:end-1) - y(1), obj.CentPokeOutTime, obj.CentPokeInTime, 'UniformOutput', false);
            inter_dur = cellfun(@(x, y) y(2:end) - x(1:end-1), obj.CentPokeOutTime, obj.CentPokeInTime, 'UniformOutput', false);

            inter_on_collect  = [];
            inter_dur_collect = [];
            trial             = [];
            port_correct      = [];
            fp                = [];
            cued              = [];
            hold_duration     = [];
            outcome           = [];

            valid_trials = find(obj.Stage);
            for t = 1:length(valid_trials)
                inter_on_collect  = [inter_on_collect  ; inter_on{t}];
                inter_dur_collect = [inter_dur_collect ; inter_dur{t}];
                trial             = [trial             ; repmat(t, length(inter_on{t}), 1)];
                port_correct      = [port_correct      ; repmat(obj.Ports(obj.PortCorrect(t)), length(inter_on{t}), 1)];
                fp                = [fp                ; repmat(obj.FP(t), length(inter_on{t}), 1)];
                cued              = [cued              ; repmat(obj.Cued(t), length(inter_on{t}), 1)];
                hold_duration     = [hold_duration     ; repmat(obj.HoldDuration(t), length(inter_on{t}), 1)];
                outcome           = [outcome           ; repmat(obj.Outcome(t), length(inter_on{t}), 1)];
            end

            if isempty(trial)
                inter_on_collect  = zeros(1,1);
                inter_dur_collect = zeros(1,1);
                trial             = zeros(1,1);
                port_correct      = zeros(1,1);
                fp                = zeros(1,1);
                cued              = zeros(1,1);
                hold_duration     = zeros(1,1);
                outcome           = zeros(1,1);
            end

            Subjects = repmat(string(obj.Subject), length(trial), 1);
            Sessions = repmat(string(obj.Session), length(trial), 1);

            inter = table(Subjects, Sessions, trial, inter_on_collect, inter_dur_collect, port_correct, fp, cued, hold_duration, outcome, ...
                'VariableNames', {'Subjects', 'Sessions', 'Trials', 'On', 'Dur', 'PortCorrect', 'FP', 'Cued', 'HoldDuration', 'Outcome'});
            inter(inter.Dur<0.001, :) = [];
            
            value = inter;
        end

        %%
        function value = get.Performance(obj)
            cueuncue = obj.CueUncue;
            if size(cueuncue, 1)==1
                cueuncue = cueuncue';
            end

            TargetPort = repmat(["L"; "R"; "Both"], length(obj.CueUncue), 1);
            Cued_this = [repmat(cueuncue(1), 3, 1); repmat(cueuncue(2), 3, 1)];
            NumTrialsSorted = zeros(length(Cued_this), 1);

            % group by foreperiod
            CorrectRatio = zeros(length(Cued_this), 1);
            PrematureRatio = zeros(length(Cued_this), 1);
            LateRatio = zeros(length(Cued_this), 1);
            WrongRatio = zeros(length(Cued_this), 1);

            % group by target port
            format short

            for i = 1:length(obj.CueUncue)
                for j = 1:3 % Target port
                    ind = j + 3*(i-1);
                    if mod(j, length(obj.CueUncue)+1)~=0 % group by FP and Target port
                        n_correct   = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & eval("obj.Ind.correct"   + obj.Ports(j)) );
                        n_premature = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & eval("obj.Ind.premature" + obj.Ports(j)) );
                        n_late      = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & eval("obj.Ind.late"      + obj.Ports(j)) );
                        n_wrong     = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & eval("obj.Ind.wrong"     + obj.Ports(j)) );
                        n_legit     = n_correct + n_premature + n_late + n_wrong;

                        NumTrialsSorted(ind) = n_legit;

                        CorrectRatio(ind)   = 100 * n_correct   / n_legit;
                        PrematureRatio(ind) = 100 * n_premature / n_legit;
                        LateRatio(ind)      = 100 * n_late      / n_legit;
                        WrongRatio(ind)     = 100 * n_wrong     / n_legit;

                    else % group by port
                        n_correct   = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & obj.Ind.correct);
                        n_premature = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & obj.Ind.premature);
                        n_late      = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & obj.Ind.late);
                        n_wrong     = sum( obj.Stage==1 & obj.Cued==obj.CueUncue(i) & obj.Ind.wrong);
                        n_legit     = n_correct + n_premature + n_late + n_wrong;

                        NumTrialsSorted(ind) = n_legit;

                        CorrectRatio(ind)   = 100 * n_correct   / n_legit;
                        PrematureRatio(ind) = 100 * n_premature / n_legit;
                        LateRatio(ind)      = 100 * n_late      / n_legit;
                        WrongRatio(ind)     = 100 * n_wrong     / n_legit;
                    end
                end
            end

            Subjects = repmat(string(obj.Subject), length(Cued_this), 1);
            Sessions = repmat(string(obj.Session), length(Cued_this), 1);

            rt_table = table(Subjects, Sessions, Cued_this, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = rt_table;
        end

        function value = get.PerformanceTrackCue(obj)

            id_cue = obj.Cued==1;
            num_trials = sum(id_cue);

            WinSize  = floor(num_trials / 5);
            StepSize = max(1, floor(WinSize / 5));

            CountStart = 1;
            WinPos     = [];

            CorrectRatio    = [];
            WrongRatio      = [];
            PrematureRatio  = [];
            LateRatio       = [];

            correct_cue = obj.Ind.correct(id_cue);
            premature_cue = obj.Ind.premature(id_cue);
            wrong_cue = obj.Ind.wrong(id_cue);
            late_cue = obj.Ind.late(id_cue);
            center_pokes = obj.TrialStartTime(id_cue) + cellfun(@(x)x(1), obj.CentPokeInTime(id_cue)); % only count the first one

            while CountStart+WinSize-1 < num_trials

                thisWin = CountStart:(CountStart+WinSize-1);

                CorrectRatio    = [CorrectRatio    ;  100 * sum(correct_cue(thisWin))   / WinSize];
                WrongRatio      = [WrongRatio      ;  100 * sum(wrong_cue(thisWin))     / WinSize];
                PrematureRatio  = [PrematureRatio  ;  100 * sum(premature_cue(thisWin)) / WinSize];
                LateRatio       = [LateRatio       ;  100 * sum(late_cue(thisWin))      / WinSize];
                WinPos          = [WinPos          ;  center_pokes(thisWin(end))];
                
                CountStart = CountStart + StepSize;
            end

            Subjects = repmat(string(obj.Subject), length(WinPos), 1);
            Sessions = repmat(string(obj.Session), length(WinPos), 1);

            perf_track = table(Subjects, Sessions, WinPos, CorrectRatio, WrongRatio, PrematureRatio, LateRatio);
            value = perf_track;
        end

        function value = get.PerformanceTrackUncue(obj)

            id_uncue = obj.Cued==0;
            num_trials = sum(id_uncue);

            WinSize  = floor(num_trials / 5);
            StepSize = max(1, floor(WinSize / 5));

            CountStart = 1;
            WinPos     = [];

            CorrectRatio    = [];
            WrongRatio      = [];
            PrematureRatio  = [];
            LateRatio       = [];

            correct_uncue = obj.Ind.correct(id_uncue);
            premature_uncue = obj.Ind.premature(id_uncue);
            wrong_uncue = obj.Ind.wrong(id_uncue);
            late_uncue = obj.Ind.late(id_uncue);
            center_pokes = obj.TrialStartTime(id_uncue) + cellfun(@(x)x(1), obj.CentPokeInTime(id_uncue)); % only count the first one

            while CountStart+WinSize-1 < num_trials

                thisWin = CountStart:(CountStart+WinSize-1);

                CorrectRatio    = [CorrectRatio    ;  100 * sum(correct_uncue(thisWin))   / WinSize];
                WrongRatio      = [WrongRatio      ;  100 * sum(wrong_uncue(thisWin))     / WinSize];
                PrematureRatio  = [PrematureRatio  ;  100 * sum(premature_uncue(thisWin)) / WinSize];
                LateRatio       = [LateRatio       ;  100 * sum(late_uncue(thisWin))      / WinSize];
                WinPos          = [WinPos          ;  center_pokes(thisWin(end))];
                
                CountStart = CountStart + StepSize;
            end

            Subjects = repmat(string(obj.Subject), length(WinPos), 1);
            Sessions = repmat(string(obj.Session), length(WinPos), 1);

            perf_track = table(Subjects, Sessions, WinPos, CorrectRatio, WrongRatio, PrematureRatio, LateRatio);
            value = perf_track;
        end

        function value = get.BehavTable(obj)
            InitInTime  = cellfun(@(x) x(1),   obj.InitPokeInTime);
            InitOutTime = cellfun(@(x) x(end), obj.InitPokeOutTime);
            CentInTime  = cellfun(@(x) x(1),   obj.CentPokeInTime);
            CentOutTime = cellfun(@(x) x(end), obj.CentPokeOutTime);

            RWthis = obj.RW;
            if isempty(RWthis)
                RWthis = nan(obj.NumTrials, 1);
            end

            SessionDate  = repmat(string(obj.Session), obj.NumTrials, 1);
            SessionStart = repmat(string(obj.SessionStartTime), obj.NumTrials, 1);

            Subjects = repmat(string(obj.Subject), obj.NumTrials, 1);

            behav_table = table( ...
                Subjects, SessionDate, SessionStart, obj.Trials, obj.TrialStartTime, obj.Stage, ...
                InitInTime, InitOutTime, CentInTime, CentOutTime, obj.ChoicePokeTime, obj.ChoiceCueTime, obj.TriggerCueTime, ...
                obj.PortCorrect, obj.PortChosen, obj.FP, RWthis, obj.Outcome, ...
                obj.ShuttleTime, obj.HoldDuration, obj.RT, obj.MovementTime, obj.ChoiceTime, obj.Cued, ...
                'VariableNames', ...
                {'Subjects', 'SessionDate', 'SessionStartTime', 'Trials', 'TrialStartTime', 'Stage', ...
                'InitInTime', 'InitOutTime', 'CentInTime', 'CentOutTime', 'ChoicePokeTime', 'ChoiceCueTime', 'TriggerCueTime', ...
                'PortCorrect', 'PortChosen', 'FP', 'RW', 'Outcome', ...
                'ShuttleTime', 'HoldDuration', 'RT', 'MovementTime', 'ChoiceTime', 'Cued'});

            value = behav_table;
        end

        %%
        function save(obj, savepath)
            if nargin<2
                savepath = obj.BpodFilePath;
            end
            save(fullfile(savepath, ['GPSSessionClass_' obj.Task '_' upper(obj.Subject) '_' obj.Session]), 'obj');
        end

        function updateANMInfo(obj, savepath)
            anm_info = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");

            session_ind = strcmp(string(anm_info.Session), obj.Session);

            anm_info.Task(session_ind)      = obj.Task;
            anm_info.BpodFile(session_ind)  = obj.BpodFile;
            

            if nargin<2
                savepath = obj.BpodFilePath;
            end
            anm_info.SessionFolder(session_ind) = savepath;
            
            session_file = fullfile(savepath, ['GPSSessionClass_' obj.Task '_' upper(obj.Subject) '_' obj.Session]);
            anm_info.SessionClassFile(session_ind) = session_file;

            anm_info.UpdateTime(session_ind) = string(datetime());

            writetable(anm_info, obj.ANMInfoFile, "Sheet", obj.Subject);
        end

        function print(obj, targetDir)
            savename = fullfile(obj.BpodFilePath, ['GPSSessionClass_' obj.Task '_' upper(obj.Subject) '_' obj.Session]);
            hf = obj.plot();
            print(hf,'-dpdf', savename, '-bestfit')
            print(hf,'-dpng', savename)
            saveas(hf, savename, 'fig')

            if nargin==2
                % check if targetDir exists
                if ~contains(targetDir, '/') && ~contains(targetDir, '\')
                    % so it is a relative path
                    if ~exist(targetDir, 'dir')
                        mkdir(targetDir)
                    end
                end
                savename = fullfile(targetDir, ['GPSSessionClass_' obj.Task '_' upper(obj.Subject) '_' obj.Session]);
                print(hf,'-dpdf', savename, '-bestfit')
                print(hf,'-dpng', savename)
                saveas(hf, savename, 'fig')
            end
        end

        function fig = plot(obj)
            try
                set_matlab_default
            catch
                disp('You do not have "set_matlab_default"' )
            end

            fig = feval([obj.Task '.plotSession'], obj);
        end

        %% Data processing

        function data_sorted = sortData(obj, variable, perf)
            % Sort input data (RT, MovementTime, ChoiceTime, HoldDuration) by (CueUncue, Port).
            % inputs:
            %       variable: string (1,1), should be a member in ["RT", "MovementTime", "ChoiceTime", "HoldDuration"]
            %       perf: string (1,1), should be a member in ["correct", "late", "wrong", "premature", "port"]; if "port" is used, all trials will be taken

            data_origin = obj.(variable);
            data_sorted = cell(length(obj.CueUncue), length(obj.Ports));

            [~, ~, indrmv] = rmoutliers_custome(data_origin);

            for port = 1:length(obj.Ports)
                for cued = 1:length(obj.CueUncue)
                    ind_this = obj.Cued==obj.CueUncue(cued) & obj.Stage==1 & eval("obj.Ind." + perf + obj.Ports(port));
                    data_this = data_origin(setdiff(find(ind_this), indrmv));

                    data_sorted{cued, port} = data_this;
                end
            end
        end % sortData

        function data_pdf = getPDF(obj, variable, perf)

            data_sorted = obj.sortData(variable, perf);

            binEdges = obj.Bins.(variable);

            data_pdf = cell(length(obj.CueUncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    data_this = data_sorted{i, j};
                    if length(data_this) > 5
                        data_pdf{i, j} = ksdensity(data_this, binEdges, 'Function', 'pdf');
                    else
                        data_pdf{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
        end % getPDF

        function data_cdf = getCDF(obj, variable, perf)

            data_sorted = obj.sortData(variable, perf);

            binEdges = obj.Bins.(variable);

            data_cdf = cell(length(obj.CueUncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    data_this = data_sorted{i, j};
                    if length(data_this) > 5
                        data_cdf{i, j} = ksdensity(data_this, binEdges, 'Function', 'cdf');
                    else
                        data_cdf{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
        end % getCDF

        function stats = getStat(obj, variable, perf)
            % Get statistics from input data (RT, MovementTime, ChoiceTime, HoldDuration) by (FP, Port).
            % inputs:
            %       variable: string (1,1), should be a member in ["RT", "MovementTime", "ChoiceTime", "HoldDuration"]
            %       perf: string (1,1), should be a member in ["correct", "late", "wrong", "premature", "port"]; if "port" is used, all trials will be taken

            data_origin = obj.(variable);
            data_sorted = obj.sortData(variable, perf);

            % find iqr and quartiles
            if perf=="port"
                data_in = data_origin(obj.Stage==1);
            else
                data_in = data_origin(obj.Stage==1 & obj.Ind.(perf));
            end
            data_in(isnan(data_in)) = [];
            [data_2575] = prctile(data_in, [25, 75]);
            interq = data_2575(2) - data_2575(1);
            c = 5; % threshold for removing outliers

            % Set estimands
            height_this = (length(obj.CueUncue)+1)*length(obj.Ports);

            count = 0;

            thisCued = zeros(height_this, 1);
            Port = strings(height_this, 1);
            N = zeros(height_this, 1);

            Mean = zeros(height_this, 1);
            STD  = zeros(height_this, 1);
            SEM  = zeros(height_this, 1);

            Median = zeros(height_this, 1);
            Median_kde = zeros(height_this, 1);
            Q1 = zeros(height_this, 1);
            Q3 = zeros(height_this, 1);

            % Only include correct trials
            % PortL for all FPs
            data_this = data_origin(obj.Stage==1 & obj.Ind.(perf + "L"));
            data_this(isnan(data_this)) = [];
            data_this(data_this>data_2575(2)+interq*c | data_this<data_2575(1)-interq*c) = [];
            stat_this = calDur(data_this*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);

            count = count + 1;

            thisCued(count) = -1;
            Port(count) = "L";
            N(count) = length(data_this);

            Mean(count) = mean(data_this, 'omitnan');
            STD(count)  = std(data_this, 'omitnan');
            SEM(count)  = STD(count) / sqrt(N(count));

            Median(count) = stat_this.median * 0.001;
            Median_kde(count) = stat_this.median_ksdensity * 0.001;
            Q1(count) = prctile(data_this, 25);
            Q3(count) = prctile(data_this, 75);

            % PortR for all FPs
            data_this = data_origin(obj.Stage==1 & obj.Ind.(perf + "R"));
            data_this(isnan(data_this)) = [];
            data_this(data_this>data_2575(2)+interq*c | data_this<data_2575(1)-interq*c) = [];
            stat_this = calDur(data_this*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);

            count = count + 1;

            thisCued(count) = -1;
            Port(count) = "R";
            N(count) = length(data_this);

            Mean(count) = mean(data_this, 'omitnan');
            STD(count)  = std(data_this, 'omitnan');
            SEM(count)  = STD(count) / sqrt(N(count));

            Median(count) = stat_this.median * 0.001;
            Median_kde(count) = stat_this.median_ksdensity * 0.001;
            Q1(count) = prctile(data_this, 25);
            Q3(count) = prctile(data_this, 75);

            % By (Port, FP)
            for port = 1:length(obj.Ports) % Port
                for cued = 1:length(obj.CueUncue) % FP
                    data_this = data_sorted{cued, port};
                    data_this(isnan(data_this)) = [];
                    data_this(data_this>data_2575(2)+interq*c | data_this<data_2575(1)-interq*c) = [];
                    stat_this = calDur(data_this*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);

                    count = count + 1;

                    thisCued(count) = obj.CueUncue(cued);
                    Port(count) = obj.Ports(port);
                    N(count) = length(data_this);

                    Mean(count) = mean(data_this, 'omitnan');
                    STD(count)  = std(data_this, 'omitnan');
                    SEM(count)  = STD(count) / sqrt(N(count));

                    Median(count) = stat_this.median * 0.001;
                    Median_kde(count) = stat_this.median_ksdensity * 0.001;
                    Q1(count) = prctile(data_this, 25);
                    Q3(count) = prctile(data_this, 75);
                end
            end

            IQR = Q3 - Q1;

            Subjects = repmat(string(obj.Subject), height_this, 1);
            Sessions = repmat(string(obj.Session), height_this, 1);

            stats = table(Subjects, Sessions, thisCued, Port, N, Mean, STD, SEM, Median, Median_kde, IQR, Q1, Q3);
        end % getStat

    end % Methods

end % GPSSessionKbClass