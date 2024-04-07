classdef GPSSessionClass
    % Turn SessionData file into a class
    %   Detailed explanation goes here
    %

    properties
        % Experiment information
        BpodFile
        BpodFilePath
        BpodFileName
        Subject % Name of the animal, extrated from data file name

        ANMInfoFile;
        ANMInfo;
        SessionInfo;

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

        RTPDF = []; % reaction time
        RTCDF = [];
        MTPDF = []; % movement time
        MTCDF = [];
        HDPDF = []; % hold duration
        HDCDF = [];
        CTPDF = []; % choice time
        CTCDF = [];
        LogSTPDF = []; % Log10 shuttle time
        LogSTCDF = [];
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
        BandWidth = .05;
    end

    properties (Dependent)
        Strain % Strain of the animal, set manually
        Gender
        Experimenter % Name (initials) of the Experimenter, set manually
        Treatment % Set manually
        Dose % Set manually
        Label
        SortLabels

        Stage % warm-up or 3FP
        Engaged
        InitFail

        Ind
        Bins

        ShuttleTime % check
        ShuttleTimeSplit
        ShuttleTimeDistribution % check
        ShuttleTimeStat

        MovementTime
        MovementTimeSorted
        MovementTimeDistribution
        MovementTimeStat

        ChoiceTime
        ChoiceTimeSorted
        ChoiceTimeDistribution
        ChoiceTimeStat

        RT % check
        RTSorted
        RTDistribution % check
        RTStat % check

        HoldDuration %
        HoldDurationSorted
        HoldDurationDistribution
        HoldDurationStat

        Interruption

        Performance % Performance table, grouped by PortCorrect (and FP)
        PerformanceTrack

        BehavTable
    end

    methods
        %% Initiate
        function obj = GPSSessionClass(BpodFile, AnmInfoFile, CalCI)
            % Process SessionData to get properties
            load(BpodFile, 'SessionData');
            obj.BpodFile = BpodFile;
            [obj.BpodFilePath, obj.BpodFileName] = fileparts(obj.BpodFile);

            % Rat name
            obj.Subject = extractBefore(obj.BpodFileName, '_');
            Protocol = extractAfter(obj.BpodFileName, [obj.Subject '_']);
            obj.Session = Protocol(end-14:end-7);
            switch Protocol(1:end-16) % e.g. 'NEW_03_Wait3FPFlash'
                case {'NEW_01_Autoshaping', 'GPS_01_Autoshaping'}
                    obj.Task = 'Autoshaping';
                case {'GPS_02_Wait1Hold', 'GPS_04_Wait1Hold'}
                    obj.Task = 'Wait1Hold';
                case {'GPS_02_Wait1HoldSRT'}
                    obj.Task = 'Wait1HoldSRT';
                case {'GPS_04_Wait1HoldFlash', 'GPS_03_Wait1HoldFlash'}
                    obj.Task = 'Wait1HoldCRT';
                case {'GPS_03_Wait2HoldSRT'}
                    obj.Task = 'Wait2HoldSRT';
                case {'GPS_05_Wait2HoldFlash', 'GPS_04_Wait2HoldFlash'}
                    obj.Task = 'Wait2HoldCRT';
                case {'GPS_05_3FPHoldFlash', 'GPS_06_3FPHoldFlash'}
                    obj.Task = 'ThreeFPHoldCRT';
                case {'GPS_06_3FPHoldSRT'}
                    obj.Task = 'ThreeFPHoldSRT';
                case {'GPS_07_3FPHoldWM'}
                    obj.Task = 'ThreeFPHoldWM';
                case {'GPS_08_KornblumHold'}
                    obj.Task = 'KornblumSRT';
            end

            obj.ANMInfoFile = AnmInfoFile;
            obj.ANMInfo = readtable(obj.ANMInfoFile, "Sheet", "ANM", "TextType", "string");
            obj.ANMInfo = obj.ANMInfo(strcmp(obj.ANMInfo.Name, obj.Subject), :);
            obj.SessionInfo = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");
            obj.SessionInfo = obj.SessionInfo(strcmp(string(obj.SessionInfo.Session), obj.Session), :);

            % Session meta-information
            obj.SessionStartTime = SessionData.Info.SessionStartTime_UTC;

            % Get paradigm specific trial information
            obj = feval([obj.Task '.getTrialInfo'], obj, SessionData);

            obj = obj.getAllKDEs(CalCI);
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

        function obj = set.Dose(obj, dose)
            % Manually set the dose of injection
            if isnumeric(dose)
                obj.Dose = dose;
            else
                error('"dose" must be a scalar')
            end
        end

        function obj = set.Label(obj, label)
            % Manually set the treatment of this session, NaN/Saline/DCZ
            if ismember(label, {'None', 'Control', 'Chemo', 'Lesion'})
                obj.Treatment = string(label);
            else
                error('Label can only be: None, Control, Chemo, Lesion')
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
            strain = obj.ANMInfo.Strain;
            value = strain;
        end

        function value = get.Gender(obj)
            gender = obj.ANMInfo.Gender;
            value = gender;
        end

        function value = get.Treatment(obj)
            treatment = obj.SessionInfo.Treatment;
            value = treatment;
        end

        function value = get.Dose(obj)
            dose = obj.SessionInfo.Dose;
            value = dose;
        end

        function value = get.Label(obj)
            label = obj.SessionInfo.Label;
            value = label;
        end

        function value = get.Experimenter(obj)
            experimenter = obj.SessionInfo.Experimenter;
            value = experimenter;
        end

        function value = get.SortLabels(obj)
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

        %% Check whether the animal was engaged in the task
        function value = get.Engaged(obj)

            median_shuttle_time = median(obj.ShuttleTime);
            iqr_shuttle_time = diff(prctile(obj.ShuttleTime, [25, 75]));

            %             thres_dur = 0.1;

            in_task = obj.ShuttleTime < median_shuttle_time + 5*iqr_shuttle_time;
            value = in_task;
        end

        function value = get.InitFail(obj)
            init_fail = obj.HoldDuration < 0.5;
            value = init_fail;
        end

        %% Time interval info (same way to get for each paradigm)
        %% Shuttle time
        function value = get.ShuttleTime(obj)
            % find the last poke out time before center poke
            shuttle_time = cellfun(@(x, y) y(1)-x(end), obj.InitPokeOutTime, obj.CentPokeInTime);
            value = shuttle_time;
        end

        function value = get.ShuttleTimeSplit(obj)
            value = obj.splitData("ShuttleTime");
        end

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
            reaction_time(reaction_time<0.05) = nan;
            value = reaction_time;
        end

        function value = get.RTSorted(obj)
            value = obj.sortData("RT", "port");
        end

        function value = get.RTStat(obj)
            value = obj.getStat("RT", "port");
        end

        function value = get.RTDistribution(obj)
            value = obj.getDistr("RT", "port");
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

        function value = get.HoldDurationStat(obj)
            value = obj.getStat("HoldDuration", "port");
        end

        function value = get.HoldDurationDistribution(obj)
            value = obj.getDistr("HoldDuration", "port");
        end

        %% Movement time
        function value = get.MovementTime(obj)
            movement_time = obj.ChoicePokeTime - cellfun(@(x) x(end), obj.CentPokeOutTime);
            value = movement_time;
        end

        function value = get.MovementTimeSorted(obj)
            value = obj.sortData("MovementTime", "correct");
        end

        function value = get.MovementTimeStat(obj)
            value = obj.getStat("MovementTime", "correct");
        end

        function value = get.MovementTimeDistribution(obj)
            value = obj.getDistr("MovementTime", "correct");
        end

        %% Choice time
        function value = get.ChoiceTime(obj)
            choice_time = obj.RT + obj.MovementTime;
            value = choice_time;
        end

        function value = get.ChoiceTimeSorted(obj)
            value = obj.sortData("ChoiceTime", "correct");
        end

        function value = get.ChoiceTimeStat(obj)
            value = obj.getStat("ChoiceTime", "correct");
        end

        function value = get.ChoiceTimeDistribution(obj)
            value = obj.getDistr("ChoiceTime", "correct");
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
            hold_duration     = [];
            outcome           = [];

            valid_trials = find(obj.Stage);
            for t = 1:length(valid_trials)
                inter_on_collect  = [inter_on_collect  ; inter_on{t}];
                inter_dur_collect = [inter_dur_collect ; inter_dur{t}];
                trial             = [trial             ; repmat(t, length(inter_on{t}), 1)];
                port_correct      = [port_correct      ; repmat(obj.Ports(obj.PortCorrect(t)) , length(inter_on{t}), 1)];
                fp                = [fp                ; repmat(obj.FP(t)           , length(inter_on{t}), 1)];
                hold_duration     = [hold_duration     ; repmat(obj.HoldDuration(t) , length(inter_on{t}), 1)];
                outcome           = [outcome           ; repmat(obj.Outcome(t)      , length(inter_on{t}), 1)];
            end

            Subjects = repmat(string(obj.Subject), length(trial), 1);
            Sessions = repmat(string(obj.Session), length(trial), 1);

            inter = table(Subjects, Sessions, trial, inter_on_collect, inter_dur_collect, port_correct, fp, hold_duration, outcome, ...
                'VariableNames', {'Subjects', 'Sessions', 'Trials', 'On', 'Dur', 'PortCorrect', 'FP', 'HoldDuration', 'Outcome'});
            inter(inter.Dur<0.001, :) = [];

            value = inter;
        end

        %%
        function value = get.Performance(obj)
            MixedFP_this = obj.MixedFP;
            if size(MixedFP_this, 1)==1
                MixedFP_this = MixedFP_this';
            end
            Foreperiod = [MixedFP_this; 0];
            TargetPort = [repmat("L", length(obj.MixedFP)+1, 1); repmat("R", length(obj.MixedFP)+1, 1); "Both"]; % 1 or 2
            Foreperiod = [repmat(Foreperiod, 2, 1); 0];
            NumTrialsSorted = zeros(length(Foreperiod), 1);

            % group by foreperiod
            CorrectRatio = zeros(length(Foreperiod), 1);
            PrematureRatio = zeros(length(Foreperiod), 1);
            LateRatio = zeros(length(Foreperiod), 1);
            WrongRatio = zeros(length(Foreperiod), 1);

            % group by target port
            format short

            for i = 1:length(obj.MixedFP)+1
                for j = 1:2 % Target port
                    ind = i + (j - 1) * (length(obj.MixedFP) + 1);
                    if mod(i, length(obj.MixedFP)+1)~=0 % group by FP and Target port
                        n_correct   = sum( obj.Stage==1 & obj.FP==obj.MixedFP(i) & eval("obj.Ind.correct"   + obj.Ports(j)) );
                        n_premature = sum( obj.Stage==1 & obj.FP==obj.MixedFP(i) & eval("obj.Ind.premature" + obj.Ports(j)) );
                        n_late      = sum( obj.Stage==1 & obj.FP==obj.MixedFP(i) & eval("obj.Ind.late"      + obj.Ports(j)) );
                        n_wrong     = sum( obj.Stage==1 & obj.FP==obj.MixedFP(i) & eval("obj.Ind.wrong"     + obj.Ports(j)) );
                        n_legit     = n_correct + n_premature + n_late + n_wrong;

                        NumTrialsSorted(ind) = n_legit;

                        CorrectRatio(ind)   = 100 * n_correct   / n_legit;
                        PrematureRatio(ind) = 100 * n_premature / n_legit;
                        LateRatio(ind)      = 100 * n_late      / n_legit;
                        WrongRatio(ind)     = 100 * n_wrong     / n_legit;

                    else % group by port
                        n_correct   = sum( obj.Stage==1 & eval("obj.Ind.correct"   + obj.Ports(j)) );
                        n_premature = sum( obj.Stage==1 & eval("obj.Ind.premature" + obj.Ports(j)) );
                        n_late      = sum( obj.Stage==1 & eval("obj.Ind.late"      + obj.Ports(j)) );
                        n_wrong     = sum( obj.Stage==1 & eval("obj.Ind.wrong"     + obj.Ports(j)) );
                        n_legit     = n_correct + n_premature + n_late + n_wrong;

                        NumTrialsSorted(ind) = n_legit;

                        CorrectRatio(ind)   = 100 * n_correct   / n_legit;
                        PrematureRatio(ind) = 100 * n_premature / n_legit;
                        LateRatio(ind)      = 100 * n_late      / n_legit;
                        WrongRatio(ind)     = 100 * n_wrong     / n_legit;
                    end
                end
            end

            % total
            ind = ind+1;
            n_correct   = sum(obj.Stage==1 & obj.Ind.correct);
            n_premature = sum(obj.Stage==1 & obj.Ind.premature);
            n_late      = sum(obj.Stage==1 & obj.Ind.late);
            n_wrong     = sum(obj.Stage==1 & obj.Ind.wrong);
            n_legit     = n_correct + n_premature + n_late + n_wrong;

            NumTrialsSorted(ind) = n_legit;

            CorrectRatio(ind)   = 100 * n_correct   / n_legit;
            PrematureRatio(ind) = 100 * n_premature / n_legit;
            LateRatio(ind)      = 100 * n_late      / n_legit;
            WrongRatio(ind)     = 100 * n_wrong     / n_legit;

            Subjects = repmat(string(obj.Subject), length(Foreperiod), 1);
            Sessions = repmat(string(obj.Session), length(Foreperiod), 1);

            rt_table = table(Subjects, Sessions, Foreperiod, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = rt_table;
        end

        function value = get.PerformanceTrack(obj)

            center_pokes = obj.TrialStartTime + cellfun(@(x)x(1), obj.CentPokeInTime); % only count the first one

            WinSize  = floor(obj.NumTrials / 5);
            StepSize = max(1, floor(WinSize / 5));

            CountStart = 1;
            WinPos     = [];

            CorrectRatio    = [];
            WrongRatio      = [];
            PrematureRatio  = [];
            LateRatio       = [];

            while CountStart+WinSize-1 < obj.NumTrials

                thisWin = CountStart:(CountStart+WinSize-1);

                CorrectRatio    = [CorrectRatio    ;  100 * sum(obj.Ind.correct(thisWin))   / WinSize];
                WrongRatio      = [WrongRatio      ;  100 * sum(obj.Ind.wrong(thisWin))     / WinSize];
                PrematureRatio  = [PrematureRatio  ;  100 * sum(obj.Ind.premature(thisWin)) / WinSize];
                LateRatio       = [LateRatio       ;  100 * sum(obj.Ind.late(thisWin))      / WinSize];
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
                obj.ShuttleTime, obj.HoldDuration, obj.RT, obj.MovementTime, obj.ChoiceTime, ...
                'VariableNames', ...
                {'Subjects', 'SessionDate', 'SessionStartTime', 'Trials', 'TrialStartTime', 'Stage', ...
                'InitInTime', 'InitOutTime', 'CentInTime', 'CentOutTime', 'ChoicePokeTime', 'ChoiceCueTime', 'TriggerCueTime', ...
                'PortCorrect', 'PortChosen', 'FP', 'RW', 'Outcome', ...
                'ShuttleTime', 'HoldDuration', 'RT', 'MovementTime', 'ChoiceTime'});

            value = behav_table;
        end

        %% Information saving
        function save(obj, savepath)
            if nargin<2
                savepath = obj.BpodFilePath;
            end
            save(fullfile(savepath, ['GPSSessionClass_' obj.Task '_' upper(obj.Subject) '_' obj.Session]), 'obj');
        end

        function updateANMInfo(obj, savepath)
            if nargin<2
                savepath = obj.BpodFilePath;
            end

            anm_info = readtable(obj.ANMInfoFile, "Sheet", obj.Subject, "TextType", "string");
            session_ind = strcmp(string(anm_info.Session), obj.Session);
            session_file = fullfile(savepath, ['GPSSessionClass_' obj.Task '_' upper(obj.Subject) '_' obj.Session]);

            anm_info.Task(session_ind)      = obj.Task;
            anm_info.BpodFile(session_ind)  = obj.BpodFile;
            anm_info.SessionFolder(session_ind) = savepath;
            anm_info.SessionClassFile(session_ind) = session_file;
            anm_info.UpdateTime(session_ind) = string(datetime());

            writetable(anm_info, obj.ANMInfoFile, "Sheet", obj.Subject);
        end

        %% Polt and Print
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
            % Sort input data (RT, MovementTime, ChoiceTime, HoldDuration) by (FP, Port).
            % inputs:
            %       variable: string (1,1), should be a member in ["RT", "MovementTime", "ChoiceTime", "HoldDuration"]
            %       perf: string (1,1), should be a member in ["correct", "late", "wrong", "premature", "port"]; if "port" is used, all trials will be taken

            data_origin = obj.(variable);
            data_sorted = cell(length(obj.MixedFP), length(obj.Ports));

            ind_valid = ~isnan(data_origin);
%             [~, ~, indrmv] = rmoutliers_custome(data_origin);

            for port = 1:length(obj.Ports)
                for fp = 1:length(obj.MixedFP)
                    ind_this = obj.FP==obj.MixedFP(fp) & obj.Stage==1 & eval("obj.Ind." + perf + obj.Ports(port));
%                     data_this = data_origin(setdiff(find(ind_this), indrmv));

                    data_sorted{fp, port} = data_origin(ind_this & ind_valid);
                end
            end
        end % sortData

        function data_splited = splitData(obj, variable)
            % Sort input data (RT, MovementTime, ChoiceTime, HoldDuration) by (FP, Port).
            % inputs:
            %       variable: string (1,1), should be a member in ["RT", "MovementTime", "ChoiceTime", "HoldDuration"]
            %       perf: string (1,1), should be a member in ["correct", "late", "wrong", "premature", "port"]; if "port" is used, all trials will be taken

            data_origin = obj.(variable);
            data_splited = cell(3, 1);

            data_this = data_origin(obj.Stage==1);
            num_data = length(data_this);
            ind = 1:num_data;
%             [~, ~, indrmv] = rmoutliers_custome(data_origin);

            for i = 1:3
                ind_range = [i-1 i] * num_data/3;
                ind_this = ind(ind>ind_range(1) & ind<=ind_range(2));
                data_splited{i} = data_origin(ind_this);
            end
        end % splitData

        %% get PDFs and CDFs
        function data_pdf = getPDF(obj, data, bin_edges, var_name, cal_ci)
            sz = size(data);
            data_pdf = cell(sz(1), sz(2));

            if cal_ci
                fprintf("... Calculate 95CI for PDF of %s ... \n", var_name);
            end
            for i = 1:sz(1)
                for j = 1:sz(2)
                    data_this = data{i, j};
                    data_pdf{i, j} = obj.calKDE(data_this, bin_edges, 'pdf', cal_ci);
                end
            end
        end % getPDF

        function data_cdf = getCDF(obj, data, bin_edges, var_name, cal_ci)
            sz = size(data);
            data_cdf = cell(sz(1), sz(2));

            if cal_ci
                fprintf("... Calculate 95CI for CDF of %s ... \n", var_name);
            end
            for i = 1:sz(1)
                for j = 1:sz(2)
                    data_this = data{i, j};
                    data_cdf{i, j} = obj.calKDE(data_this, bin_edges, 'cdf', cal_ci);
                end
            end
        end % getCDF

        function data_kde = calKDE(obj, data_this, bin_edges, func, cal_ci)
            kde = @(x) ksdensity(x, bin_edges, 'Function', func, 'Bandwidth', obj.BandWidth);

            data_kde.x = bin_edges;
            data_kde.ci = [];
            if length(data_this) > 5
                data_kde.f = kde(data_this);
                if cal_ci
                    data_kde.ci = bootci(1000, {kde, data_this}, 'type', 'cper', 'alpha', .05);
                end
            else
                data_kde.f = zeros(1, length(bin_edges));
            end
        end % calKDE

        function obj = getAllKDEs(obj, cal_ci)

            obj.RTPDF = obj.getPDF(obj.RTSorted, obj.Bins.RT, "RT", cal_ci);
            obj.RTCDF = obj.getCDF(obj.RTSorted, obj.Bins.RT, "RT", cal_ci);
            obj.MTPDF = obj.getPDF(obj.MovementTimeSorted, obj.Bins.MovementTime, "MT", cal_ci);
            obj.MTCDF = obj.getCDF(obj.MovementTimeSorted, obj.Bins.MovementTime, "MT", cal_ci);
            obj.HDPDF = obj.getPDF(obj.HoldDurationSorted, obj.Bins.HoldDuration, "HD", cal_ci);
            obj.HDCDF = obj.getCDF(obj.HoldDurationSorted, obj.Bins.HoldDuration, "HD", cal_ci);
            obj.CTPDF = obj.getPDF(obj.ChoiceTimeSorted, obj.Bins.ChoiceTime, "CT", cal_ci);
            obj.CTCDF = obj.getCDF(obj.ChoiceTimeSorted, obj.Bins.ChoiceTime, "CT", cal_ci);
            logST = cellfun(@log10, obj.ShuttleTimeSplit, 'UniformOutput', false);
            obj.LogSTPDF = obj.getPDF(logST, obj.Bins.ShuttleTimeLog, "LogST", cal_ci);
            obj.LogSTCDF = obj.getCDF(logST, obj.Bins.ShuttleTimeLog, "LogST", cal_ci);

        end % getAllKDEs

        %% Get statistics
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
            height_this = (length(obj.MixedFP)+1)*length(obj.Ports);
            stats = array2table(zeros(height_this, 13), 'VariableNames', ...
                {'Subjects', 'Sessions', 'thisFP', 'Port', 'N', ...
                'Mean', 'STD', 'SEM', ...
                'Median', 'Median_kde', 'Q1', 'Q3', 'IQR'});
            stats.Port = strings(height_this, 1);
            count = 0;

            % By (Port, FP)
            for port = 1:length(obj.Ports) % Port
                data_this = data_origin(obj.Stage==1 & obj.Ind.(perf + obj.Ports(port)));
                data_this(isnan(data_this)) = [];
                data_this(data_this>data_2575(2)+interq*c | data_this<data_2575(1)-interq*c) = [];
                stat_this = calDur(data_this*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);

                count = count + 1;

                stats.thisFP(count) = 0;
                stats.Port(count) = obj.Ports(port);
                stats.N(count) = sum(~isnan(data_this));

                stats.Mean(count) = stat_this.mean * .001;
                stats.STD(count)  = stat_this.std * .001;
                stats.SEM(count)  = stat_this.sem * .001;

                stats.Median(count) = stat_this.median * .001;
                stats.Median_kde(count) = stat_this.median_ksdensity * .001;
                stats.Q1(count) = stat_this.q1 * .001;
                stats.Q3(count) = stat_this.q3 * .001;

                for fp = 1:length(obj.MixedFP) % FP
                    data_this = data_sorted{fp, port};
                    data_this(isnan(data_this)) = [];
                    data_this(data_this>data_2575(2)+interq*c | data_this<data_2575(1)-interq*c) = [];
                    stat_this = calDur(data_this*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);

                    count = count + 1;

                    stats.thisFP(count) = obj.MixedFP(fp);
                    stats.Port(count) = obj.Ports(port);
                    stats.N(count) = sum(~isnan(data_this));

                    stats.Mean(count) = stat_this.mean * .001;
                    stats.STD(count)  = stat_this.std * .001;
                    stats.SEM(count)  = stat_this.sem * .001;

                    stats.Median(count) = stat_this.median * .001;
                    stats.Median_kde(count) = stat_this.median_ksdensity * .001;
                    stats.Q1(count) = stat_this.q1 * .001;
                    stats.Q3(count) = stat_this.q3 * .001;
                end
            end

            stats.IQR = stats.Q3 - stats.Q1;

            stats.Subjects = repmat(string(obj.Subject), height_this, 1);
            stats.Sessions = repmat(string(obj.Session), height_this, 1);

        end % getStat

        function data_distr = getDistr(obj, variable, perf)

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
            c = 5;

            bin_edges = obj.Bins.(variable);
            bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

            % Only include correct trials
            data_portL = data_origin(obj.Ind.(perf+"L"));
            data_portL(data_portL>data_2575(2)+interq*c | data_portL<data_2575(1)-interq*c) = [];
            data_portR = data_origin(obj.Ind.(perf+"R"));
            data_portR(data_portR>data_2575(2)+interq*c | data_portR<data_2575(1)-interq*c) = [];

            % compute RT from two ports separately
            num_counts = histcounts(data_portL, bin_edges);
            counts_portL = num_counts';
            num_counts = histcounts(data_portR, bin_edges);
            counts_portR = num_counts';
            RT_centers = bin_centers';

            switch obj.Task
                case {'Wait1Hold', 'Wait1HoldCRT', 'Wait2HoldCRT'}
                    data_distr = table(RT_centers, counts_portL, counts_portR);
                case {'ThreeFPHoldCRT', 'ThreeFPHoldSRT', 'ThreeFPHoldWM'}
                    FP_types = ["S", "M", "L"];

                    % compute RT from two ports and three FPs separately
                    for i = 1:length(obj.MixedFP)
                        for j = 1:2
                            iRTs = data_sorted{i, j};
                            iRTs(isnan(iRTs)) = [];
                            % removeOutliers
                            iRTs(iRTs>data_2575(2)+interq*c | iRTs<data_2575(1)-interq*c) = [];
                            num_counts = histcounts(iRTs, bin_edges);
                            num_counts = num_counts';
                            eval("counts_" + FP_types(i) + "_port" + obj.Ports(j) + "= num_counts;");
                        end
                    end
                    data_distr = table(RT_centers, counts_portL, counts_portR, counts_S_portL, counts_M_portL, counts_L_portL, ...
                        counts_S_portR, counts_M_portR, counts_L_portR);
            end
        end % getDistr

    end % Methods

end % GPSSessionClass