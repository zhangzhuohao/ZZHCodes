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
        ANMInfoFile = "E:\YuLab\Work\GPS\Data\ANMInfo.xlsx";
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
        MovementTimeDistribution
        MovementTimeStat

        ChoiceTime
        ChoiceTimeSorted
        ChoiceTimePDF
        ChoiceTimeCDF
        ChoiceTimeDistribution
        ChoiceTimeStat

        RT % check
        RTSorted
        RTPDF
        RTCDF
        RTDistribution % check
        RTStat % check

        HoldDuration %
        HoldDurationSorted
        HoldDurationPDF
        HoldDurationCDF
        HoldDurationDistribution
        HoldDurationStat

        Interruption

        Performance % Performance table, grouped by PortCorrect (and FP)
        PerformanceTrackCue
        PerformanceTrackUncue

        BehavTable
    end

    methods
        %% Initiate
        function obj = GPSSessionKbClass(BpodFile)
            % Process SessionData to get properties
            load(BpodFile, 'SessionData');
            obj.BpodFile = BpodFile;
            [obj.BpodFilePath, obj.BpodFileName] = fileparts(obj.BpodFile);

            % Rat name
            obj.Subject = extractBefore(obj.BpodFileName, '_');
            Protocol = extractAfter(obj.BpodFileName, [obj.Subject '_']);
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

        %% Time interval info (same way to get for each paradigm)
        %% Shuttle time
        function value = get.ShuttleTime(obj)
            % find the last poke out time before center poke
            shuttle_time = cellfun(@(x, y) y(1)-x(end), obj.InitPokeOutTime, obj.CentPokeInTime);
            value = shuttle_time;
        end

        %% Reaction time
        function value = get.RT(obj)
            % find the last poke out time before center poke
            reaction_time = cellfun(@(x) x(end), obj.CentPokeOutTime) - obj.TriggerCueTime;
            value = reaction_time;
        end

        %% Movement time
        function value = get.MovementTime(obj)
            movement_time = obj.ChoicePokeTime - cellfun(@(x) x(end), obj.CentPokeOutTime);
            value = movement_time;
        end

        function value = get.ChoiceTime(obj)
            choice_time = obj.RT + obj.MovementTime;
            value = choice_time;
        end

        %% Hold duration
        function value = get.HoldDuration(obj)
            % find the last poke out time before center poke
            hold_duration = cellfun(@(x, y) y(end) - x(1), obj.CentPokeInTime, obj.CentPokeOutTime);
            value = hold_duration;
        end

        %%
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

        %% Statistics
        %% Shuttle time statistics
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

        %% Reaction time statistics
        function value = get.RTSorted(obj)

            [~, ~, indrmv] = rmoutliers_custome(obj.RT);

            RTCollected = cell(length(obj.CueUncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    ind_thisFP = obj.Cued==obj.CueUncue(i) & obj.Stage==1 & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.RT(setdiff(find(ind_thisFP), indrmv));
                    RTCollected{i, j} = iRTs;
                end
            end
            value = RTCollected;
        end

        function value = get.RTPDF(obj)
            cueuncue = obj.CueUncue;
            binEdges = obj.Bins.RT;

            RT_PDF = cell(length(cueuncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iRTs = obj.RTSorted{i, j};
                    if length(iRTs) > 5
                        RT_PDF{i, j} = ksdensity(iRTs, binEdges, 'Function', 'pdf');
                    else
                        RT_PDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = RT_PDF;
        end

        function value = get.RTCDF(obj)
            cueuncue = obj.CueUncue;
            binEdges = obj.Bins.RT;

            RT_CDF = cell(length(cueuncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iRTs = obj.RTSorted{i, j};
                    if length(iRTs) > 5
                        RT_CDF{i, j} = ksdensity(iRTs, binEdges, 'Function', 'cdf');
                    else
                        RT_CDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = RT_CDF;
        end

        function value = get.RTStat(obj)

            datain = obj.RT(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            thisFP = [];
            thisCued = [];
            N = [];
            Port = [];
            Median = [];
            Median_kde = [];
            Q1 = [];
            Q3 = [];

            % Only include correct trials
            for j = 1:2 % Port
                for i = 1:length(obj.CueUncue) % FP
                    iRTs = obj.RTSorted{i, j};
                    iRTs(isnan(iRTs)) = [];
                    % removeOutliers
                    iRTs(iRTs>data2575(2)+interq*c | iRTs<data2575(1)-interq*c) = [];
                    
                    thisFP = [thisFP; unique(obj.FP)];
                    thisCued = [thisCued; obj.CueUncue(i)];
                    N = [N; length(iRTs)];
                    Port = [Port; obj.Ports(j)];
                    iRTOut = calDur(iRTs*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);
                    Median = [Median; iRTOut.median*0.001];
                    Median_kde = [Median_kde; iRTOut.median_ksdensity*0.001];
                    Q1 = [Q1; prctile(iRTs, 25)];
                    Q3 = [Q3; prctile(iRTs, 75)];
                end
            end
            IQR = Q3 - Q1;

            Subjects = repmat(string(obj.Subject), length(thisFP), 1);
            Sessions = repmat(string(obj.Session), length(thisFP), 1);

            value = table(Subjects, Sessions, thisFP, thisCued, Port, N, Median, Median_kde, IQR, Q1, Q3);
        end

        function value = get.RTDistribution(obj)
            datain = obj.RT(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            binEdges = obj.Bins.RT;
            binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

            % Only include correct trials
            RT_PortL = obj.RT(obj.Ind.correctL);
            RT_PortL(RT_PortL>data2575(2)+interq*c | RT_PortL<data2575(1)-interq*c) = [];
            RT_PortR = obj.RT(obj.Ind.correctR);
            RT_PortR(RT_PortR>data2575(2)+interq*c | RT_PortR<data2575(1)-interq*c) = [];

            % compute RT from two ports separately
            Ncounts = histcounts(RT_PortL, binEdges);
            Counts_PortL = Ncounts';
            Ncounts = histcounts(RT_PortR, binEdges);
            Counts_PortR = Ncounts';
            RT_Centers = binCenters';

            switch obj.Task
                case {'Wait1Hold', 'Wait1HoldCRT', 'Wait2HoldCRT'}
                    value = table(RT_Centers, Counts_PortL, Counts_PortR);
                case {'ThreeFPHoldCRT', 'ThreeFPHoldSRT'}
                    FP_types = ["S", "M", "L"];

                    % compute RT from two ports and three FPs separately
                    for i = 1:length(obj.MixedFP)
                        for j = 1:2
                            iRTs = obj.RTSorted{i, j};
                            iRTs(isnan(iRTs)) = [];
                            % removeOutliers
                            iRTs(iRTs>data2575(2)+interq*c | iRTs<data2575(1)-interq*c) = [];
                            Ncounts = histcounts(iRTs, binEdges);
                            Ncounts = Ncounts';
                            eval("Counts_" + FP_types(i) + "_Port" + obj.Ports(j) + "= Ncounts;");
                        end
                    end
                    value = table(RT_Centers, Counts_PortL, Counts_PortR, Counts_S_PortL, Counts_M_PortL, Counts_L_PortL, ...
                        Counts_S_PortR, Counts_M_PortR, Counts_L_PortR);
            end
        end

        %% Hold duration
        function value = get.HoldDurationSorted(obj)
            cueuncue = obj.CueUncue;

            [~, ~, indrmv] = rmoutliers_custome(obj.HoldDuration);

            HDCollected = cell(length(cueuncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    ind_thisFP = obj.Cued==obj.CueUncue(i) & obj.Stage==1 & eval("obj.Ind.port"+obj.Ports(j));
                    iHDs = obj.HoldDuration(setdiff(find(ind_thisFP), indrmv));
                    HDCollected{i, j} = iHDs;
                end
            end
            value = HDCollected;
        end

        function value = get.HoldDurationPDF(obj)
            cueuncue = obj.CueUncue;

            HDPDF = cell(length(cueuncue), length(obj.Ports));
            binEdges = obj.Bins.HoldDuration;
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iHDs = obj.HoldDurationSorted{i, j};
                    if length(iHDs) > 5
                        HDPDF{i, j} = ksdensity(iHDs, binEdges, 'Function', 'pdf');
                    else
                        HDPDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = HDPDF;
        end

        function value = get.HoldDurationCDF(obj)
            cueuncue = obj.CueUncue;

            HDCDF = cell(length(cueuncue), length(obj.Ports));
            bins = obj.Bins.HoldDuration;
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iHDs = obj.HoldDurationSorted{i, j};
                    if length(iHDs) > 5
                        HDCDF{i, j} = ksdensity(iHDs, bins, 'Function', 'cdf');
                    else
                        HDCDF{i, j} = zeros(1, length(bins));
                    end
                end
            end
            value = HDCDF;
        end

        function value = get.HoldDurationStat(obj)

            datain = obj.HoldDuration(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            thisFP = [];
            thisCued = [];
            N = [];
            Port = [];
            Median = [];
            Median_kde = [];
            Q1 = [];
            Q3 = [];

            % Only include correct trials
            for j = 1:2 % Port
                for i = 1:length(obj.CueUncue) % FP
                    iHDs = obj.HoldDurationSorted{i, j};
                    iHDs(isnan(iHDs)) = [];
                    % removeOutliers
                    iHDs(iHDs>data2575(2)+interq*c | iHDs<data2575(1)-interq*c) = [];
                    
                    thisFP = [thisFP; unique(obj.FP)];
                    thisCued = [thisCued; obj.CueUncue(i)];
                    N = [N; length(iHDs)];
                    Port = [Port; obj.Ports(j)];
                    iHDOut = calDur(iHDs*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);
                    Median = [Median; iHDOut.median*0.001];
                    Median_kde = [Median_kde; iHDOut.median_ksdensity*0.001];
                    Q1 = [Q1; prctile(iHDs, 25)];
                    Q3 = [Q3; prctile(iHDs, 75)];
                end
            end
            IQR = Q3 - Q1;

            Subjects = repmat(string(obj.Subject), length(thisFP), 1);
            Sessions = repmat(string(obj.Session), length(thisFP), 1);

            value = table(Subjects, Sessions, thisFP, thisCued, Port, N, Median, Median_kde, IQR, Q1, Q3);
        end

        function value = get.HoldDurationDistribution(obj)

            datain = obj.HoldDuration(obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            binEdges = obj.Bins.HoldDuration;
            binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

            % Only include correct trials
            HD_PortL = obj.HoldDuration(obj.Ind.portL);
            HD_PortL(HD_PortL>data2575(2)+interq*c | HD_PortL<data2575(1)-interq*c) = [];
            HD_PortR = obj.HoldDuration(obj.Ind.portR);
            HD_PortR(HD_PortR>data2575(2)+interq*c | HD_PortR<data2575(1)-interq*c) = [];

            % compute RT from two ports separately
            Ncounts = histcounts(HD_PortL, binEdges);
            Counts_PortL = Ncounts';
            Ncounts = histcounts(HD_PortR, binEdges);
            Counts_PortR = Ncounts';
            HD_Centers = binCenters';

            switch obj.Task
                case {'Wait1Hold', 'Wait1HoldCRT', 'Wait2HoldCRT'}
                    value = table(HD_Centers, Counts_PortL, Counts_PortR);
                case {'ThreeFPHoldCRT', 'ThreeFPHoldSRT'}
                    FP_types = ["S", "M", "L"];

                    % compute RT from two ports and three FPs separately
                    for j = 1:length(obj.Ports)
                        for i = 1:length(obj.MixedFP)
                            iHDs = obj.HoldDurationSorted{i, j};
                            iHDs(isnan(iHDs)) = [];
                            % removeOutliers
                            iHDs(iHDs>data2575(2)+interq*c | iHDs<data2575(1)-interq*c) = [];
                            Ncounts = histcounts(iHDs, binEdges);
                            Ncounts = Ncounts';
                            eval("Counts_" + FP_types(i) + "_Port" + obj.Ports(j) + "= Ncounts;");
                        end
                    end
                    value = table(HD_Centers, Counts_PortL, Counts_PortR, Counts_S_PortL, Counts_M_PortL, Counts_L_PortL, ...
                        Counts_S_PortR, Counts_M_PortR, Counts_L_PortR);
            end
        end

        %% Movement time
        function value = get.MovementTimeSorted(obj)
            cueuncue = obj.CueUncue;

            [~, ~, indrmv] = rmoutliers_custome(obj.MovementTime);

            movement_time_sorted = cell(length(cueuncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    ind_thisFP = obj.Cued==obj.CueUncue(i) & obj.Stage==1 & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.MovementTime(setdiff(find(ind_thisFP), indrmv));
                    movement_time_sorted{i, j} = iMTs;
                end
            end
            value = movement_time_sorted;
        end

        function value = get.MovementTimePDF(obj)
            cueuncue = obj.CueUncue;

            MT_PDF = cell(length(cueuncue), length(obj.Ports));
            binEdges = obj.Bins.MovementTime;
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iMTs = obj.MovementTimeSorted{i, j};
                    if length(iMTs) > 5
                        MT_PDF{i, j} = ksdensity(iMTs, binEdges, 'Function', 'pdf');
                    else
                        MT_PDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = MT_PDF;
        end

        function value = get.MovementTimeCDF(obj)
            cueuncue = obj.CueUncue;

            MT_CDF = cell(length(cueuncue), length(obj.Ports));
            binEdges = obj.Bins.MovementTime;
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iMTs = obj.MovementTimeSorted{i, j};
                    if length(iMTs) > 5
                        MT_CDF{i, j} = ksdensity(iMTs, binEdges, 'Function', 'cdf');
                    else
                        MT_CDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = MT_CDF;
        end

        function value = get.MovementTimeStat(obj)
            datain = obj.MovementTime(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            thisFP = [];
            thisCued = [];
            N = [];
            Port = [];
            Median = [];
            Median_kde = [];
            Q1 = [];
            Q3 = [];

            % Only include correct trials
            for j = 1:2 % Port
                for i = 1:length(obj.CueUncue) % FP
                    iMTs = obj.MovementTimeSorted{i, j};
                    iMTs(isnan(iMTs)) = [];
                    % removeOutliers
                    iMTs(iMTs>data2575(2)+interq*c | iMTs<data2575(1)-interq*c) = [];
                    
                    thisFP = [thisFP; unique(obj.FP)];
                    thisCued = [thisCued; obj.CueUncue(i)];
                    N = [N; length(iMTs)];
                    Port = [Port; obj.Ports(j)];
                    iMTOut = calDur(iMTs*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);
                    Median = [Median; iMTOut.median*0.001];
                    Median_kde = [Median_kde; iMTOut.median_ksdensity*0.001];
                    Q1 = [Q1; prctile(iMTs, 25)];
                    Q3 = [Q3; prctile(iMTs, 75)];
                end
            end
            IQR = Q3 - Q1;

            Subjects = repmat(string(obj.Subject), length(thisFP), 1);
            Sessions = repmat(string(obj.Session), length(thisFP), 1);

            value = table(Subjects, Sessions, thisFP, thisCued, Port, N, Median, Median_kde, IQR, Q1, Q3);
        end

        function value = get.MovementTimeDistribution(obj)
            datain = obj.MovementTime(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            binEdges = 0:0.05:3;
            binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

            % Only include correct trials
            MT_PortL = obj.MovementTime(obj.Stage==1 & obj.Ind.correctL);
            MT_PortL(MT_PortL>data2575(2)+interq*c | MT_PortL<data2575(1)-interq*c) = [];
            MT_PortR = obj.MovementTime(obj.Stage==1 & obj.Ind.correctR);
            MT_PortR(MT_PortR>data2575(2)+interq*c | MT_PortR<data2575(1)-interq*c) = [];

            % compute RT from two ports separately
            Ncounts = histcounts(MT_PortL, binEdges);
            Counts_PortL = Ncounts';
            Ncounts = histcounts(MT_PortR, binEdges);
            Counts_PortR = Ncounts';
            MT_Centers = binCenters';

            switch obj.Task
                case {'Wait1Hold', 'Wait1HoldCRT', 'Wait2HoldCRT'}
                    movement_time_distr = table(MT_Centers, Counts_PortL, Counts_PortR);
                case {'ThreeFPHoldCRT', 'ThreeFPHoldSRT'}
                    FP_types = ["S", "M", "L"];

                    % compute RT from two ports and three FPs separately
                    for j = 1:length(obj.Ports)
                        for i = 1:length(obj.MixedFP)
                            iMTs = obj.MovementTimeSorted{i, j};
                            iMTs(isnan(iMTs)) = [];
                            % removeOutliers
                            iMTs(iMTs>data2575(2)+interq*c | iMTs<data2575(1)-interq*c) = [];
                            Ncounts = histcounts(iMTs, binEdges);
                            Ncounts = Ncounts';
                            eval("Counts_" + FP_types(i) + "_Port" + obj.Ports(j) + "= Ncounts;");
                        end
                    end
                    movement_time_distr = table(MT_Centers, Counts_PortL, Counts_PortR, Counts_S_PortL, Counts_M_PortL, Counts_L_PortL, ...
                        Counts_S_PortR, Counts_M_PortR, Counts_L_PortR);
            end
            value = movement_time_distr;
        end

        %% Choice time
        function value = get.ChoiceTimeSorted(obj)
            cueuncue = obj.CueUncue;

            [~, ~, indrmv] = rmoutliers_custome(obj.ChoiceTime);

            movement_time_sorted = cell(length(cueuncue), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    ind_thisFP = obj.Cued==obj.CueUncue(i) & obj.Stage==1 & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.ChoiceTime(setdiff(find(ind_thisFP), indrmv));
                    movement_time_sorted{i, j} = iMTs;
                end
            end
            value = movement_time_sorted;
        end

        function value = get.ChoiceTimePDF(obj)
            cueuncue = obj.CueUncue;

            MT_PDF = cell(length(cueuncue), length(obj.Ports));
            binEdges = obj.Bins.ChoiceTime;
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iMTs = obj.ChoiceTimeSorted{i, j};
                    if length(iMTs) > 5
                        MT_PDF{i, j} = ksdensity(iMTs, binEdges, 'Function', 'pdf');
                    else
                        MT_PDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = MT_PDF;
        end

        function value = get.ChoiceTimeCDF(obj)
            cueuncue = obj.CueUncue;

            MT_CDF = cell(length(cueuncue), length(obj.Ports));
            binEdges = obj.Bins.ChoiceTime;
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.CueUncue)
                    iMTs = obj.ChoiceTimeSorted{i, j};
                    if length(iMTs) > 5
                        MT_CDF{i, j} = ksdensity(iMTs, binEdges, 'Function', 'cdf');
                    else
                        MT_CDF{i, j} = zeros(1, length(binEdges));
                    end
                end
            end
            value = MT_CDF;
        end

        function value = get.ChoiceTimeStat(obj)
            datain = obj.ChoiceTime(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            thisFP = [];
            thisCued = [];
            N = [];
            Port = [];
            Median = [];
            Median_kde = [];
            Q1 = [];
            Q3 = [];

            % Only include correct trials
            for j = 1:2 % Port
                for i = 1:length(obj.CueUncue) % FP
                    iCTs = obj.ChoiceTimeSorted{i, j};
                    iCTs(isnan(iCTs)) = [];
                    % removeOutliers
                    iCTs(iCTs>data2575(2)+interq*c | iCTs<data2575(1)-interq*c) = [];
                    
                    thisFP = [thisFP; unique(obj.FP)];
                    thisCued = [thisCued; obj.CueUncue(i)];
                    N = [N; length(iCTs)];
                    Port = [Port; obj.Ports(j)];
                    iCTOut = calDur(iCTs*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);
                    Median = [Median; iCTOut.median*0.001];
                    Median_kde = [Median_kde; iCTOut.median_ksdensity*0.001];
                    Q1 = [Q1; prctile(iCTs, 25)];
                    Q3 = [Q3; prctile(iCTs, 75)];
                end
            end
            IQR = Q3 - Q1;

            Subjects = repmat(string(obj.Subject), length(thisFP), 1);
            Sessions = repmat(string(obj.Session), length(thisFP), 1);

            value = table(Subjects, Sessions, thisFP, thisCued, Port, N, Median, Median_kde, IQR, Q1, Q3);
        end

        function value = get.ChoiceTimeDistribution(obj)
            datain = obj.ChoiceTime(obj.Ind.correct & obj.Stage==1);
            datain(isnan(datain)) = [];
            [data2575] = prctile(datain, [25, 75]);
            interq = data2575(2) - data2575(1);
            c = 5;

            binEdges = 0:0.05:3;
            binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

            % Only include correct trials
            MT_PortL = obj.ChoiceTime(obj.Stage==1 & obj.Ind.correctL);
            MT_PortL(MT_PortL>data2575(2)+interq*c | MT_PortL<data2575(1)-interq*c) = [];
            MT_PortR = obj.ChoiceTime(obj.Stage==1 & obj.Ind.correctR);
            MT_PortR(MT_PortR>data2575(2)+interq*c | MT_PortR<data2575(1)-interq*c) = [];

            % compute RT from two ports separately
            Ncounts = histcounts(MT_PortL, binEdges);
            Counts_PortL = Ncounts';
            Ncounts = histcounts(MT_PortR, binEdges);
            Counts_PortR = Ncounts';
            MT_Centers = binCenters';

            switch obj.Task
                case {'Wait1Hold', 'Wait1HoldCRT', 'Wait2HoldCRT'}
                    movement_time_distr = table(MT_Centers, Counts_PortL, Counts_PortR);
                case {'ThreeFPHoldCRT', 'ThreeFPHoldSRT'}
                    FP_types = ["S", "M", "L"];

                    % compute RT from two ports and three FPs separately
                    for j = 1:length(obj.Ports)
                        for i = 1:length(obj.MixedFP)
                            iMTs = obj.ChoiceTimeSorted{i, j};
                            iMTs(isnan(iMTs)) = [];
                            % removeOutliers
                            iMTs(iMTs>data2575(2)+interq*c | iMTs<data2575(1)-interq*c) = [];
                            Ncounts = histcounts(iMTs, binEdges);
                            Ncounts = Ncounts';
                            eval("Counts_" + FP_types(i) + "_Port" + obj.Ports(j) + "= Ncounts;");
                        end
                    end
                    movement_time_distr = table(MT_Centers, Counts_PortL, Counts_PortR, Counts_S_PortL, Counts_M_PortL, Counts_L_PortL, ...
                        Counts_S_PortR, Counts_M_PortR, Counts_L_PortR);
            end
            value = movement_time_distr;
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
    end
end