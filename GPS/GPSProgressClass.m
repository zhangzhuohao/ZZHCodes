classdef GPSProgressClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Subject
        Strain
        Sessions
        Task
        NumSessions
        NumTrialsPerSession
        
        Treatment
        Dose
        Label

        MixedFP
        Ports

        BehavTable
        Performance
        PerformanceTrack

        STStat
        RTStat
        MTStat
        HDStat
        CTStat

        RTSorted
        MTSorted
        HDSorted
        CTSorted

        RTPDF
        MTPDF
        HDPDF
        CTPDF

        RTCDF
        MTCDF
        HDCDF
        CTCDF

        Interruption
    end

    properties (Constant)
        PhaseCount = 50; % Trial number to determine the early or late training phase
    end

    properties (Dependent)
        Ind
        Bins

        PerformanceAll
        PerformanceControl
        PerformanceChemo

        RTStatAll
        RTStatControl
        RTStatChemo

        RTSortedAll
        RTSortedNone
        RTSortedSaline
        RTSortedDCZ
        RTSortedControl
        RTSortedChemo

        MTStatAll
        MTStatControl
        MTStatChemo

        MTSortedAll
        MTSortedNone
        MTSortedSaline
        MTSortedDCZ
        MTSortedControl
        MTSortedChemo

        CTStatControl
        CTStatChemo

        CTSortedAll
        CTSortedNone
        CTSortedSaline
        CTSortedDCZ
        CTSortedControl
        CTSortedChemo
        
        HDStatAll
        HDStatControl
        HDStatChemo

        HDSortedAll
        HDSortedNone
        HDSortedSaline
        HDSortedDCZ
        HDSortedControl
        HDSortedChemo

        InterruptionNone
        InterruptionSaline
        InterruptionDCZ
        InterruptionControl
        InterruptionChemo
    end

    methods
        function obj = GPSProgressClass(SessionClassAll)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Subject             = string(unique(cellfun(@(x) x.Subject, SessionClassAll, 'UniformOutput', false)));
            obj.Strain              = string(unique(cellfun(@(x) char(x.Strain), SessionClassAll, 'UniformOutput', false)));
            obj.Sessions            = string(cellfun(@(x) x.Session, SessionClassAll, 'UniformOutput', false));
            obj.Task                = string(unique(cellfun(@(x) x.Task, SessionClassAll, 'UniformOutput', false)));
            obj.NumSessions         = length(SessionClassAll);
            obj.NumTrialsPerSession = cellfun(@(x) x.NumTrials, SessionClassAll);

            obj.Treatment           = string(cellfun(@(x) x.Treatment, SessionClassAll, 'UniformOutput', false));
            obj.Dose                = cellfun(@(x) x.Dose, SessionClassAll);
            obj.Label               = string(cellfun(@(x) x.Label, SessionClassAll, 'UniformOutput', false));

            MixedFP                 = cell2mat(cellfun(@(x) x.MixedFP, SessionClassAll , 'UniformOutput', false)');
            obj.MixedFP             = [unique(MixedFP(:, 1)) unique(MixedFP(:, 2)) unique(MixedFP(:, 3))];

            obj.Ports = ["L", "R"];

            % Behavior table & Performance tracking
            allTables               = cellfun(@(x)x.BehavTable, SessionClassAll, 'UniformOutput', false);
            for i = 1:obj.NumSessions
                thisTable = allTables{i};
                if i == 1
                    TrialStartTimeProgress  = thisTable.TrialStartTime;
                    TrialProgress           = thisTable.Trials;

                    thisTable               = addvars(thisTable, TrialStartTimeProgress, 'After', "TrialStartTime");
                    thisTable               = addvars(thisTable, TrialProgress, 'After', "Trials");
                    thisTable               = addvars(thisTable, repmat(obj.Treatment(i), height(thisTable), 1), 'After', "SessionDate", 'NewVariableNames', "Treatment");
                    thisTable               = addvars(thisTable, repmat(obj.Dose(i)     , height(thisTable), 1), 'After', "Treatment"  , 'NewVariableNames', "Dose");
                    thisTable               = addvars(thisTable, repmat(obj.Label(i)    , height(thisTable), 1), 'After', "Dose"       , 'NewVariableNames', "Label");

                    obj.BehavTable          = thisTable;
                else
                    TrialStartTimeProgress  = thisTable.TrialStartTime + obj.BehavTable.TrialStartTimeProgress(end);
                    TrialProgress           = thisTable.Trials + obj.BehavTable.TrialProgress(end);

                    thisTable               = addvars(thisTable, TrialStartTimeProgress, 'After', "TrialStartTime");
                    thisTable               = addvars(thisTable, TrialProgress, 'After', "Trials");
                    thisTable               = addvars(thisTable, repmat(obj.Treatment(i), height(thisTable), 1), 'After', "SessionDate", 'NewVariableNames', "Treatment");
                    thisTable               = addvars(thisTable, repmat(obj.Dose(i)     , height(thisTable), 1), 'After', "Treatment"  , 'NewVariableNames', "Dose");
                    thisTable               = addvars(thisTable, repmat(obj.Label(i)    , height(thisTable), 1), 'After', "Dose"       , 'NewVariableNames', "Label");

                    obj.BehavTable          = [obj.BehavTable; thisTable];
                end
            end

            % Performance
            allPerfs                = cellfun(@(x) x.Performance, SessionClassAll, 'UniformOutput', false);
            obj.Performance         = [];
            for i = 1:obj.NumSessions
                obj.Performance = [obj.Performance; allPerfs{i}];
            end

            % Performance tracking
            allPerfTracks           = cellfun(@(x) x.PerformanceTrack, SessionClassAll, 'UniformOutput', false);
            for i = 1:obj.NumSessions
                thisPerfTrack = allPerfTracks{i};
                if i == 1
                    WinPosProgress          = thisPerfTrack.WinPos;
                    thisPerfTrack           = addvars(thisPerfTrack, WinPosProgress, 'After', "WinPos");
                    obj.PerformanceTrack    = thisPerfTrack;
                else
                    WinPosProgress          = thisPerfTrack.WinPos + obj.BehavTable.TrialStartTimeProgress(find(obj.BehavTable.SessionDate==obj.Sessions(i-1), 1, 'last'));
                    thisPerfTrack           = addvars(thisPerfTrack, WinPosProgress, 'After', "WinPos");
                    obj.PerformanceTrack    = [obj.PerformanceTrack; thisPerfTrack];
                end
            end

            % Statistics
            allRTStats              = cellfun(@(x) x.RTStat, SessionClassAll, 'UniformOutput', false);
            obj.RTStat              = [];
            for i = 1:length(allTables)
                obj.RTStat = [obj.RTStat; allRTStats{i}];
            end

            allMTStats              = cellfun(@(x) x.MovementTimeStat, SessionClassAll, 'UniformOutput', false);
            obj.MTStat              = [];
            for i = 1:length(allTables)
                obj.MTStat = [obj.MTStat; allMTStats{i}];
            end

            allCTStats              = cellfun(@(x) x.ChoiceTimeStat, SessionClassAll, 'UniformOutput', false);
            obj.CTStat              = [];
            for i = 1:length(allTables)
                obj.CTStat = [obj.CTStat; allCTStats{i}];
            end

            allSTStats              = cellfun(@(x) x.ShuttleTimeStat, SessionClassAll, 'UniformOutput', false);
            obj.STStat              = [];
            for i = 1:length(allTables)
                obj.STStat = [obj.STStat; allSTStats{i}];
            end

            allHDStats              = cellfun(@(x) x.HoldDurationStat, SessionClassAll, 'UniformOutput', false);
            obj.HDStat              = [];
            for i = 1:length(allTables)
                obj.HDStat = [obj.HDStat; allHDStats{i}];
            end

            obj.RTSorted            = cellfun(@(x) x.RTSorted, SessionClassAll, 'UniformOutput', false);
            obj.HDSorted            = cellfun(@(x) x.HoldDurationSorted, SessionClassAll, 'UniformOutput', false);
            obj.MTSorted            = cellfun(@(x) x.MovementTimeSorted, SessionClassAll, 'UniformOutput', false);
            obj.CTSorted            = cellfun(@(x) x.ChoiceTimeSorted, SessionClassAll, 'UniformOutput', false);

            obj.RTPDF               = cellfun(@(x) x.RTPDF, SessionClassAll, 'UniformOutput', false);
            obj.HDPDF               = cellfun(@(x) x.HoldDurationPDF, SessionClassAll, 'UniformOutput', false);
            obj.MTPDF               = cellfun(@(x) x.MovementTimePDF, SessionClassAll, 'UniformOutput', false);
            obj.CTPDF               = cellfun(@(x) x.ChoiceTimePDF, SessionClassAll, 'UniformOutput', false);

            obj.RTCDF               = cellfun(@(x) x.RTCDF, SessionClassAll, 'UniformOutput', false);
            obj.HDCDF               = cellfun(@(x) x.HoldDurationCDF, SessionClassAll, 'UniformOutput', false);
            obj.MTCDF               = cellfun(@(x) x.MovementTimeCDF, SessionClassAll, 'UniformOutput', false);
            obj.CTCDF               = cellfun(@(x) x.ChoiceTimeCDF, SessionClassAll, 'UniformOutput', false);

            allInterruption         = cellfun(@(x) x.Interruption, SessionClassAll, 'UniformOutput', false);
            obj.Interruption        = [];
            for i = 1:obj.NumSessions
                thisTable = allInterruption{i};
                if i == 1
                    TrialProgress           = thisTable.Trials;
                    thisTable               = addvars(thisTable, TrialProgress, 'After', "Trials");
                    obj.Interruption        = thisTable;
                else
                    TrialProgress           = thisTable.Trials + obj.BehavTable.TrialProgress(find(obj.BehavTable.SessionDate==obj.Sessions(i), 1, "last"));
                    thisTable               = addvars(thisTable, TrialProgress, 'After', "Trials");
                    obj.Interruption        = [obj.Interruption; thisTable];
                end
            end
        end

        %%
        function value = get.Ind(obj)

            ind = feval(obj.Task + ".getProgressInd", obj);
            value = ind;
        end

        function value = get.Bins(obj)

            bins = feval(obj.Task + ".getBins");
            value = bins;
        end

        function value = get.PerformanceAll(obj)
            thisMixedFP = obj.MixedFP;
            if size(thisMixedFP, 1)==1
                thisMixedFP = thisMixedFP';
            end
            TargetPort = [repmat("L", length(obj.MixedFP)+1, 1); repmat("R", length(obj.MixedFP)+1, 1); "Both"]; % Left or Right or Both
            Foreperiod = [repmat([thisMixedFP; 0], 2, 1); 0];
            NumTrialsSorted = zeros(length(Foreperiod), 1);

            % group by foreperiod
            CorrectRatio    = zeros(length(Foreperiod), 1);
            PrematureRatio  = zeros(length(Foreperiod), 1);
            LateRatio       = zeros(length(Foreperiod), 1);
            WrongRatio      = zeros(length(Foreperiod), 1);

            % group by target port

            for j = 1:length(obj.Ports) % Target port
                for i = 1:length(obj.MixedFP)+1
                    ind = i + (j - 1) * (length(obj.MixedFP) + 1);
                    if i <= length(obj.MixedFP)
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==obj.MixedFP(i);
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    else
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==0;
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    end
                end
            end

            ind = length(Foreperiod);
            ind_this = obj.Performance.TargetPort=="Both" & obj.Performance.Foreperiod==0;
                        
            NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
            CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);

            Subjects = repmat(string(obj.Subject), length(Foreperiod), 1);

            perf_table = table(Subjects, Foreperiod, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = perf_table;
        end
        
        function value = get.PerformanceControl(obj)
            thisMixedFP = obj.MixedFP;
            if size(thisMixedFP, 1)==1
                thisMixedFP = thisMixedFP';
            end
            TargetPort = [repmat("L", length(obj.MixedFP)+1, 1); repmat("R", length(obj.MixedFP)+1, 1); "Both"]; % Left or Right or Both
            Foreperiod = [repmat([thisMixedFP; 0], 2, 1); 0];
            NumTrialsSorted = zeros(length(Foreperiod), 1);

            % group by foreperiod
            CorrectRatio    = zeros(length(Foreperiod), 1);
            PrematureRatio  = zeros(length(Foreperiod), 1);
            LateRatio       = zeros(length(Foreperiod), 1);
            WrongRatio      = zeros(length(Foreperiod), 1);

            % group by target port

            for j = 1:length(obj.Ports) % Target port
                for i = 1:length(obj.MixedFP)+1
                    ind = i + (j - 1) * (length(obj.MixedFP) + 1);
                    if i <= length(obj.MixedFP)
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==obj.MixedFP(i) & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Control"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    else
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==0 & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Control"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    end
                end
            end

            ind = length(Foreperiod);
            ind_this = obj.Performance.TargetPort=="Both" & obj.Performance.Foreperiod==0 & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Control"));
                        
            NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
            CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);

            Subjects = repmat(string(obj.Subject), length(Foreperiod), 1);

            perf_table = table(Subjects, Foreperiod, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = perf_table;
        end

        function value = get.PerformanceChemo(obj)
            thisMixedFP = obj.MixedFP;
            if size(thisMixedFP, 1)==1
                thisMixedFP = thisMixedFP';
            end
            TargetPort = [repmat("L", length(obj.MixedFP)+1, 1); repmat("R", length(obj.MixedFP)+1, 1); "Both"]; % Left or Right or Both
            Foreperiod = [repmat([thisMixedFP; 0], 2, 1); 0];
            NumTrialsSorted = zeros(length(Foreperiod), 1);

            % group by foreperiod
            CorrectRatio    = zeros(length(Foreperiod), 1);
            PrematureRatio  = zeros(length(Foreperiod), 1);
            LateRatio       = zeros(length(Foreperiod), 1);
            WrongRatio      = zeros(length(Foreperiod), 1);

            % group by target port

            for j = 1:length(obj.Ports) % Target port
                for i = 1:length(obj.MixedFP)+1
                    ind = i + (j - 1) * (length(obj.MixedFP) + 1);
                    if i <= length(obj.MixedFP)
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==obj.MixedFP(i) & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Chemo"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    else
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==0 & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Chemo"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    end
                end
            end

            ind = length(Foreperiod);
            ind_this = obj.Performance.TargetPort=="Both" & obj.Performance.Foreperiod==0 & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Chemo"));
                        
            NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this));
            CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);

            Subjects = repmat(string(obj.Subject), length(Foreperiod), 1);

            perf_table = table(Subjects, Foreperiod, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = perf_table;
        end

        %%
        function value = get.RTStatAll(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.RTStat.thisFP==thisFP(i) & obj.RTStat.Port==Port(i);
                N(i)          = mean(obj.RTStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.RTStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.RTStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.RTStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.RTStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.RTStat.Q3(ind_this), 'omitnan');
                Median_sem(i) = std(obj.RTStat.Median(ind_this), 'omitnan') / sqrt(sum(ind_this, 'omitnan'));
            end

            rt_stat_all = table(Subjects, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = rt_stat_all;
        end

        function value = get.RTStatControl(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Control", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.RTStat.thisFP==thisFP(i) & obj.RTStat.Port==Port(i) & ismember(obj.RTStat.Sessions, obj.Sessions(obj.Label=="Control"));
                N(i)          = mean(obj.RTStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.RTStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.RTStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.RTStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.RTStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.RTStat.Q3(ind_this), 'omitnan');
                Median_sem(i) = std(obj.RTStat.Median(ind_this), 'omitnan') / sqrt(sum(~isnan(Median)));
            end

            rt_stat_control = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = rt_stat_control;
        end

        function value = get.RTStatChemo(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Chemo", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.RTStat.thisFP==thisFP(i) & obj.RTStat.Port==Port(i) & ismember(obj.RTStat.Sessions, obj.Sessions(obj.Label=="Chemo"));
                N(i)          = mean(obj.RTStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.RTStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.RTStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.RTStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.RTStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.RTStat.Q3(ind_this), 'omitnan');
                Median_sem(i) = std(obj.RTStat.Median(ind_this), 'omitnan') / sqrt(sum(~isnan(Median)));
            end

            rt_stat_control = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = rt_stat_control;
        end

        function value = get.RTSortedAll(obj)

            rt_sorted_all = cell(length(obj.MixedFP), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.BehavTable.RT(find(ind_thisFP));
                    rt_sorted_all{i, j} = iRTs;
                end
            end

            value = rt_sorted_all;
        end

        function value = get.RTSortedNone(obj)

            rt_sorted_none = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="None"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.BehavTable.RT(find(ind_thisFP));
                    rt_sorted_none{i, j} = iRTs;
                end
            end

            value = rt_sorted_none;
        end

        function value = get.RTSortedSaline(obj)

            rt_sorted_saline = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind     = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="Saline"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.BehavTable.RT(find(ind_thisFP));
                    rt_sorted_saline{i, j} = iRTs;
                end
            end

            value = rt_sorted_saline;
        end

        function value = get.RTSortedDCZ(obj)

            rt_sorted_dcz = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind  = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="DCZ"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.BehavTable.RT(find(ind_thisFP));
                    rt_sorted_dcz{i, j} = iRTs;
                end
            end

            value = rt_sorted_dcz;
        end

        function value = get.RTSortedControl(obj)

            rt_sorted_Control = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind  = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Control"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.BehavTable.RT(find(ind_thisFP));
                    rt_sorted_Control{i, j} = iRTs;
                end
            end

            value = rt_sorted_Control;
        end

        function value = get.RTSortedChemo(obj)

            rt_sorted_chemo = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind  = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Chemo"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iRTs = obj.BehavTable.RT(find(ind_thisFP));
                    rt_sorted_chemo{i, j} = iRTs;
                end
            end

            value = rt_sorted_chemo;
        end

        %%
        function value = get.HDStatAll(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.HDStat.thisFP==thisFP(i) & obj.HDStat.Port==Port(i);
                N(i)          = mean(obj.HDStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.HDStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.HDStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.HDStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.HDStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.HDStat.Q3(ind_this), 'omitnan');
            end

            hd_stat_all = table(Subjects, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3);
            value = hd_stat_all;
        end

        function value = get.HDStatControl(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Control", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.HDStat.thisFP==thisFP(i) & obj.HDStat.Port==Port(i) & ismember(obj.HDStat.Sessions, obj.Sessions(obj.Label=="Control"));
                N(i)          = mean(obj.HDStat.N(ind_this));
                Median(i)     = mean(obj.HDStat.Median(ind_this));
                Median_kde(i) = mean(obj.HDStat.Median_kde(ind_this));
                IQR(i)        = mean(obj.HDStat.IQR(ind_this));
                Q1(i)         = mean(obj.HDStat.Q1(ind_this));
                Q3(i)         = mean(obj.HDStat.Q3(ind_this));
            end

            hd_stat_control = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3);
            value = hd_stat_control;
        end

        function value = get.HDStatChemo(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Chemo", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.HDStat.thisFP==thisFP(i) & obj.HDStat.Port==Port(i) & ismember(obj.HDStat.Sessions, obj.Sessions(obj.Label=="Chemo"));
                N(i)          = mean(obj.HDStat.N(ind_this));
                Median(i)     = mean(obj.HDStat.Median(ind_this));
                Median_kde(i) = mean(obj.HDStat.Median_kde(ind_this));
                IQR(i)        = mean(obj.HDStat.IQR(ind_this));
                Q1(i)         = mean(obj.HDStat.Q1(ind_this));
                Q3(i)         = mean(obj.HDStat.Q3(ind_this));
            end

            hd_stat_chemo = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3);
            value = hd_stat_chemo;
        end

        function value = get.HDSortedAll(obj)

            hd_sorted_all = cell(length(obj.MixedFP), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.port" + obj.Ports(j));
                    iHDs = obj.BehavTable.HoldDuration(find(ind_thisFP));
                    hd_sorted_all{i, j} = iHDs;
                end
            end

            value = hd_sorted_all;
        end

        function value = get.HDSortedNone(obj)

            hd_sorted_none = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="None"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.port" + obj.Ports(j));
                    iHDs = obj.BehavTable.HoldDuration(find(ind_thisFP));
                    hd_sorted_none{i, j} = iHDs;
                end
            end

            value = hd_sorted_none;
        end

        function value = get.HDSortedSaline(obj)

            hd_sorted_saline = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind     = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="Saline"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.port" + obj.Ports(j));
                    iHDs = obj.BehavTable.HoldDuration(find(ind_thisFP));
                    hd_sorted_saline{i, j} = iHDs;
                end
            end

            value = hd_sorted_saline;
        end

        function value = get.HDSortedDCZ(obj)

            hd_sorted_dcz = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind  = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="DCZ"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.port" + obj.Ports(j));
                    iHDs = obj.BehavTable.HoldDuration(find(ind_thisFP));
                    hd_sorted_dcz{i, j} = iHDs;
                end
            end

            value = hd_sorted_dcz;
        end

        function value = get.HDSortedControl(obj)

            hd_sorted_Control = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Control"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.port" + obj.Ports(j));
                    iHDs = obj.BehavTable.HoldDuration(find(ind_thisFP));
                    hd_sorted_Control{i, j} = iHDs;
                end
            end

            value = hd_sorted_Control;
        end

        function value = get.HDSortedChemo(obj)

            hd_sorted_chemo = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind    = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Chemo"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.port" + obj.Ports(j));
                    iHDs = obj.BehavTable.HoldDuration(find(ind_thisFP));
                    hd_sorted_chemo{i, j} = iHDs;
                end
            end

            value = hd_sorted_chemo;
        end

        %%
        function value = get.MTStatAll(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.MTStat.thisFP==thisFP(i) & obj.MTStat.Port==Port(i);
                N(i)          = mean(obj.MTStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.MTStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.MTStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.MTStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.MTStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.MTStat.Q3(ind_this), 'omitnan');
                Median_sem(i) = std(obj.MTStat.Median(ind_this), 'omitnan') / sqrt(sum(~isnan(Median)));
            end

            mt_stat_all = table(Subjects, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = mt_stat_all;
        end

        function value = get.MTStatControl(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Control", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.MTStat.thisFP==thisFP(i) & obj.MTStat.Port==Port(i) & ismember(obj.MTStat.Sessions, obj.Sessions(obj.Label=="Control"));
                N(i)          = mean(obj.MTStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.MTStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.MTStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.MTStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.MTStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.MTStat.Q3(ind_this), 'omitnan');
                Median_sem(i) = std(obj.MTStat.Median(ind_this), 'omitnan') / sqrt(sum(~isnan(Median)));
            end

            mt_stat_control = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = mt_stat_control;
        end

        function value = get.MTStatChemo(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Chemo", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.MTStat.thisFP==thisFP(i) & obj.MTStat.Port==Port(i) & ismember(obj.MTStat.Sessions, obj.Sessions(obj.Label=="Chemo"));
                N(i)          = mean(obj.MTStat.N(ind_this), 'omitnan');
                Median(i)     = mean(obj.MTStat.Median(ind_this), 'omitnan');
                Median_kde(i) = mean(obj.MTStat.Median_kde(ind_this), 'omitnan');
                IQR(i)        = mean(obj.MTStat.IQR(ind_this), 'omitnan');
                Q1(i)         = mean(obj.MTStat.Q1(ind_this), 'omitnan');
                Q3(i)         = mean(obj.MTStat.Q3(ind_this), 'omitnan');
                Median_sem(i) = std(obj.MTStat.Median(ind_this), 'omitnan') / sqrt(sum(~isnan(Median)));
            end

            mt_stat_chemo = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = mt_stat_chemo;
        end

        function value = get.MTSortedAll(obj)

            mt_sorted_all = cell(length(obj.MixedFP), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.BehavTable.MovementTime(find(ind_thisFP));
                    mt_sorted_all{i, j} = iMTs;
                end
            end

            value = mt_sorted_all;
        end

        function value = get.MTSortedNone(obj)

            mt_sorted_none = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="None"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.BehavTable.MovementTime(find(ind_thisFP));
                    mt_sorted_none{i, j} = iMTs;
                end
            end

            value = mt_sorted_none;
        end

        function value = get.MTSortedSaline(obj)

            mt_sorted_saline = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind     = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="Saline"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.BehavTable.MovementTime(find(ind_thisFP));
                    mt_sorted_saline{i, j} = iMTs;
                end
            end

            value = mt_sorted_saline;
        end

        function value = get.MTSortedDCZ(obj)

            mt_sorted_dcz = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind  = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="DCZ"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.BehavTable.MovementTime(find(ind_thisFP));
                    mt_sorted_dcz{i, j} = iMTs;
                end
            end

            value = mt_sorted_dcz;
        end

        function value = get.MTSortedControl(obj)

            mt_sorted_Control = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Control"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.BehavTable.MovementTime(find(ind_thisFP));
                    mt_sorted_Control{i, j} = iMTs;
                end
            end

            value = mt_sorted_Control;
        end

        function value = get.MTSortedChemo(obj)

            mt_sorted_chemo = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Chemo"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iMTs = obj.BehavTable.MovementTime(find(ind_thisFP));
                    mt_sorted_chemo{i, j} = iMTs;
                end
            end

            value = mt_sorted_chemo;
        end

        %%
        function value = get.CTStatControl(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Control", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.CTStat.thisFP==thisFP(i) & obj.CTStat.Port==Port(i) & ismember(obj.CTStat.Sessions, obj.Sessions(obj.Label=="Control"));
                N(i)          = mean(obj.CTStat.N(ind_this));
                Median(i)     = mean(obj.CTStat.Median(ind_this));
                Median_kde(i) = mean(obj.CTStat.Median_kde(ind_this));
                IQR(i)        = mean(obj.CTStat.IQR(ind_this));
                Q1(i)         = mean(obj.CTStat.Q1(ind_this));
                Q3(i)         = mean(obj.CTStat.Q3(ind_this));
                Median_sem(i) = std(obj.CTStat.Median(ind_this)) / sqrt(sum(ind_this));
            end

            mt_stat_control = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = mt_stat_control;
        end

        function value = get.CTStatChemo(obj)
            
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.MixedFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.MixedFP), 1), [], 1)];

            num_entry = length(thisFP);
            Subjects  = repmat(obj.Subject, num_entry, 1);
            Labels    = repmat("Chemo", num_entry, 1);
            
            N          = zeros(num_entry, 1);
            Median     = zeros(num_entry, 1);
            Median_kde = zeros(num_entry, 1);
            IQR        = zeros(num_entry, 1);
            Q1         = zeros(num_entry, 1);
            Q3         = zeros(num_entry, 1);
            Median_sem = zeros(num_entry, 1);
            for i = 1:num_entry
                ind_this   = obj.CTStat.thisFP==thisFP(i) & obj.CTStat.Port==Port(i) & ismember(obj.CTStat.Sessions, obj.Sessions(obj.Label=="Chemo"));
                N(i)          = mean(obj.CTStat.N(ind_this));
                Median(i)     = mean(obj.CTStat.Median(ind_this));
                Median_kde(i) = mean(obj.CTStat.Median_kde(ind_this));
                IQR(i)        = mean(obj.CTStat.IQR(ind_this));
                Q1(i)         = mean(obj.CTStat.Q1(ind_this));
                Q3(i)         = mean(obj.CTStat.Q3(ind_this));
                Median_sem(i) = std(obj.CTStat.Median(ind_this)) / sqrt(sum(ind_this));
            end

            mt_stat_chemo = table(Subjects, Labels, thisFP, Port, N, Median, Median_kde, IQR, Q1, Q3, Median_sem);
            value = mt_stat_chemo;
        end

        function value = get.CTSortedAll(obj)

            ct_sorted_all = cell(length(obj.MixedFP), length(obj.Ports));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iCTs = obj.BehavTable.ChoiceTime(find(ind_thisFP));
                    ct_sorted_all{i, j} = iCTs;
                end
            end

            value = ct_sorted_all;
        end

        function value = get.CTSortedNone(obj)

            mt_sorted_none = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="None"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iCTs = obj.BehavTable.ChoiceTime(find(ind_thisFP));
                    mt_sorted_none{i, j} = iCTs;
                end
            end

            value = mt_sorted_none;
        end

        function value = get.CTSortedSaline(obj)

            mt_sorted_saline = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind     = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="Saline"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iCTs = obj.BehavTable.ChoiceTime(find(ind_thisFP));
                    mt_sorted_saline{i, j} = iCTs;
                end
            end

            value = mt_sorted_saline;
        end

        function value = get.CTSortedDCZ(obj)

            mt_sorted_dcz = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind  = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Treatment=="DCZ"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iCTs = obj.BehavTable.ChoiceTime(find(ind_thisFP));
                    mt_sorted_dcz{i, j} = iCTs;
                end
            end

            value = mt_sorted_dcz;
        end

        function value = get.CTSortedControl(obj)

            mt_sorted_Control = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Control"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iCTs = obj.BehavTable.ChoiceTime(find(ind_thisFP));
                    mt_sorted_Control{i, j} = iCTs;
                end
            end

            value = mt_sorted_Control;
        end

        %%
        function value = get.CTSortedChemo(obj)

            mt_sorted_chemo = cell(length(obj.MixedFP), length(obj.Ports));
            sessions_ind   = ismember(obj.BehavTable.SessionDate, obj.Sessions(obj.Label=="Chemo"));
            for j = 1:length(obj.Ports)
                for i = 1:length(obj.MixedFP)
                    ind_thisFP = sessions_ind & obj.BehavTable.Stage==1 & obj.BehavTable.FP==obj.MixedFP(i) & eval("obj.Ind.correct" + obj.Ports(j));
                    iCTs = obj.BehavTable.ChoiceTime(find(ind_thisFP));
                    mt_sorted_chemo{i, j} = iCTs;
                end
            end

            value = mt_sorted_chemo;
        end

        %%
        function value = get.InterruptionControl(obj)

            sessions_ind = ismember(obj.Interruption.Sessions, obj.Sessions(obj.Label=="Control"));
            interruption_control = obj.Interruption(sessions_ind, :);
            
            value = interruption_control;
        end

        function value = get.InterruptionChemo(obj)

            sessions_ind = ismember(obj.Interruption.Sessions, obj.Sessions(obj.Label=="Chemo"));
            interruption_chemo = obj.Interruption(sessions_ind, :);
            
            value = interruption_chemo;
        end

        function value = get.InterruptionNone(obj)

            sessions_ind = ismember(obj.Interruption.Sessions, obj.Sessions(obj.Treatment=="None"));
            interruption_none = obj.Interruption(sessions_ind, :);
            
            value = interruption_none;
        end
        
        function value = get.InterruptionSaline(obj)

            sessions_ind = ismember(obj.Interruption.Sessions, obj.Sessions(obj.Treatment=="Saline"));
            interruption_saline = obj.Interruption(sessions_ind, :);
            
            value = interruption_saline;
        end

        function value = get.InterruptionDCZ(obj)

            sessions_ind = ismember(obj.Interruption.Sessions, obj.Sessions(obj.Treatment=="DCZ"));
            interruption_dcz = obj.Interruption(sessions_ind, :);
            
            value = interruption_dcz;
        end

        %%
        function save(obj, savePath)

            saveName = fullfile(savePath, "GPSProgressClass_" + obj.Task + "_" + upper(obj.Subject) + ".mat");
            save(saveName, 'obj');
        end

        function print(obj, Func, targetDir)
            
            if nargin==1
                Func = "Progress";
            end

            switch lower(Func)
                case {'progress'}
                    hf = obj.plotProgress();
                case {'chemoeffect'}
                    hf = obj.plotChemoEffect();
                case {'show'}
                    hf = obj.plotShow();
            end

            if nargin==3
                % check if targetDir exists
                if ~contains(targetDir, '/') && ~contains(targetDir, '\')
                    % so it is a relative path
                    if ~exist(targetDir, 'dir')
                        mkdir(targetDir)
                    end
                end
                savename = fullfile(targetDir, "GPSProgressClass_" + obj.Task + "_" + Func + "_" + upper(obj.Subject));
                print(hf, '-dpdf', savename, '-bestfit')
                print(hf, '-dpng', savename)
                saveas(hf, savename, 'fig')
                exportgraphics(hf, savename+".jpg", 'Resolution', 600);
            end
        end

        function fig = plotProgress(obj)

            try
                set_matlab_default
            catch
                disp('You do not have "set_matlab_default"');
            end

            fig = feval(obj.Task + ".plotProgress", obj);
        end

        function fig = plotChemoEffect(obj)

            try
                set_matlab_default
            catch
                disp('You do not have "set_matlab_default"');
            end

            fig = feval(obj.Task + ".plotChemoEffect", obj);
        end

        function fig = plotShow(obj)

            try
                set_matlab_default
            catch
                disp('You do not have "set_matlab_default"');
            end

            fig = feval(obj.Task + ".plotShow", obj);
        end

    end
end
