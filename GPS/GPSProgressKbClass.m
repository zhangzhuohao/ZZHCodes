classdef GPSProgressKbClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Subject
        Strain
        Sessions
        Task
        TaskName
        NumSessions
        NumTrialsPerSession
        
        Treatment
        Dose
        Label

        TargetFP
        Ports
        Cued
        Guided
        Filled

        BehavTable
        Performance
        PerformanceTrackCue
        PerformanceTrackUncue
        PerformanceTrackUncueL
        PerformanceTrackUncueR

        STStat
        RTStat
        MTStat
        HDStat
        CTStat

        RTSorted
        MTSorted
        HDSorted
        CTSorted
        STSplit

        RTPDF
        MTPDF
        HDPDF
        CTPDF
        LogSTPDF

        RTCDF
        MTCDF
        HDCDF
        CTCDF
        LogSTCDF

        Interruption
    end

    properties (Constant)
        PhaseCount = 50; % Trial number to determine the early or late training phase
        CueUncue = [1, 0];
        BandWidth = .05;
    end

    properties (Dependent)
        Ind
        Bins

        PerformanceAll
        PerformanceControl
        PerformanceChemo

        InterruptionNone
        InterruptionSaline
        InterruptionDCZ
        InterruptionControl
        InterruptionChemo
    end

    methods
        function obj = GPSProgressKbClass(SessionClassAll)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Subject             = string(unique(cellfun(@(x) x.Subject, SessionClassAll, 'UniformOutput', false)));
            obj.Strain              = string(unique(cellfun(@(x) char(x.Strain), SessionClassAll, 'UniformOutput', false)));
            obj.Sessions            = string(cellfun(@(x) x.Session, SessionClassAll, 'UniformOutput', false));
            obj.Task                = string(unique(cellfun(@(x) x.Task, SessionClassAll, 'UniformOutput', false)));
            obj.TaskName            = string(unique(cellfun(@(x) x.TaskName, SessionClassAll, 'UniformOutput', false)));
            obj.NumSessions         = length(SessionClassAll);
            obj.NumTrialsPerSession = cellfun(@(x) x.NumTrials, SessionClassAll);

            obj.Treatment           = string(cellfun(@(x) x.Treatment, SessionClassAll, 'UniformOutput', false));
            obj.Dose                = cellfun(@(x) x.Dose, SessionClassAll);
            obj.Label               = string(cellfun(@(x) x.Label, SessionClassAll, 'UniformOutput', false));
            
            obj.Cued                = cell2mat(cellfun(@(x) x.Cued, SessionClassAll , 'UniformOutput', false));
            obj.Guided              = cell2mat(cellfun(@(x) x.Guided, SessionClassAll , 'UniformOutput', false));
            obj.Filled              = cell2mat(cellfun(@(x) x.Filled, SessionClassAll , 'UniformOutput', false));

            TargetFP                = cell2mat(cellfun(@(x) x.TargetFP, SessionClassAll , 'UniformOutput', false)');
            obj.TargetFP            = unique(TargetFP(:, 1));

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
            allPerfTracksCue           = cellfun(@(x) x.PerformanceTrackCue, SessionClassAll, 'UniformOutput', false);
            allPerfTracksUncue         = cellfun(@(x) x.PerformanceTrackUncue, SessionClassAll, 'UniformOutput', false);
            allPerfTracksUncueL        = cellfun(@(x) x.PerformanceTrackUncueL, SessionClassAll, 'UniformOutput', false);
            allPerfTracksUncueR        = cellfun(@(x) x.PerformanceTrackUncueR, SessionClassAll, 'UniformOutput', false);
            for i = 1:obj.NumSessions
                thisPerfTrackCue = allPerfTracksCue{i};
                thisPerfTrackUncue = allPerfTracksUncue{i};
                thisPerfTrackUncueL = allPerfTracksUncueL{i};
                thisPerfTrackUncueR = allPerfTracksUncueR{i};
                if i == 1
                    WinPosProgressCue           = thisPerfTrackCue.WinPos;
                    WinPosProgressUncue         = thisPerfTrackUncue.WinPos;
                    WinPosProgressUncueL        = thisPerfTrackUncueL.WinPos;
                    WinPosProgressUncueR        = thisPerfTrackUncueR.WinPos;
                    thisPerfTrackCue            = addvars(thisPerfTrackCue, WinPosProgressCue, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    thisPerfTrackUncue          = addvars(thisPerfTrackUncue, WinPosProgressUncue, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    thisPerfTrackUncueL         = addvars(thisPerfTrackUncueL, WinPosProgressUncueL, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    thisPerfTrackUncueR         = addvars(thisPerfTrackUncueR, WinPosProgressUncueR, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    obj.PerformanceTrackCue     = thisPerfTrackCue;
                    obj.PerformanceTrackUncue   = thisPerfTrackUncue;
                    obj.PerformanceTrackUncueL  = thisPerfTrackUncueL;
                    obj.PerformanceTrackUncueR  = thisPerfTrackUncueR;
                else
                    WinEndCue    = obj.BehavTable.TrialStartTimeProgress(find(obj.BehavTable.SessionDate==obj.Sessions(i-1), 1, 'last'));
                    WinEndUncue  = obj.BehavTable.TrialStartTimeProgress(find(obj.BehavTable.SessionDate==obj.Sessions(i-1), 1, 'last'));
                    WinPosProgressCue           = thisPerfTrackCue.WinPos + max([WinEndCue WinEndUncue]);
                    WinPosProgressUncue         = thisPerfTrackUncue.WinPos + max([WinEndCue WinEndUncue]);
                    WinPosProgressUncueL        = thisPerfTrackUncueL.WinPos + max([WinEndCue WinEndUncue]);
                    WinPosProgressUncueR        = thisPerfTrackUncueR.WinPos + max([WinEndCue WinEndUncue]);
                    thisPerfTrackCue            = addvars(thisPerfTrackCue, WinPosProgressCue, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    thisPerfTrackUncue          = addvars(thisPerfTrackUncue, WinPosProgressUncue, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    thisPerfTrackUncueL         = addvars(thisPerfTrackUncueL, WinPosProgressUncueL, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    thisPerfTrackUncueR         = addvars(thisPerfTrackUncueR, WinPosProgressUncueR, 'After', "WinPos", 'NewVariableNames', "WinPosProgress");
                    obj.PerformanceTrackCue     = [obj.PerformanceTrackCue; thisPerfTrackCue];
                    obj.PerformanceTrackUncue   = [obj.PerformanceTrackUncue; thisPerfTrackUncue];
                    obj.PerformanceTrackUncueL  = [obj.PerformanceTrackUncueL; thisPerfTrackUncueL];
                    obj.PerformanceTrackUncueR  = [obj.PerformanceTrackUncueR; thisPerfTrackUncueR];
                end
            end

            % Statistics
            allRTStats              = cellfun(@(x) x.RTStat, SessionClassAll, 'UniformOutput', false);
            obj.RTStat.Session      = [];
            for i = 1:length(allTables)
                obj.RTStat.Session = [obj.RTStat.Session; allRTStats{i}];
            end

            allMTStats              = cellfun(@(x) x.MovementTimeStat, SessionClassAll, 'UniformOutput', false);
            obj.MTStat.Session      = [];
            for i = 1:length(allTables)
                obj.MTStat.Session = [obj.MTStat.Session; allMTStats{i}];
            end

            allCTStats              = cellfun(@(x) x.ChoiceTimeStat, SessionClassAll, 'UniformOutput', false);
            obj.CTStat.Session      = [];
            for i = 1:length(allTables)
                obj.CTStat.Session = [obj.CTStat.Session; allCTStats{i}];
            end

            allSTStats              = cellfun(@(x) x.ShuttleTimeStat, SessionClassAll, 'UniformOutput', false);
            obj.STStat.Session      = [];
            for i = 1:length(allTables)
                obj.STStat.Session = [obj.STStat.Session; allSTStats{i}];
            end

            allHDStats              = cellfun(@(x) x.HoldDurationStat, SessionClassAll, 'UniformOutput', false);
            obj.HDStat.Session      = [];
            for i = 1:length(allTables)
                obj.HDStat.Session = [obj.HDStat.Session; allHDStats{i}];
            end

            % Sorted data
            obj.RTSorted.Session    = cellfun(@(x) x.RTSorted, SessionClassAll, 'UniformOutput', false);
            obj.HDSorted.Session    = cellfun(@(x) x.HoldDurationSorted, SessionClassAll, 'UniformOutput', false);
            obj.MTSorted.Session    = cellfun(@(x) x.MovementTimeSorted, SessionClassAll, 'UniformOutput', false);
            obj.CTSorted.Session    = cellfun(@(x) x.ChoiceTimeSorted, SessionClassAll, 'UniformOutput', false);

            obj = obj.gatherAllSorted();

            % Splited shuttle time
            obj.STSplit.Session     = cellfun(@(x) x.ShuttleTimeSplit, SessionClassAll, 'UniformOutput', false);
            obj = obj.gatherAllSplit();

            % KDEs
            obj.RTPDF.Session       = cellfun(@(x) x.RTPDF, SessionClassAll, 'UniformOutput', false);
            obj.HDPDF.Session       = cellfun(@(x) x.HDPDF, SessionClassAll, 'UniformOutput', false);
            obj.MTPDF.Session       = cellfun(@(x) x.MTPDF, SessionClassAll, 'UniformOutput', false);
            obj.CTPDF.Session       = cellfun(@(x) x.CTPDF, SessionClassAll, 'UniformOutput', false);
            obj.LogSTPDF.Session    = cellfun(@(x) x.LogSTPDF, SessionClassAll, 'UniformOutput', false);

            obj.RTCDF.Session       = cellfun(@(x) x.RTCDF, SessionClassAll, 'UniformOutput', false);
            obj.HDCDF.Session       = cellfun(@(x) x.HDCDF, SessionClassAll, 'UniformOutput', false);
            obj.MTCDF.Session       = cellfun(@(x) x.MTCDF, SessionClassAll, 'UniformOutput', false);
            obj.CTCDF.Session       = cellfun(@(x) x.CTCDF, SessionClassAll, 'UniformOutput', false);
            obj.LogSTPDF.Session    = cellfun(@(x) x.LogSTPDF, SessionClassAll, 'UniformOutput', false);

            obj = obj.getAllKDEs(0);

            % Interruptions
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
            thisTargetFP = obj.TargetFP;
            if size(thisTargetFP, 1)==1
                thisTargetFP = thisTargetFP';
            end
            TargetPort = [repmat("L", length(obj.TargetFP)+1, 1); repmat("R", length(obj.TargetFP)+1, 1); "Both"]; % Left or Right or Both
            Foreperiod = [repmat([thisTargetFP; 0], 2, 1); 0];
            NumTrialsSorted = zeros(length(Foreperiod), 1);

            % group by foreperiod
            CorrectRatio    = zeros(length(Foreperiod), 1);
            PrematureRatio  = zeros(length(Foreperiod), 1);
            LateRatio       = zeros(length(Foreperiod), 1);
            WrongRatio      = zeros(length(Foreperiod), 1);

            % group by target port

            for j = 1:length(obj.Ports) % Target port
                for i = 1:length(obj.TargetFP)+1
                    ind = i + (j - 1) * (length(obj.TargetFP) + 1);
                    if i <= length(obj.TargetFP)
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Foreperiod==obj.TargetFP(i);
                        
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

            TargetPort = [repmat("L", length(obj.CueUncue)+1, 1); repmat("R", length(obj.CueUncue)+1, 1); "Both"]; % Left or Right or Both
            thisCued = [repmat([obj.CueUncue'; -1], 2, 1); -1];
            NumTrialsSorted = zeros(length(thisCued), 1);

            % group by foreperiod
            CorrectRatio    = zeros(length(thisCued), 1);
            PrematureRatio  = zeros(length(thisCued), 1);
            LateRatio       = zeros(length(thisCued), 1);
            WrongRatio      = zeros(length(thisCued), 1);

            % group by target port

            for j = 1:length(obj.Ports) % Target port
                for i = 1:length(obj.CueUncue)+1
                    ind = i + (j - 1) * (length(obj.CueUncue) + 1);
                    if i <= length(obj.CueUncue)
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Cued_this==obj.CueUncue(i) & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Control"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this), "omitnan");
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    else
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Control"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this), "omitnan");
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    end
                end
            end

            ind = length(thisCued);
            ind_this = obj.Performance.TargetPort=="Both" & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Control"));
                        
            NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this), "omitnan");
            CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);

            Subjects = repmat(string(obj.Subject), length(thisCued), 1);

            perf_table = table(Subjects, thisCued, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = perf_table;
        end

        function value = get.PerformanceChemo(obj)
            TargetPort = [repmat("L", length(obj.CueUncue)+1, 1); repmat("R", length(obj.CueUncue)+1, 1); "Both"]; % Left or Right or Both
            thisCued = [repmat([obj.CueUncue'; -1], 2, 1); -1];
            NumTrialsSorted = zeros(length(thisCued), 1);

            % group by foreperiod
            CorrectRatio    = zeros(length(thisCued), 1);
            PrematureRatio  = zeros(length(thisCued), 1);
            LateRatio       = zeros(length(thisCued), 1);
            WrongRatio      = zeros(length(thisCued), 1);

            % group by target port

            for j = 1:length(obj.Ports) % Target port
                for i = 1:length(obj.CueUncue)+1
                    ind = i + (j - 1) * (length(obj.CueUncue) + 1);
                    if i <= length(obj.CueUncue)
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & obj.Performance.Cued_this==obj.CueUncue(i) & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Chemo"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this), "omitnan");
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    else
                        ind_this = obj.Performance.TargetPort==obj.Ports(j) & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Chemo"));
                        
                        NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this), "omitnan");
                        CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                        WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
                    end
                end
            end

            ind = length(thisCued);
            ind_this = obj.Performance.TargetPort=="Both" & ismember(obj.Performance.Sessions, obj.Sessions(obj.Label=="Chemo"));
                        
            NumTrialsSorted(ind) = sum(obj.Performance.NumTrialsSorted(ind_this), "omitnan");
            CorrectRatio(ind)    = obj.Performance.CorrectRatio(ind_this)'   * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            PrematureRatio(ind)  = obj.Performance.PrematureRatio(ind_this)' * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            LateRatio(ind)       = obj.Performance.LateRatio(ind_this)'      * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);
            WrongRatio(ind)      = obj.Performance.WrongRatio(ind_this)'     * obj.Performance.NumTrialsSorted(ind_this) / NumTrialsSorted(ind);

            Subjects = repmat(string(obj.Subject), length(thisCued), 1);

            perf_table = table(Subjects, thisCued, TargetPort, NumTrialsSorted, CorrectRatio, PrematureRatio, LateRatio, WrongRatio);
            value = perf_table;
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
        function value = getPerformanceTrack(~, Perf, Ind)

            num_trials = height(Perf);

            WinSize  = floor(num_trials / 4);
            StepSize = max(1, floor(WinSize / 5));

            CountStart = 1;
            WinPos     = [];

            CorrectRatio    = [];
            WrongRatio      = [];
            PrematureRatio  = [];
            LateRatio       = [];

            center_pokes = Perf.TrialStartTime + Perf.CentInTime; % only count the first one
            while CountStart+WinSize-1 < num_trials

                thisWin = CountStart:(CountStart+WinSize-1);

                CorrectRatio    = [CorrectRatio    ;  100 * sum(Ind.correct(thisWin))   / WinSize];
                WrongRatio      = [WrongRatio      ;  100 * sum(Ind.wrong(thisWin))     / WinSize];
                PrematureRatio  = [PrematureRatio  ;  100 * sum(Ind.premature(thisWin)) / WinSize];
                LateRatio       = [LateRatio       ;  100 * sum(Ind.late(thisWin))      / WinSize];
                WinPos          = [WinPos          ;  center_pokes(thisWin(end))];
                
                CountStart = CountStart + StepSize;
            end

            perf_track = table(WinPos, CorrectRatio, WrongRatio, PrematureRatio, LateRatio);
            value = perf_track;
        end

        %%
        function save(obj, savePath)

            saveName = fullfile(savePath, "GPSProgressClass_" + obj.TaskName + "_" + upper(obj.Subject) + ".mat");
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
                savename = fullfile(targetDir, "GPSProgressClass_" + obj.TaskName + "_" + Func + "_" + upper(obj.Subject));
                print(hf, '-dpdf', savename, '-bestfit')
                print(hf, '-dpng', savename)
                saveas(hf, savename, 'fig')
                exportgraphics(hf, strcat(savename, ".jpg"), 'Resolution', 600);
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

                %% Data processing
        function obj = gatherAllSorted(obj)
            Vars = ["RT", "MT", "HD", "CT"];
            Labels = ["All", "Control", "Chemo"];

            for v = 1:length(Vars)
                for l = 1:length(Labels)
                    if Labels(l)=="All"
                        session_this = obj.Label~="Chemo";
                    else
                        session_this = obj.Label==Labels(l);
                    end
                    if any(session_this)
                        obj.(Vars(v)+"Sorted").(Labels(l)) = obj.spliceDataCell(obj.(Vars(v)+"Sorted").Session(session_this));
                    end
                end
            end
        end % gatherSorted

        function obj = gatherAllSplit(obj)
            Vars = "ST";
            Labels = ["All", "Control", "Chemo"];

            for v = 1:length(Vars)
                for l = 1:length(Labels)
                    if Labels(l)=="All"
                        session_this = obj.Label~="Chemo";
                    else
                        session_this = obj.Label==Labels(l);
                    end
                    if any(session_this)
                        obj.(Vars(v)+"Split").(Labels(l)) = obj.spliceDataCell(obj.(Vars(v)+"Split").Session(session_this));
                    end
                end
            end
        end

        function data_spliced = spliceDataCell(~, data_cell)
            sz = size(data_cell{1});
            data_spliced = cell(sz(1), sz(2));
            for i = 1:sz(1)
                for j = 1:sz(2)
                    data_this = cellfun(@(x) x{i, j}(:), data_cell, 'UniformOutput', false);
                    data_spliced{i, j} = cell2mat(data_this);
                end
            end
        end

        %%
        function stat = gatherStat(obj, var, lb)
            thisFP    = [zeros(length(obj.Ports), 1); repmat(obj.TargetFP', length(obj.Ports), 1)];
            Port      = [obj.Ports'; reshape(repmat(obj.Ports, length(obj.TargetFP), 1), [], 1)];
            num_entry = length(thisFP);
            
            stat_session = obj.(var+"Stat").Session;
            stat = array2table(zeros(num_entry, 10), "VariableNames", {'Subjects', 'thisFP', 'Port', 'N', 'Median', 'Median_kde', 'IQR', 'Q1', 'Q3', 'Median_sem'});
            stat.Subjects   = repmat(obj.Subject, num_entry, 1);
            stat.thisFP     = thisFP;
            stat.Port       = Port;

            for i = 1:num_entry
                if lb=="All"
                    ind_this = stat_session.thisFP==thisFP(i) & stat_session.Port==Port(i) & ismember(stat_session.Sessions, obj.Sessions(obj.Label~="Chemo"));
                else
                    ind_this = stat_session.thisFP==thisFP(i) & stat_session.Port==Port(i) & ismember(stat_session.Sessions, obj.Sessions(obj.Label==lb));
                end
                stat.N(i)          = mean(stat_session.N(ind_this), 'omitnan');
                stat.Median(i)     = mean(stat_session.Median(ind_this), 'omitnan');
                stat.Median_kde(i) = mean(stat_session.Median_kde(ind_this), 'omitnan');
                stat.IQR(i)        = mean(stat_session.IQR(ind_this), 'omitnan');
                stat.Q1(i)         = mean(stat_session.Q1(ind_this), 'omitnan');
                stat.Q3(i)         = mean(stat_session.Q3(ind_this), 'omitnan');
                stat.Median_sem(i) = std(stat_session.Median(ind_this), 'omitnan') / sqrt(sum(~isnan(stat_session.Median)));
            end
        end

        function obj = gatherAllStats(obj)
            Vars = ["RT", "MT", "HD", "CT"];
            Labels = ["All", "Control", "Chemo"];
            for v = 1:length(Vars)
                for l = 1:length(Labels)
                    obj.(Vars(v)+"Stat").(Labels(l)) = obj.gatherStat(Vars(v), Labels(l));
                end
            end
        end

        %% get PDFs and CDFs
        function data_pdf = getPDF(obj, data, bin_edges, var_name, lb, cal_ci)
            sz = size(data);
            data_pdf = cell(sz(1), sz(2));

            if cal_ci
                fprintf("... Calculate 95CI for PDF of %s at %s condition ... \n", var_name, lb);
            end
            for i = 1:sz(1)
                for j = 1:sz(2)
                    data_this = data{i, j};
                    data_pdf{i, j} = obj.calKDE(data_this, bin_edges, 'pdf', cal_ci);
                end
            end
        end % getPDF

        function data_cdf = getCDF(obj, data, bin_edges, var_name, lb, cal_ci)
            sz = size(data);
            data_cdf = cell(sz(1), sz(2));

            if cal_ci
                fprintf("... Calculate 95CI for CDF of %s at %s condition ... \n", var_name, lb);
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

            Vars = ["MT", "HD"];
            VarsB = ["MovementTime", "HoldDuration"];
            Labels = ["Control", "Chemo"];

            for l = 1:length(Labels)
                if ~contains(Labels(l), obj.Label)
                    continue;
                end
                for v = 1:length(Vars)
                    obj.(Vars(v)+"PDF").(Labels(l)) = obj.getPDF(obj.(Vars(v)+"Sorted").(Labels(l)), obj.Bins.(VarsB(v)), Vars(v), Labels(l), cal_ci);
                    obj.(Vars(v)+"CDF").(Labels(l)) = obj.getCDF(obj.(Vars(v)+"Sorted").(Labels(l)), obj.Bins.(VarsB(v)), Vars(v), Labels(l), cal_ci);
                end
                logST = cellfun(@log10, obj.STSplit.(Labels(l)), 'UniformOutput', false);
                obj.LogSTPDF.(Labels(l)) = obj.getPDF(logST, obj.Bins.ShuttleTimeLog, "LogST", Labels(l), cal_ci);
                obj.LogSTCDF.(Labels(l)) = obj.getCDF(logST, obj.Bins.ShuttleTimeLog, "LogST", Labels(l), cal_ci);
            end

        end % getAllKDEs
    end
end
