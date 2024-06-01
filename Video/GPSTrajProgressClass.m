classdef GPSTrajProgressClass < GPSTrajClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % Class object information
        ProtocolDir
        UpdateTime

        % Animal information
        Subject

        % Task information
        Task
        Protocol
        TargetFP

        % Session information
        Sessions
        NumSessions
        NumTrialsSession
        
        % Experiment information
        Treatment
        Dose
        Label
        
        % Trial behavior information
        TrialInfo

        % Trace time points
        TimeFromIn
        TimeFromOut
        TimeFromCue
        TimeWarpHD

        % Port location
        PortVec
        PortCent

        % Head postion
        AngleHead
        PosXHead
        PosYHead

        % Head speed
        AngSpeedHead
        SpeedXHead
        SpeedYHead
        SpeedHead
        SpeedDirHead

        % Head acceleration
        AngAccHead
        AccXHead
        AccYHead
        AccHead
        dPhiHead

        %
        TraceMatrix
        TraceInterp
        TraceMedian

        % Distance matrix
        DistMatLw
        DistMatDtw
    end

    properties (Dependent)
        SaveName
    end

    methods
        function obj = GPSTrajProgressClass(TrajSessionClassAll, ProtocolDir)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.ProtocolDir = ProtocolDir;

            % Animal information
            obj.Subject = unique(cellfun(@(x) string(x.Subject), TrajSessionClassAll));
            if length(obj.Subject)>1
                error("More than one animal");
            end

            % Task information
            obj.Task     = unique(cellfun(@(x) string(x.Task), TrajSessionClassAll));
            obj.Protocol = unique(cellfun(@(x) string(x.Protocol), TrajSessionClassAll));
            if length(obj.Protocol)>1
                error("More than one protocol");
            end
            obj.TargetFP = unique(cell2mat(cellfun(@(x) x.TargetFP, TrajSessionClassAll, 'UniformOutput', false)), "rows");
            if size(obj.TargetFP, 1)>1
                error("More than one set of target fp");
            end

            % Session information
            obj.Sessions            = cellfun(@(x) string(x.Session), TrajSessionClassAll);
            obj.NumSessions         = length(TrajSessionClassAll);
            obj.NumTrialsSession    = cellfun(@(x) x.NumTrials, TrajSessionClassAll);
            fprintf("\n%s: from %s to %s\n", obj.SaveName, obj.Sessions(1), obj.Sessions(end));

            % Experiment information
            obj.Treatment = cellfun(@(x) string(x.Treatment), TrajSessionClassAll);
            obj.Dose      = cellfun(@(x) x.Dose, TrajSessionClassAll);
            obj.Label     = cellfun(@(x) string(x.Label), TrajSessionClassAll);

            % Trial behavior information
            allTables = cellfun(@(x) {x.TrialInfo}, TrajSessionClassAll, 'UniformOutput', false);
            obj.TrialInfo = obj.splice_data(allTables);
            obj.TrialInfo = addvars(obj.TrialInfo{1}, (1:sum(obj.NumTrialsSession))', 'After', 'Trials', 'NewVariableNames', "Index");

            % Gather timepoints, trajectories and ports
            obj.gather_all_progress(TrajSessionClassAll);

            obj.gather_all_matrix();

        end

        %%
        function gather_all_progress(obj, TrajSessionClassAll)
            vars_to_gather = [
                "TimeFromIn" % time points
                "TimeFromOut"
                "TimeFromCue"
                "TimeWarpHD"
                "PortVec" % port location
                "PortCent"
                "AngleHead" % head position
                "PosXHead"
                "PosYHead"
                "AngSpeedHead" % head speed
                "SpeedXHead"
                "SpeedYHead"
                "SpeedHead"
                "SpeedDirHead"
                "AngAccHead" % head acceleration
                "AccXHead"
                "AccYHead"
                "AccHead"
                "dPhiHead"
                ];
            for i = 1:length(vars_to_gather)
                var_this = vars_to_gather(i);
                data_all = cellfun(@(x) {x.(var_this)}, TrajSessionClassAll, 'UniformOutput', false);
                data_spliced = obj.splice_data(data_all);
                data_spliced = data_spliced{1};
                obj.(var_this) = data_spliced;
            end
        end % gather_all_progress

        function gather_all_matrix(obj)
            
            info = obj.TrialInfo;

            Labels     = ["All", "Chemo", "Control", "PreControl", "LesionEarly", "LesionExtensive"];
            Aligns     = ["In", "Out", "Cue", "HD"];
            TimeTrace  = ["TimeFromIn", "TimeFromOut", "TimeFromCue", "TimeWarpHD"];
            TimeMatrix = ["TimeMatIn", "TimeMatOut", "TimeMatCue", "TimeMatHD"];

            for lb = 1:length(Labels)
                lb_this = Labels(lb);
                if lb_this=="All"
                    ind_this = ismember(info.Label, ["None", "Saline", "Control", "PreLesion", "PreControl"]);
                else
                    ind_this = info.Label==lb_this;
                end

                for a = 1:length(Aligns)
                    a_this = Aligns(a);
                    t_trace = obj.(TimeTrace(a));
                    t_matrix = obj.(TimeMatrix(a));
                    obj.TraceMatrix.(lb_this).(a_this) = obj.get_matrix(obj.Features, t_trace, t_matrix, ind_this);
                end
            end
        end

        %% 
        function save_name = get.SaveName(obj)
            save_name = sprintf("GPSTrajProgressClass_%s_%s", obj.Protocol, upper(obj.Subject));
        end

        %% Save
        function save(obj, copy_dir)
            save_path = fullfile(obj.SessionFolder, obj.SaveName);
            save(save_path, 'obj');

            if nargin==2
                copy_path = fullfile(copy_dir, obj.SaveName);
                copyfile(save_path+".mat", copy_path+".mat");
            end
        end % save

        %% Plots
        function print(obj, copy_dir)
            save_path = fullfile(obj.SessionFolder, obj.SaveName);
            fig = obj.plot();
            exportgraphics(fig, save_path+".png", 'Resolution', 600);
            exportgraphics(fig, save_path+".pdf", 'ContentType', 'vector');
            saveas(fig, save_path, 'fig');

            if nargin==2
                if ~isfolder(copy_dir)
                    mkdir(copy_dir);
                end
                copy_path = fullfile(copy_dir, obj.SaveName);
                copyfile(save_path+".pdf", copy_path+".pdf");
                copyfile(save_path+".png", copy_path+".png");
                copyfile(save_path+".fig", copy_path+".fig");
            end
        end % print

        function fig = plot(obj)
            try
                set_matlab_default
            catch
                disp('You do not have "set_matlab_default"' )
            end

            fig = feval(obj.Task+".plotTrajSession", obj);
        end % plot

    end
end