clear;

%% Choose the target tasks to generate ProgressClass
Tasks = [
    "Kornblum1000SRT";
    "Kornblum1000SRTSelf";
    "Kornblum1500SRTSelf";
    "Kornblum2000SRTSelf";
    "Kornblum2000SRTEmpSelf";
    "Kornblum2000SRTUnguideSelf";
    "Kornblum2000SRTMixSelf"
    ];

[TaskInd, tf] = listdlg("ListString", Tasks, "ListSize", [200, 200]);
if ~tf
    return
end

%%
AnmInfoFile = "F:\YuLab\Work\GPS\Data\ANMInfo.xlsx";

% Use the uigetdir function to open a dialog box and allow the user to select the parent directory
ParentDir = 'F:\YuLab\Work\GPS\Data';

Entry = dir(ParentDir);
ANMFolders = [];
for e = 1:length(Entry)
    if Entry(e).isdir && ~any(strcmp(Entry(e).name, {'.', '..'}))
        ANMFolders = [ANMFolders; string(fullfile(Entry(e).folder, Entry(e).name))];
    end
end
[~, ANMs] = arrayfun(@(x) fileparts(x), ANMFolders);

[ANMInd, tf] = listdlg("ListString", ANMs, "ListSize", [200, 200]);
if ~tf
    return
end

ANMFolders = ANMFolders(ANMInd);

%%
ProtocolFolders = [];
for f = 1:length(ANMFolders)
    ProtocolFolders = [ProtocolFolders; get_folders(ANMFolders(f), "FolderType", 'Protocol', "Tasks", Tasks(TaskInd))];
end

%%
SessionFolders = [];
for f = 1:length(ProtocolFolders)
    SessionFolders = [SessionFolders; get_folders(ProtocolFolders(f), "FolderType", 'Session')];
end

%%
[~, SessionsAll] = arrayfun(@(x) fileparts(x), SessionFolders);
Sessions = unique(SessionsAll);

[SessionInd, tf] = listdlg("ListString", Sessions, "ListSize", [200, 200]);
if ~tf
    return
end

%%
Folders = SessionFolders(ismember(SessionsAll, Sessions(SessionInd)));

%%
for d = 1:length(Folders)
    % Call the recursive function to get a list of all the .mat files in the parent directory and its subdirectories
    ProtocolDir = Folders(d);
    Files = get_mat_files(ProtocolDir, "FileType", 'Bpod');

    for f = 1:length(Files)
        file = Files{f};
        file_info = split(file, filesep);
        [file_path, file_name, ext] = fileparts(file);
        [protocol_folder, file_date] = fileparts(file_path);

        SessionDataFigsFolder = fullfile(protocol_folder, 'SessionDataFigs');
        if ~isfolder(SessionDataFigsFolder)
            mkdir(SessionDataFigsFolder);
        end

        SessionClass = GPSSessionKbClass(file, AnmInfoFile, 0);
        SessionClass.save();
        SessionClass.updateANMInfo();
        SessionClass.print(SessionDataFigsFolder);

        BehavCsvName = ['GPSSessionTable_' SessionClass.TaskName '_' SessionClass.Subject, '.csv'];
        writetable(SessionClass.BehavTable, fullfile(file_path, BehavCsvName));

        InterCsvName = ['GPSInterruptTable_' SessionClass.TaskName '_' SessionClass.Subject, '.csv'];
        writetable(SessionClass.Interruption, fullfile(file_path, InterCsvName));

    end
end

% close all;
