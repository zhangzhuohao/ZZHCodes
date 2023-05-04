clear;

%% Choose the target tasks to generate ProgressClass
Tasks = [
    "Autoshaping";
    "Wait1Hold";
    "Wait1HoldSRT";
    "Wait1HoldCRT";
    "Wait2HoldSRT";
    "Wait2HoldCRT";
    "ThreeFPHoldCRT";
    "ThreeFPHoldSRT";
    "KornblumSRT";
    ];

[indx, tf] = listdlg("ListString", Tasks, "ListSize", [200, 200]);
if ~tf
    return
end

%%
% Use the uigetdir function to open a dialog box and allow the user to select the parent directory
ParentDir = uigetdir('E:\YuLab\Work\GPS\Data', 'Select parent directory');
if ~ParentDir
    return
end

%%
Entry = dir(ParentDir);

haveSubFolder = 0;
for e = 1:length(Entry)
    if isfolder(fullfile(Entry(e).folder, Entry(e).name)) && ~any(strcmp(Entry(e).name, {'.', '..'}))
        haveSubFolder = 1;
        break;
    end
end

if ~haveSubFolder
    Folders = string(ParentDir);
    infos   = split(Folders, ["_", filesep]);
    if ~any(strcmp(infos, Tasks(indx)))
        return
    end
else
    Folders = get_folders(ParentDir, "FolderType", 'Protocol', "Tasks", Tasks(indx));
end

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

        SessionClass = GPSSessionClass(file);
        SessionClass.save();
        SessionClass.updateANMInfo();
        SessionClass.print(SessionDataFigsFolder);

        BehavCsvName = ['GPSSessionTable' '_' SessionClass.Task '_' SessionClass.Subject, '.csv'];
        writetable(SessionClass.BehavTable, fullfile(file_path, BehavCsvName));

        InterCsvName = ['GPSSessionInter' '_' SessionClass.Task '_' SessionClass.Subject, '.csv'];
        writetable(SessionClass.Interruption, fullfile(file_path, InterCsvName));

    end
end

% close all;
