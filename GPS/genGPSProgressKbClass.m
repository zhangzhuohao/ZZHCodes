clear;

%% Choose the target tasks to generate ProgressClass 
Tasks = [
    "Kornblum1000SRT";
    "Kornblum1000SRTSelf";
    "Kornblum1500SRTSelf";
    "Kornblum2000SRTSelf"
    ];

[indx, tf] = listdlg("ListString", Tasks, "ListSize", [200, 200]);
if ~tf
    return
end

%% Use the uigetdir function to open a dialog box and allow the user to select the parent directory
ParentDir = uigetdir("E:\YuLab\Work\GPS\Data", "Select parent directory");
if ~ParentDir
    return
end

% indx = 5;
% ParentDir = "E:\YuLab\Work\GPS\Data\Gena";
Folders = get_folders(ParentDir, "FolderType", 'Protocol', "Tasks", Tasks(indx));

%%
for d = 1:length(Folders)
    % Call the recursive function to get a list of all the .mat files in the parent directory and its subdirectories
    ProtocolDir = Folders(d);
    Files = get_mat_files(ProtocolDir, "FileType", 'SessionClass');

    if isempty(Files)
        continue
    end
    
    %%
    SessionClassAll = cell(length(Files), 1);
    for f = 1:length(Files)
        file = Files{f};
        load(file);
        SessionClassAll{f} = obj;
    end
    clear obj;
    
    % Make sure that those session classes are sorted by their date
    SessionDate = cell2mat(cellfun(@(x) str2double(x.Session), SessionClassAll, 'UniformOutput', false));
    [SessionDate, SortID] = sort(SessionDate);
    SessionClassAll = SessionClassAll(SortID);

    %%
    ProgressClass = GPSProgressKbClass(SessionClassAll);
    ProgressClass.save(ProtocolDir);
    
    ProgressClass.print("Progress", ProtocolDir);
%     ProgressClass.print("Show", ProtocolDir);
    
    if all(ismember(["Control", "Chemo"], ProgressClass.Label))
        ProgressClass.print("ChemoEffect", ProtocolDir);
    end

    BehavCsvName = fullfile(ProtocolDir, "GPSProgressClass_" + ProgressClass.Task + "_" + upper(ProgressClass.Subject) + ".csv");
    writetable(ProgressClass.BehavTable, BehavCsvName);
end
