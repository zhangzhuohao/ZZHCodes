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
    "ThreeFPHoldWM";
    "KornblumHold1000SRT";
    "KornblumHold1000SRTSelf";
    "KornblumHold1500SRTSelf";
    "KornblumHold2000SRTSelf";
    "KornblumHold2000SRTEmpSelf";
    "KornblumHold2000SRTUnguideSelf";
    "KornblumHold2000SRTMixSelf"
    ];

[indx, tf] = listdlg("ListString", Tasks, "ListSize", [200, 200]);
if ~tf
    return
end

%% Use the uigetdir function to open a dialog box and allow the user to select the parent directory
Drives = char('A':'Z');
Drives = string(Drives(:));
for i = 1:length(Drives)
    if isfolder(Drives(i)+":\OneDrive")
        DataFolder  = Drives(i) + ":\OneDrive\YuLab\Work\GPS\Data\";
        break;
    end
end

ParentDir = uigetdir(DataFolder, "Select parent directory");
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
    Files = get_mat_files(ProtocolDir, "FileType", 'TrajSessionClass');

    if isempty(Files)
        continue
    end
    
    %
    TrajSessionClassAll = cell(length(Files), 1);
    for f = 1:length(Files)
        file = Files{f};
        load(file);
        TrajSessionClassAll{f} = obj;
    end
    clear obj;
    
    % Make sure that those session classes are sorted by their date
    SessionDate = cell2mat(cellfun(@(x) str2double(x.Session), TrajSessionClassAll, 'UniformOutput', false));
    [SessionDate, SortID] = sort(SessionDate);
    TrajSessionClassAll = TrajSessionClassAll(SortID);

    %
    TrajProgressClass = GPSTrajProgressClass(TrajSessionClassAll, ProtocolDir);
%     TrajProgressClass.save();
    
%     fig_progress = TrajProgressClass.plotProgress();
%     TrajProgressClass.print(fig_progress);
% %     ProgressClass.print("Show", ProtocolDir);
% 
%     BehavCsvName = fullfile(ProtocolDir, "GPSBehProgressClass_" + TrajProgressClass.Protocol + "_" + upper(TrajProgressClass.Subject) + ".csv");
%     writetable(TrajProgressClass.BehavTable, BehavCsvName);
% 
%     %
%     if all(ismember(["Control", "Chemo"], ProgressClass.Label))
%         ProgressClass.print("ChemoEffect", ProtocolDir);
%     end
end
