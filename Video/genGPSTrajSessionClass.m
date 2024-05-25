%%
clear; close all;

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

[TaskInd, tf] = listdlg("ListString", Tasks, "ListSize", [200, 200]);
if ~tf
    return
end

%%
View = "Top";

Drives = char('A':'Z');
Drives = string(Drives(:));
for i = 1:length(Drives)
    if isfolder(Drives(i)+":\OneDrive")
        ParentDir   = Drives(i) + ":\OneDrive\YuLab\Work\GPS\Data\";
        break;
    end
end

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
for i = 1:length(SessionFolders)
    tracking_file = dir(fullfile(SessionFolders(i), "*DLCTrackingOut_*.mat"));
    if isempty(tracking_file)
        SessionFolders(i) = "";
    end
end
SessionFolders(SessionFolders=="") = [];

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
for i = 1:length(Folders)

    tracking_file = dir(fullfile(Folders(i), "*DLCTrackingOut_*.mat"));
    behclass_file = dir(fullfile(Folders(i), "*BehSessionClass_*.mat"));

    if isempty(tracking_file)
        continue;
    end

    TrajSessionClass = GPSTrajSessionClass(fullfile(Folders(i), tracking_file.name), fullfile(Folders(i), behclass_file.name));

    fprintf("\n%s\n", TrajSessionClass.Session);

    TrajSessionClass.save();
end
