clear;
%%
View = "Front";
CheckEvent = "CentOut";
%%
TaskFolder = uigetdir('D:\YuLab\Work\GPS\Video\');
if ~TaskFolder
    return
end

%%
GatherFolder = fullfile(TaskFolder, "BadTrials");
if ~isfolder(GatherFolder)
    mkdir(GatherFolder);
end

GatherFolderLeft = fullfile(GatherFolder, "Left");
if ~isfolder(GatherFolderLeft)
    mkdir(GatherFolderLeft);
end
GatherFolderRight = fullfile(GatherFolder, "Right");
if ~isfolder(GatherFolderRight)
    mkdir(GatherFolderRight);
end

%%
SessionFolders = get_folders(TaskFolder, "FolderType", 'Session');

[~, SessionsAll] = arrayfun(@(x) fileparts(x), SessionFolders);
Sessions = unique(SessionsAll);

[SessionInd, tf] = listdlg("ListString", Sessions, "ListSize", [200, 200]);
if ~tf
    return
end

SessionFolders = SessionFolders(ismember(SessionsAll, Sessions(SessionInd)));

%%D:\YuLab\Work\GPS\Video\Morad\GPS_13_ThreeFPHoldSRT
for s = 1:length(SessionFolders)
%%
    ClipFolder = fullfile(SessionFolders(s), View, "Clips");
    if ~isfolder(ClipFolder)
        fprintf("no clip folder in %s", fullfile(SessionFolders(s), View));
    end

%%
    CheckFrameFolder = fullfile(ClipFolder, "CheckFrame");
    EventFolder = fullfile(CheckFrameFolder, CheckEvent);
    LeftFolder = fullfile(EventFolder, "Left");
    RightFolder = fullfile(EventFolder, "Right");

    if ~isfolder(EventFolder) || ~isfolder(CheckFrameFolder) || ~isfolder(LeftFolder) || ~isfolder(RightFolder)
        fprintf("check frames not extracted\n")
    end

%%
    bad_folder_left = fullfile(LeftFolder, "Bad");
    if isfolder(bad_folder_left)
        copyfile(bad_folder_left, GatherFolderLeft);
    end

    bad_folder_right = fullfile(RightFolder, "Bad");
    if isfolder(bad_folder_right)
        copyfile(bad_folder_right, GatherFolderRight);
    end
%%
end