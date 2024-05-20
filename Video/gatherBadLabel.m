clear;
%%
View = "Top";

%%
TaskFolder = uigetdir('D:\YuLab\Work\GPS\Video\');
if ~TaskFolder
    return
end

%%
GatherFolder = fullfile(TaskFolder, "BadLabels");
if ~isfolder(GatherFolder)
    mkdir(GatherFolder);
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
%%
for s = 1:length(SessionFolders)
%
    ClipFolder = fullfile(SessionFolders(s), View, "Clips");
    if ~isfolder(ClipFolder)
        fprintf("no clip folder in %s", fullfile(SessionFolders(s), View));
    end

%
    bad_folder = fullfile(ClipFolder, "BadLabels");

    if ~isfolder(bad_folder)
%         fprintf("bad labels not extracted\n");
        continue
    end

%
    if isfolder(bad_folder)
        copyfile(bad_folder, GatherFolder);
    end

%
end

%%
bad_label_frames = dir(GatherFolder+"\*.jpg");
n_bad_labels = length(bad_label_frames);
fprintf("\n%d bad label frames have been gathered.\n", n_bad_labels);

