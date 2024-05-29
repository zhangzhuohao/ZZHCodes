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
        continue;
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
fprintf("\n%d bad label frames have been gathered to %s.\n", n_bad_labels, GatherFolder);

% %% Create videos from video
% 
% VidClipName = fullfile(GatherFolder, "RefineVideo.avi");
% 
% writerObj = VideoWriter(VidClipName);
% writerObj.FrameRate = 25; % this is 2 x slower
% 
% % set the seconds per image
% % open the video writer
% 
% open(writerObj);
% for ii=1:n_bad_labels
%     currentfilename = fullfile(GatherFolder, bad_label_frames(ii).name);
%     currentimage = imread(currentfilename);
%     writeVideo(writerObj, currentimage);
% end
% close(writerObj);
% 
% fprintf("\nRefine video created.\n");

