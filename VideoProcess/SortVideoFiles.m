clear; clc;

fprintf("\nVideo sorting ...\n");

%%
CamTopID    =   "00G40298598";
CamFrontID  =   "00G40298581";

%%
VideoFolder = uigetdir("E:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolder
    return;
end

VidSessionFolders = get_folders(VideoFolder, "FolderType", "Session");

%%
for s = 1:length(VidSessionFolders)

    VidTopFiles   = dir(VidSessionFolders(s) + "\*" + CamTopID   + "*");
    if ~isempty(VidTopFiles)
        TopFolder = fullfile(VidSessionFolders(s), "Top");
        if ~isfolder(TopFolder)
            mkdir(TopFolder);
        end
        for f = 1:length(VidTopFiles)
            vid_file = fullfile(VidTopFiles(f).folder, VidTopFiles(f).name);
            movefile(vid_file, TopFolder);
        end
    end

    VidFrontFiles = dir(VidSessionFolders(s) + "\*" + CamFrontID + "*");
    if ~isempty(VidFrontFiles)
        FrontFolder = fullfile(VidSessionFolders(s), "Front");
        if ~isfolder(FrontFolder)
            mkdir(FrontFolder);
        end
        for f = 1:length(VidFrontFiles)
            vid_file = fullfile(VidFrontFiles(f).folder, VidFrontFiles(f).name);
            movefile(vid_file, FrontFolder);
        end
    end

end

