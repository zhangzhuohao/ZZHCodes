clear; clc;

fprintf("\nVideo sorting ...\n");

%%
%               GPS_1           GPS_2
CamTopID    =   ["00G40298598", "L16636102"];
CamFrontID  =   ["00G40298581", "00D41933012"];
CamFieldID  =   ["G44627565"  , "DA0069619"];
CamInitID   =   ["00G40298619", "L16636084"];

%%
VideoFolder = uigetdir("F:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolder
    return;
end

VidSessionFolders = get_folders(VideoFolder, "FolderType", "Session");

%%
for b = 1:length(CamTopID)
    for s = 1:length(VidSessionFolders)
    
        VidTopFiles   = dir(VidSessionFolders(s) + "\*" + CamTopID(b)   + "*");
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
    
        VidFrontFiles = dir(VidSessionFolders(s) + "\*" + CamFrontID(b) + "*");
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

        VidFieldFiles = dir(VidSessionFolders(s) + "\*" + CamFieldID(b) + "*");
        if ~isempty(VidFieldFiles)
            FieldFolder = fullfile(VidSessionFolders(s), "Field");
            if ~isfolder(FieldFolder)
                mkdir(FieldFolder);
            end
            for f = 1:length(VidFieldFiles)
                vid_file = fullfile(VidFieldFiles(f).folder, VidFieldFiles(f).name);
                movefile(vid_file, FieldFolder);
            end
        end

        VidInitFiles = dir(VidSessionFolders(s) + "\*" + CamInitID(b) + "*");
        if ~isempty(VidInitFiles)
            InitFolder = fullfile(VidSessionFolders(s), "Init");
            if ~isfolder(InitFolder)
                mkdir(InitFolder);
            end
            for f = 1:length(VidInitFiles)
                vid_file = fullfile(VidInitFiles(f).folder, VidInitFiles(f).name);
                movefile(vid_file, InitFolder);
            end
        end
    
    end
end

%%
fprintf("\nDone\n")
