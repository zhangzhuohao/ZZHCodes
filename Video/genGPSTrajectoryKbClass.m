%%
clear; close all;

View = "Top";

AnmInfoFile = 'D:\YuLab\Work\GPS\Data\ANMInfo.xlsx';

ParentDir = uigetdir("F:\YuLab\Work\GPS\Video", "Select parent directory");
if ~ParentDir
    return
end

TrajSessionFolder = fullfile(ParentDir, "TrajSessionData");
if ~isfolder(TrajSessionFolder)
    mkdir(TrajSessionFolder);
end

%%
% indx = 5;
% ParentDir = "E:\YuLab\Work\GPS\Data\Gena";
SessionFolders = get_folders(ParentDir, "FolderType", 'Session');

[~, Sessions] = arrayfun(@(x) fileparts(x), SessionFolders);
[indx, tf] = listdlg("ListString", Sessions, "ListSize", [200, 200]);
if ~tf
    return
end

%%
ClipFolders = arrayfun(@(x) x+"\Top\Clips", SessionFolders(indx));

for i = 1:length(ClipFolders)

    tracking_file = dir(fullfile(ClipFolders(i), "*DLCTrackingOutAuto.mat"));

    if isempty(tracking_file)
        continue;
    end

    TrajectoryClass = GPSTrajectoryKbClass(fullfile(ClipFolders(i), tracking_file.name), AnmInfoFile);

    fprintf("\n%s\n", TrajectoryClass.Session);

    TrajectoryClass.print("HeatMap", TrajSessionFolder);
    TrajectoryClass.print("Trace", TrajSessionFolder);

%     TrajectoryClass.AngleHeadTraceInTest = TrajectoryClass.testTrace("In", 2000);
%     TrajectoryClass.AngleHeadTraceOutTest = TrajectoryClass.testTrace("Out", 2000);

    TrajectoryClass.save(TrajSessionFolder);

end
