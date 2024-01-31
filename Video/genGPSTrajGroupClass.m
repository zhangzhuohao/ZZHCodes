%%
clear; close all;

View = "Top";

ParentDir = uigetdir("F:\YuLab\Work\GPS\Video", "Select parent directory");
if ~ParentDir
    return
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

ClipFolders = arrayfun(@(x) x+"\Top\Clips", SessionFolders(indx));

%%
TrajClassFiles = arrayfun(@(x) dir(fullfile(x, "GPSTrajectoryClass*.mat")), ClipFolders);
AllTrajClass = arrayfun(@(x, y) load(fullfile(x, y.name), "obj"), ClipFolders, TrajClassFiles);

TrajGroupClass = GPSTrajGroupClass(AllTrajClass);
TrajGroupClass.save(ParentDir);

TrajGroupClass.print("HeatMap", ParentDir);
% TrajGroupClass.print("Trace", ParentDir);

