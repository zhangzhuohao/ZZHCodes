Drives = char('A':'Z');
Drives = string(Drives(:));
for i = 1:length(Drives)
    if isfolder(Drives(i)+":\OneDrive")
        ParentDir   = Drives(i) + ":\OneDrive\YuLab\Work\GPS\Data\";
        break;
    end
end

View = "Top";

VideoFolder = "Z:\YuLab\Work\GPS\Video\Kennard\GPS_08_KornblumHold1500SRTSelf";

ProtocolFolder = extractAfter(VideoFolder, 'Video\');
Anm = extractBefore(ProtocolFolder, "\");

DataFolder = fullfile(ParentDir, ProtocolFolder);

fprintf("\nCopy DLCTrackingOut files of %s\nfrom %s\nto   %s\n\n", Anm, VideoFolder, DataFolder);

VideoSessionFolders = get_folders(VideoFolder, "FolderType", "Session");
[~, Sessions] = arrayfun(@(x) fileparts(x), VideoSessionFolders);

for i = 1:length(Sessions)
    DLCTrackingFile = fullfile(VideoSessionFolders(i), View, "Clips", "DLCTrackingOutAuto.mat");
    if ~exist(DLCTrackingFile, 'file')
        continue;
    end

    load(DLCTrackingFile, "DLCTrackingOut");
    DataSessionFolder = fullfile(DataFolder, Sessions(i));

    copy_name = sprintf("DLCTrackingOut_%s_%s_%s.mat", View, Anm, Sessions(i));
    disp(copy_name);

    copy_path = fullfile(DataSessionFolder, copy_name);
    save(copy_path, "DLCTrackingOut");
end