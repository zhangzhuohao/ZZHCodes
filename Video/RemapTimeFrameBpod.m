clear;

%%
fprintf("\nRemap tFrames2Bpodms ...\n");

VideoFolderParent = uigetdir("Z:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolderParent
    return;
end

Sessions = get_folders(VideoFolderParent, "FolderType", "Session");

%%
for s = 1:length(Sessions)

    clear("FrameTable", "FrameTableFront", "FrameTableField");

    VideoFolder = Sessions(s);
    disp(VideoFolder);

    SaveNameUsedData = fullfile(VideoFolder, "UsedData.mat");
    if ~exist(SaveNameUsedData, 'file')
        continue;
    end
    load(SaveNameUsedData);

    VideoFolderTop   = fullfile(VideoFolder, 'Top');
    VideoFolderFront = fullfile(VideoFolder, 'Front');
    VideoFolderField = fullfile(VideoFolder, 'Field');

    ClipInfo.LagEvent = "InitInTime";

    %%
    tMarkingEvent  = BehTable.(ClipInfo.MarkingEvent)  + BehTable.TrialStartTime; % in seconds
    tMarkingLag    = BehTable.(ClipInfo.MarkingEvent)  - BehTable.(ClipInfo.LagEvent); % in second
    tCheckingEvent = BehTable.(ClipInfo.CheckingEvent) + BehTable.TrialStartTime; % in seconds

    %%
    figure(99); clf(99); imshow(imread((fullfile(VideoFolder, "LED_On.png"))));

    [tLEDon, FigLEDon] = FindLEDonGPS(FrameTable.tsROI, FrameTable.SummedROI);
    SaveNameFigLEDon = fullfile(VideoFolder, "LED_On_remap");
    print(FigLEDon, '-dpng', SaveNameFigLEDon);

    [IndOut, FigAlign] = findseqmatchrev(tMarkingEvent*1000, tLEDon, 1);
    SaveNameFigAlign = fullfile(VideoFolder, "Alignment_remap");
    print(FigAlign, '-dpng', SaveNameFigAlign);

    tLEDon(isnan(IndOut)) = []; % remove them
    IndOut(isnan(IndOut)) = []; % at this point, each LEDon time can be mapped to a trigger time in b (Indout)

    tFrames2Bpodms = MapVidFrameTime2Bpod(tLEDon, tMarkingEvent(IndOut)*1000, tMarkingLag(IndOut)*1000, FrameTable.tsROI);

    tFrames2Bpodms_old = FrameTable.tFrames2Bpodms;
    FrameTable.tFrames2Bpodms = tFrames2Bpodms;

    SaveNameFrameTable = fullfile(VideoFolder, "FrameInfo.csv");
    writetable(FrameTable, SaveNameFrameTable);

    %%
    FrameInfo.tLEDon = tLEDon;
    FrameInfo.Indout = IndOut;
    FrameInfo.tLEDonBpod = tMarkingEvent(IndOut); % LED timing in Bpod's world
    FrameInfo.tFramesInBpod = tFrames2Bpodms;

    SaveNameFrameInfo = fullfile(VideoFolder, "FrameInfo.mat");
    save(SaveNameFrameInfo, "FrameInfo");

    %%
    view_front = 0;
    if exist("FrameTableFront", "var")
        view_front = 1;
        tsROIFront        = zeros(length(tFrames2Bpodms), 1);

        for i = 1:length(tFrames2Bpodms)
            [~, id] = min(abs(FrameTableFront.tsROI - FrameTable.tsROI(i)));
            tsROIFront(i)          =   FrameTableFront.tsROI(id);
        end

        tFrames2BpodmsFront = tFrames2Bpodms;

        [tsROIFront, id] = unique(tsROIFront);
        tFrames2BpodmsFront = tFrames2BpodmsFront(id);

        FrameTableFront.tFrames2Bpodms = tFrames2BpodmsFront;

        SaveNameFrameTableFront = fullfile(VideoFolder, "FrameInfoFront.csv");
        writetable(FrameTableFront, SaveNameFrameTableFront);
    end
    %%
    view_field = 0;
    if exist("FrameTableField", "var")
        view_field = 1;
        tsROIField        = zeros(length(tFrames2Bpodms), 1);
        for i = 1:length(tFrames2Bpodms)
            [~, id] = min(abs(FrameTableField.tsROI - FrameTable.tsROI(i)));
            tsROIField(i)          =   FrameTableField.tsROI(id);
        end

        tFrames2BpodmsField = tFrames2Bpodms;

        [tsROIField, id] = unique(tsROIField);
        tFrames2BpodmsField = tFrames2BpodmsField(id);

        FrameTableField.tFrames2Bpodms = tFrames2BpodmsField;

        SaveNameFrameTableField = fullfile(VideoFolder, "FrameInfoField.csv");
        writetable(FrameTableField, SaveNameFrameTableField);
    end

    %%
    switch view_front
        case 1
            switch view_field
                case 1
                    save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'FrameTableFront', 'FrameTableField', 'SessionInfo', 'ClipInfo', 'ClipInfoField', 'ScnScale');
                case 0
                    save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'FrameTableFront', 'SessionInfo', 'ClipInfo', 'ScnScale');
            end
        case 0
            save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'SessionInfo', 'ClipInfo', 'ScnScale');
    end

    %% remap vid-meta file - TopView
    VidMetaFilesTop = dir(fullfile(VideoFolderTop, "Clips", "*HoldTrial*.mat"));
    disp("Top");
    for i = 1:length(VidMetaFilesTop)
        VidMetaPath = fullfile(VidMetaFilesTop(i).folder, VidMetaFilesTop(i).name);
        load(VidMetaPath, "VidMeta");

        [~, id_start] = min(abs(FrameTable.tFrames2Bpodms - VidMeta.FrameTimesB(1)));
        id_clips = id_start:id_start+length(VidMeta.FrameTimesB)-1;

        fprintf("\nTime diff before: ");
        time_diff_before = unique(round(diff(VidMeta.FrameTimesB)));
        for j = 1:length(time_diff_before)
            fprintf("%d ", time_diff_before(j));
        end

        VidMeta.FrameTimesB = FrameTable.tFrames2Bpodms(id_clips);
        VidMeta.RemapOn     = date();
        if length(VidMeta.FrameTimesB)~=length(VidMeta.FrameIndx)
            error('not match');
        end

        fprintf(";\t after: ");
        time_diff_after = unique(round(diff(VidMeta.FrameTimesB)));
        for j = 1:length(time_diff_after)
            fprintf("%d ", time_diff_after(j));
        end

        save(VidMetaPath, "VidMeta");
    end

    fprintf("\n\n------------------------------\n");
    fprintf("\n------------------------------\n\n");

end
%%
