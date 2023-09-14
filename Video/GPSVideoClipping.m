%%  File information
% Use this func to find video files and associated timestamp files.
% go to folder that includes video files, use 'ShowAviFiles' to get vid
% file names and folder name

% revised by ZZH, 5/5/2023

% run ShowAviFiles.m, select any vid file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Fill in the following information %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

ScnScale = 1.25;

%%
fprintf("\nGPS video clipping ...\n");

VideoFolder = uigetdir("D:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolder
    return;
end

SaveNameUsedData = fullfile(VideoFolder, "UsedData.mat");
if exist(SaveNameUsedData, "file")

    load(SaveNameUsedData);

    fprintf("\nUsed data already exist");
    fprintf("\n----------------------------------------");
    fprintf("\n----------------------------------------");

    ExportVideoClipFromAvi(BehTable, IntTable, FrameTable, SessionInfo, ClipInfo, ScnScale, 0, "Top", 0);

    ExportVideoClipFromAvi(BehTable, IntTable, FrameTableFront, SessionInfo, ClipInfo, ScnScale, 0, "Front", 1);

    fprintf("\n----------------------------------------");
    fprintf("\n----------------------------------------");

    fprintf("\nDone.\n\n");
    return
end

%%
ViewFolders     =   dir(VideoFolder);

view_front      =   0;
view_top        =   0;

for i = 1:length(ViewFolders)
    entry = fullfile(ViewFolders(i).folder, ViewFolders(i).name);
    if isfolder(entry)
        switch upper(ViewFolders(i).name)
            case {'FRONT'}
                view_front          =   1;
                VideoFolderFront    =   entry;
            case {'TOP'}
                view_top            =   1;
                VideoFolderTop      =   entry;
            case {'FIELD'}
                view_field          =   1;
                VideoFolderField    =   entry;
        end
    end
end
clear i entry

if ~view_top
    fprintf("\nFind no top-view vedios.\n");
    return
else
    fprintf("\nTop-view vedios folder: \n%s\n", VideoFolderTop);
    [vidFilesTop, tsFilesTop, ~] = ShowAviFiles(VideoFolderTop, 'Top');
    fprintf(" Videos\t\tTimestamps\n");
    cellfun(@(vid, ts) fprintf(" %s\t%s\n", vid, ts), vidFilesTop, tsFilesTop);
end

if ~view_front
    fprintf("\nFind no front-view vedios.\n");
else
    fprintf("\nFront-view vedios folder: \n%s\n", VideoFolderFront);
    [vidFilesFront, tsFilesFront, ~] = ShowAviFiles(VideoFolderFront, 'Front');
    fprintf(" Videos\t\tTimestamps\n");
    cellfun(@(vid, ts) fprintf(" %s\t%s\n", vid, ts), vidFilesFront, tsFilesFront);
end

if ~view_field
    fprintf("\nFind no field-view vedios.\n");
else
    fprintf("\nField-view vedios folder: \n%s\n", VideoFolderField);
    [vidFilesField, tsFilesField, ~] = ShowAviFiles(VideoFolderField, 'Field');
    fprintf(" Videos\n");
    cellfun(@(vid) fprintf(" %s\n", vid), vidFilesField);
    fprintf("Timestamps\n");
    cellfun(@(ts) fprintf(" %s\n", ts), tsFilesField);
end

%%
vidInfo         =   split(VideoFolder, filesep);
ANM             =   string(vidInfo{end-1});
Session         =   string(vidInfo{end});

ANMInfoFile     =   "D:\YuLab\Work\GPS\Data\ANMInfo.xlsx";
ANMInfo         =   readtable(ANMInfoFile, "Sheet", ANM, "TextType", "string");
ANMInfo.Session =   string(ANMInfo.Session);

if isempty(ANMInfo.Session==Session)
    fprintf("\nCannot find session information.\n");
    return
elseif ANMInfo(ANMInfo.Session==Session, :).Task == "None"
    fprintf("\nThe session class has not been generated.\n");
    return
else
    SessionInfo = ANMInfo(ANMInfo.Session==Session, :);
    SessionInfo = addvars(SessionInfo, ANM, 'Before', "Session");
    fprintf("\n");
    disp(SessionInfo(:, 1:6));
end

%%
MarkingEvent    =   'CentInTime';       % use this time to align bpod and video ts
CheckingEvent   =   'ChoiceCueTime_1';  % use this time to align bpod and video ts
Pre             =   1000; % ms pre event time for video clips
Post            =   3000; % ms post event time for video clips
switch SessionInfo.Task
    case {'Autoshaping'}
        VideoEvent  =   'CentOutTime';      % make video clips based on this event
    otherwise
        VideoEvent  =   'CentInTime';       % make video clips based on this event
end

ClipInfo                    =   struct();
ClipInfo.VideoFolder        =   VideoFolder;
ClipInfo.VideoFolderTop     =   VideoFolderTop;
ClipInfo.VideoFolderFront   =   VideoFolderFront;
ClipInfo.MarkingEvent       =   MarkingEvent;
ClipInfo.CheckingEvent      =   CheckingEvent;
ClipInfo.VideoEvent         =   VideoEvent;
ClipInfo.Pre                =   Pre;
ClipInfo.Post               =   Post;

ClipInfoField                   =   ClipInfo;
ClipInfoField.Pre               =   10000;
ClipInfoField.VideoFolderField  =   VideoFolderField;
ClipInfoField.Post              =   0;

%%
csvFileBehavior     =   dir(SessionInfo.SessionFolder+"\*SessionTable*.csv");
csvFileInterrupt    =   dir(SessionInfo.SessionFolder+"\*InterruptTable*.csv");

if isempty(csvFileBehavior) || isempty(csvFileInterrupt)
    fprintf("Create the SessionClass and generate Behav and Interrput table in advance.");
    return
end

BehTableName        =   csvFileBehavior(1).name;
BehTableFile        =   fullfile(SessionInfo.SessionFolder, BehTableName);
BehTable            =   readtable(BehTableFile);

IntTableName        =   csvFileInterrupt(1).name;
IntTableFile        =   fullfile(SessionInfo.SessionFolder, IntTableName);
IntTable            =   readtable(IntTableFile);

tMarkingEvent       =   BehTable.(MarkingEvent)  + BehTable.TrialStartTime; % in seconds
tCheckingEvent      =   BehTable.(CheckingEvent) + BehTable.TrialStartTime; % in seconds

%%  ROI extraction
% Define mask
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
[Mask, FigMask] = ExtractMask(fullfile(VideoFolderTop, vidFilesTop{1}), [1000 2000]);
% save this fig
SaveNameFigMask = fullfile(VideoFolder, "ROI_Mask");
print(FigMask, '-dpng', SaveNameFigMask);

%% Based on "mask", in top view, extract pixel intensity in ROI from all frames
tsROI           =   [];
SummedROI       =   [];
AviFrameIndx    =   [];
AviFileIndx     =   [];

tic

fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
fprintf("\nExtracting marking events ...\n");

for i = 1:length(vidFilesTop)

    fileID          =   fopen(fullfile(VideoFolderTop, tsFilesTop{i}), 'r');
    formatSpec      =   '%f' ;
    NumOuts         =   fscanf(fileID, formatSpec); % this contains frame time (in ms) and frame index
    fclose(fileID);

    ind_brk         =   find(NumOuts==0);
    FrameTs         =   NumOuts(1:ind_brk-1);   % frame timestamps
    FrameIdx        =   NumOuts(ind_brk+1:end); % frame index

    filename = fullfile(VideoFolderTop, vidFilesTop{i});
    vidObj   = VideoReader(filename);

    for k = 1:vidObj.NumFrames
        thisFrame = read(vidObj, k);
        thisFrame = thisFrame(:, :, 1);
        roi_k     = sum(thisFrame(Mask));

        tsROI           =   [tsROI; FrameTs(k)];
        SummedROI       =   [SummedROI; roi_k];
        AviFileIndx     =   [AviFileIndx; i];
        AviFrameIndx    =   [AviFrameIndx; k];
    end
    toc
end
clear i k roi_k
tsStart         = tsROI(1);
tsROI           = tsROI - tsStart; % onset normalized to 0

FrameInfo                   =   struct();
FrameInfo.tframe            =   tsROI;
FrameInfo.mask              =   Mask;
FrameInfo.ROI               =   SummedROI;
FrameInfo.AviFile           =   vidFilesTop;
FrameInfo.AviFileIndx       =   AviFileIndx;
FrameInfo.AviFrameIndx      =   AviFrameIndx;

% Save for now because it takes a long time to get tsROI
SaveNameFrameInfo = fullfile(VideoFolder, "FrameInfo.mat");
save(SaveNameFrameInfo, "FrameInfo");

MyVidFiles      = cell(length(AviFileIndx), 1);
for i = 1:length(MyVidFiles)
    MyVidFiles{i} = vidFilesTop{AviFileIndx(i)};
end
clear i

%%
if view_front
    tsROIFront_org          =   [];
    AviFileIndxFront_org    =   [];
    AviFrameIndxFront_org   =   [];
    for i = 1:length(vidFilesFront)
        fileID          =   fopen(fullfile(VideoFolderFront, tsFilesFront{i}), 'r');
        formatSpec      =   '%f' ;
        NumOuts         =   fscanf(fileID, formatSpec); % this contains frame time (in ms) and frame index
        fclose(fileID);

        ind_brk         =   find(NumOuts==0);
        FrameTs         =   NumOuts(1:ind_brk-1);   % frame timestamps
        FrameIdx        =   NumOuts(ind_brk+1:end); % frame index

        filename = fullfile(VideoFolderFront, vidFilesFront{i});
        vidObj   = VideoReader(filename);

        for k = 1:vidObj.NumFrames
            tsROIFront_org          =   [tsROIFront_org; FrameTs(k)];
            AviFileIndxFront_org    =   [AviFileIndxFront_org; i];
            AviFrameIndxFront_org   =   [AviFrameIndxFront_org; k];
        end
    end
    clear i
    tsROIFront_org  = tsROIFront_org - tsStart;

    tsROIFront          =   zeros(length(AviFileIndx), 1);
    AviFileIndxFront    =   zeros(length(AviFileIndx), 1);
    AviFrameIndxFront   =   zeros(length(AviFileIndx), 1);
    for i = 1:length(AviFileIndx)
        [~, id] = min(abs(tsROIFront_org - FrameInfo.tframe(i)));
        tsROIFront(i)          =   tsROIFront_org(id);
        AviFileIndxFront(i)    =   AviFileIndxFront_org(id);
        AviFrameIndxFront(i)   =   AviFrameIndxFront_org(id);
    end
    clear i id

    FrameInfoFront              =   struct();
    FrameInfoFront.tframe       =   tsROIFront;
    FrameInfoFront.AviFileIndx  =   AviFileIndxFront;
    FrameInfoFront.AviFrameIndx =   AviFrameIndxFront;

    SaveNameFrameInfoFront      =   fullfile(VideoFolder, "FrameInfoFront.mat");
    save(SaveNameFrameInfoFront, "FrameInfoFront");

    MyVidFilesFront = cell(length(AviFileIndx), 1);
    for i = 1:length(MyVidFiles)
        MyVidFilesFront{i} = vidFilesFront{AviFileIndxFront(i)};
    end
    clear i
end

%%
if view_field
    tsROIField_org          =   [];
    AviFileIndxField_org    =   [];
    AviFrameIndxField_org   =   [];

    for i = 1:length(tsFilesField)

        fileID          =   fopen(fullfile(VideoFolderField, tsFilesField{i}), 'r');
        formatSpec      =   '%f' ;
        NumOuts         =   fscanf(fileID, formatSpec); % this contains frame time (in ms) and frame index
        fclose(fileID);

        ind_brk         =   find(NumOuts==0);
        FrameTs         =   NumOuts(1:ind_brk-1);   % frame timestamps
        FrameIdx        =   NumOuts(ind_brk+1:end); % frame index
    end

    j = 0;
    for i = 1:length(vidFilesField)
        
        filename = fullfile(VideoFolderField, vidFilesField{i});
        vidObj   = VideoReader(filename);
        
        for k = 1:vidObj.NumFrames
            j = j + 1;
            tsROIField_org          =   [tsROIField_org; FrameTs(j)];
            AviFileIndxField_org    =   [AviFileIndxField_org; i];
            AviFrameIndxField_org   =   [AviFrameIndxField_org; k];
        end
    end
    clear i
    tsROIField_org  = tsROIField_org - tsStart;

    tsROIField          =   zeros(length(AviFileIndx), 1);
    AviFileIndxField    =   zeros(length(AviFileIndx), 1);
    AviFrameIndxField   =   zeros(length(AviFileIndx), 1);
    for i = 1:length(AviFileIndx)
        [~, id] = min(abs(tsROIField_org - tsROI(i)));
        tsROIField(i)          =   tsROIField_org(id);
        AviFileIndxField(i)    =   AviFileIndxField_org(id);
        AviFrameIndxField(i)   =   AviFrameIndxField_org(id);
    end
    clear i id

    FrameInfoField              =   struct();
    FrameInfoField.tframe       =   tsROIField;
    FrameInfoField.AviFileIndx  =   AviFileIndxField;
    FrameInfoField.AviFrameIndx =   AviFrameIndxField;

    SaveNameFrameInfoField      =   fullfile(VideoFolder, "FrameInfoField.mat");
    save(SaveNameFrameInfoField, "FrameInfoField");

    MyVidFilesField = cell(length(AviFileIndx), 1);
    for i = 1:length(MyVidFiles)
        MyVidFilesField{i} = vidFilesField{AviFileIndxField(i)};
    end
    clear i
end

%% Alignment: tsROI and MarkingEvent 
% use threshold method to find LED onset
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
[tLEDon, FigLEDon] = FindLEDonGPS(tsROI, SummedROI);
SaveNameFigLEDon = fullfile(VideoFolder, "LED_On");
print(FigLEDon, '-dpng', SaveNameFigLEDon);

%% 
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------\n");
[IndOut, FigAlign] = findseqmatchrev(tMarkingEvent*1000, tLEDon);
SaveNameFigAlign = fullfile(VideoFolder, "Alignment");
print(FigAlign, '-dpng', SaveNameFigAlign);

% these LEDon times are the ones that cannot be matched to trigger. It must be a false positive signal that was picked up by mistake in "tLEDon = FindLEDon(tsROI, SummedROI);"
tLEDon(isnan(IndOut)) = []; % remove them
IndOut(isnan(IndOut)) = []; % at this point, each LEDon time can be mapped to a trigger time in b (Indout)

FrameInfo.tLEDon = tLEDon;
FrameInfo.Indout = IndOut;
FrameInfo.tLEDonBpod = tMarkingEvent(IndOut); % LED timing in Bpod's world

%% Now, let's redefine the frame time. Each frame time should be re-mapped to the timespace in bpod.
% all frame times are here: FrameInfo.tframe
tFrames2Bpodms = MapVidFrameTime2Bpod(tLEDon, tMarkingEvent(IndOut)*1000, tsROI);
FrameTable = table(tsROI, tFrames2Bpodms, SummedROI, AviFileIndx, AviFrameIndx, MyVidFiles);

SaveNameFrameTable = fullfile(VideoFolder, "FrameInfo.csv");
writetable(FrameTable, SaveNameFrameTable);

FrameInfo.tFramesInBpod = tFrames2Bpodms;

%%
if view_front
    FrameTableFront = table(tsROIFront, tFrames2Bpodms, AviFileIndxFront, AviFrameIndxFront, MyVidFilesFront, ...
        'VariableNames', {'tsROI', 'tFrames2Bpodms', 'AviFileIndx', 'AviFrameIndx', 'MyVidFiles'});
    [~, ia] = unique(FrameTableFront.tsROI);
    FrameTableFront = FrameTableFront(ia, :);

    SaveNameFrameTableFront = fullfile(VideoFolder, "FrameInfoFront.csv");
    writetable(FrameTableFront, SaveNameFrameTableFront);
end

%%
if view_field
    FrameTableField = table(tsROIField, tFrames2Bpodms, AviFileIndxField, AviFrameIndxField, MyVidFilesField, ...
        'VariableNames', {'tsROI', 'tFrames2Bpodms', 'AviFileIndx', 'AviFrameIndx', 'MyVidFiles'});
    [~, ia] = unique(FrameTableField.tsROI);
    FrameTableField = FrameTableField(ia, :);

    SaveNameFrameTableField = fullfile(VideoFolder, "FrameInfoField.csv");
    writetable(FrameTableField, SaveNameFrameTableField);
end

%% Check if another event is also aligned
% Previously, we have defined checking event:
% CheckingEvent = 'CenterLightTime'; % use this time to align bpod and video ts
% tCheckingEvent = eval(['thisTable.' CheckingEvent]); % in seconds
CheckAnotherEvent = 1;
if CheckAnotherEvent % still working on this
end

%% Make video clips
% Extrct           =       "ExportVideoFiles(b, FrameInfo, 'Event', 'Press', 'TimeRange', [2000 3000], 'SessionName',strrep(fileinfo.MED(1:10), '-', ''), 'RatName', MEDfile(27:strfind(MEDfile, '.')-1), 'Remake', 0)";
% MakeSheet   =      "MakeSpreadSheet(b, FrameInfo, 'Event', 'Press', 'TimeRange', [2000 3000], 'SessionName', strrep(fileinfo.MED(1:10), '-', ''), 'RatName', MEDfile(27:strfind(MEDfile, '.')-1))";

% make video clip based on two tables: a behavioral data table and a frame
% info table
% 

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
%%
clearvars -except SaveNameUsedData;
load(SaveNameUsedData);

%%
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");

ExportVideoClipFromAvi(BehTable, IntTable, FrameTable, SessionInfo, ClipInfo, ScnScale, 0, "Top", 0);

ExportVideoClipFromAvi(BehTable, IntTable, FrameTableFront, SessionInfo, ClipInfo, ScnScale, 0, "Front", 1);

ExportVideoClipFromAvi_Field(BehTable, IntTable, FrameTableField, SessionInfo, ClipInfoField, ScnScale, 0, "Field", 0);

fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
fprintf("\nDone.\n\n");

