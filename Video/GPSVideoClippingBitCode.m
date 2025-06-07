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

ScnScl = 1.5;

%%
fprintf("\nGPS video clipping ...\n");

VideoFolder = uigetdir("E:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolder
    return;
end

SaveNameUsedData = fullfile(VideoFolder, "UsedData.mat");
if exist(SaveNameUsedData, "file")

    load(SaveNameUsedData);

    fprintf("\nUsed data already exist");
    fprintf("\n----------------------------------------");
    fprintf("\n----------------------------------------");

    ExportVideoClipFromAvi(BehTable, IntTable, FrameTable, SessionInfo, ClipInfo, ScnScl, 1, "Top", 0);

    %     ExportVideoClipFromAvi(BehTable, IntTable, FrameTableFront, SessionInfo, ClipInfo, ScnScale, 0, "Front", 1);

    ExportVideoClipFromAvi_Field(BehTable, FrameTableField, SessionInfo, ClipInfoField, ScnScl, 1, 0);
    ExportVideoClipFromAvi_Init(BehTable, FrameTableInit, SessionInfo, ClipInfoInit, ScnScl, 1, 0);

    fprintf("\n----------------------------------------");
    fprintf("\n----------------------------------------");

    fprintf("\nDone.\n\n");
    return
end

%%
ViewFolders     =   dir(VideoFolder);

view_front      =   0;
view_top        =   0;
view_field      =   0;
view_init       =   0;

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
            case {'INIT'}
                view_init           =   1;
                VideoFolderInit     =   entry;
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
    VideoFolderFront = [];
else
    fprintf("\nFront-view vedios folder: \n%s\n", VideoFolderFront);
    [vidFilesFront, tsFilesFront, ~] = ShowAviFiles(VideoFolderFront, 'Front');
    fprintf(" Videos\t\tTimestamps\n");
    cellfun(@(vid, ts) fprintf(" %s\t%s\n", vid, ts), vidFilesFront, tsFilesFront);
end

if ~view_field
    fprintf("\nFind no field-view vedios.\n");
    VideoFolderField = [];
else
    fprintf("\nField-view vedios folder: \n%s\n", VideoFolderField);
    [vidFilesField, tsFilesField, ~] = ShowAviFiles(VideoFolderField, 'Field');
    fprintf(" Videos\n");
    cellfun(@(vid) fprintf(" %s\n", vid), vidFilesField);
    fprintf("Timestamps\n");
    cellfun(@(ts) fprintf(" %s\n", ts), tsFilesField);
end

if ~view_init
    fprintf("\nFind no init-view vedios.\n");
    VideoFolderInit = [];
else
    fprintf("\nInit-view vedios folder: \n%s\n", VideoFolderInit);
    [vidFilesInit, tsFilesInit, ~] = ShowAviFiles(VideoFolderInit, 'Init');
    fprintf(" Videos\n");
    cellfun(@(vid) fprintf(" %s\n", vid), vidFilesInit);
    fprintf("Timestamps\n");
    cellfun(@(ts) fprintf(" %s\n", ts), tsFilesInit);
end

%%
Drives = string(char('A':'Z')');
for i = 1:length(Drives)
    if isfolder(Drives(i)+":\OneDrive")
        DataFolder  = Drives(i) + ":\OneDrive\YuLab\Work\GPS\Data\";
        break;
    end
end

vidInfo         = split(VideoFolder, filesep);
ANM             = string(vidInfo{end-2});
Session         = string(vidInfo{end});

ANMInfoFile     = fullfile(DataFolder, "ANMInfo.xlsx");
ANMInfo         = readtable(ANMInfoFile, "Sheet", ANM, "TextType", "string");
ANMInfo.Session = string(ANMInfo.Session);

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

SessionFolderPart = extractAfter(SessionInfo.SessionFolder, "Data\");
SessionFolder = fullfile(DataFolder, SessionFolderPart);

%%
MarkingEvent    =   'CentInTime';   % use this time to match led-on time and align bpod and video ts
LagEvent        =   'InitInTime';   % use this time to align bpod and video ts
CheckingEvent   =   'CentOutTime';  % use this time to check ts mapping
Pre             =   1000; % ms pre event time for video clips
Post            =   3000; % ms post event time for video clips
% switch SessionInfo.Task
%     case {'Autoshaping'}
%         VideoEvent  =   'CentOutTime';      % make video clips based on this event
%     otherwise
VideoEvent  =   'CentInTime';       % make video clips based on this event
% end

ClipInfo                    =   struct();
ClipInfo.VideoFolder        =   VideoFolder;
ClipInfo.VideoFolderTop     =   VideoFolderTop;
ClipInfo.VideoFolderFront   =   VideoFolderFront;
ClipInfo.MarkingEvent       =   MarkingEvent;
ClipInfo.LagEvent           =   LagEvent;
ClipInfo.CheckingEvent      =   CheckingEvent;
ClipInfo.VideoEvent         =   VideoEvent;
ClipInfo.Pre                =   Pre;
ClipInfo.Post               =   Post;

ClipInfoField                   =   ClipInfo;
ClipInfoField.Pre               =   10000;
ClipInfoField.VideoFolderField  =   VideoFolderField;
ClipInfoField.Post              =   500;

ClipInfoInit                   =   ClipInfo;
ClipInfoInit.Pre               =   2000;
ClipInfoInit.VideoFolderInit   =   VideoFolderInit;
ClipInfoInit.VideoEvent        =   'InitOutTime';
ClipInfoInit.Post              =   2000;

%%
csvFileBehavior     =   dir(SessionFolder+"\*SessionTable*.csv");
csvFileInterrupt    =   dir(SessionFolder+"\*InterruptTable*.csv");

if isempty(csvFileBehavior) || isempty(csvFileInterrupt)
    error("Create the SessionClass and generate Behav and Interrput table in advance.");
elseif length(csvFileBehavior)>1 || length(csvFileInterrupt)>1
    error("More than one set of Behav and Interrput table.");
end

BehTableName        =   csvFileBehavior(1).name;
BehTableFile        =   fullfile(SessionFolder, BehTableName);
BehTable            =   readtable(BehTableFile);

IntTableName        =   csvFileInterrupt(1).name;
IntTableFile        =   fullfile(SessionFolder, IntTableName);
IntTable            =   readtable(IntTableFile);

tMarkingEvent       =   BehTable.(MarkingEvent)  + BehTable.TrialStartTime; % in seconds
tMarkingLag         =   BehTable.(MarkingEvent)  - BehTable.(LagEvent) - 1; % in second
tCheckingEvent      =   BehTable.(CheckingEvent) + BehTable.TrialStartTime; % in seconds

%%  ROI extraction
% Define mask
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
[Mask, FigMask] = ExtractMask(fullfile(VideoFolderTop, vidFilesTop{1}), [3000 5000]);
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

%% Alignment: tsROI and MarkingEvent 
% use threshold method to find LED onset
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
[tLEDon, tLEDoff, FigLEDon] = FindLEDonoffGPS(tsROI, SummedROI);
SaveNameFigLEDon = fullfile(VideoFolder, "LED_On");
print(FigLEDon, '-dpng', SaveNameFigLEDon);

%% 
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------\n");
% [IndOut, FigAlign] = findseqmatchrev(tMarkingEvent*1000, tLEDon, 1);
% SaveNameFigAlign = fullfile(VideoFolder, "Alignment");
% print(FigAlign, '-dpng', SaveNameFigAlign);

[IndOut, IndValid] = decodeBitMarker(tLEDon, tLEDoff, 'fs', vidObj.FrameRate);
% these LEDon times are the ones that cannot be matched to trigger. It must be a false positive signal that was picked up by mistake in "tLEDon = FindLEDon(tsROI, SummedROI);"
tLEDon = tLEDon(IndValid); % remove them
% IndOut(isnan(IndOut)) = []; % at this point, each LEDon time can be mapped to a trigger time in b (Indout)

FrameInfo.tLEDon = tLEDon;
FrameInfo.Indout = IndOut;
FrameInfo.tLEDonBpod = tMarkingEvent(IndOut); % LED timing in Bpod's world

%% Now, let's redefine the frame time. Each frame time should be re-mapped to the timespace in bpod.
% all frame times are here: FrameInfo.tframe
tFrames2Bpodms = MapVidFrameTime2Bpod(tLEDon, tMarkingEvent(IndOut)*1000, tMarkingLag(IndOut)*1000, tsROI);
FrameTable = table(tsROI, tFrames2Bpodms, SummedROI, AviFileIndx, AviFrameIndx, MyVidFiles);

SaveNameFrameTable = fullfile(VideoFolder, "FrameInfo.csv");
writetable(FrameTable, SaveNameFrameTable);

FrameInfo.tFramesInBpod = tFrames2Bpodms;

save(SaveNameFrameInfo, "FrameInfo");

%%
if view_front
    tsROIFront        = [];
    AviFileIndxFront  = [];
    AviFrameIndxFront = [];
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
            tsROIFront          =   [tsROIFront; FrameTs(k)];
            AviFileIndxFront    =   [AviFileIndxFront; i];
            AviFrameIndxFront   =   [AviFrameIndxFront; k];
        end
    end
    tsROIFront  = tsROIFront - tsStart;

    tFrames2BpodmsFront = nan(length(tsROIFront), 1);
    for i = 1:length(tsROIFront)
        [~, id] = min(abs(tsROIFront(i) - tsROI));
        tFrames2BpodmsFront(i) = tFrames2Bpodms(id) + tsROIFront(i) - tsROI(id);
    end
    clear i id

    FrameInfoFront              =   struct();
    FrameInfoFront.tframe       =   tsROIFront;
    FrameInfoFront.AviFileIndx  =   AviFileIndxFront;
    FrameInfoFront.AviFrameIndx =   AviFrameIndxFront;

    SaveNameFrameInfoFront      =   fullfile(VideoFolder, "FrameInfoFront.mat");
    save(SaveNameFrameInfoFront, "FrameInfoFront");

    MyVidFilesFront = cell(length(AviFileIndxFront), 1);
    for i = 1:length(MyVidFilesFront)
        MyVidFilesFront{i} = vidFilesFront{AviFileIndxFront(i)};
    end
    clear i

    FrameTableFront = table(tsROIFront, tFrames2BpodmsFront, AviFileIndxFront, AviFrameIndxFront, MyVidFilesFront, ...
        'VariableNames', {'tsROI', 'tFrames2Bpodms', 'AviFileIndx', 'AviFrameIndx', 'MyVidFiles'});
    [~, ia] = unique(FrameTableFront.tsROI);
    FrameTableFront = FrameTableFront(ia, :);

    SaveNameFrameTableFront = fullfile(VideoFolder, "FrameInfoFront.csv");
    writetable(FrameTableFront, SaveNameFrameTableFront);
end

%%
if view_field
    tsROIField        = [];
    AviFileIndxField  = [];
    AviFrameIndxField = [];

    FrameTs           = [];
    FrameIdx          = [];
    for i = 1:length(tsFilesField)

        fileID          =   fopen(fullfile(VideoFolderField, tsFilesField{i}), 'r');
        formatSpec      =   '%f' ;
        NumOuts         =   fscanf(fileID, formatSpec); % this contains frame time (in ms) and frame index
        fclose(fileID);

        ind_brk         =   find(NumOuts==0);
        FrameTs         =   [FrameTs; NumOuts(1:ind_brk-1)];   % frame timestamps
        FrameIdx        =   [FrameIdx; NumOuts(ind_brk+1:end)]; % frame index
    end

    j = 0;
    for i = 1:length(vidFilesField)
        
        filename = fullfile(VideoFolderField, vidFilesField{i});
        vidObj   = VideoReader(filename);
        
        for k = 1:vidObj.NumFrames
            j = j + 1;
            tsROIField        = [tsROIField; FrameTs(j)];
            AviFileIndxField  = [AviFileIndxField; i];
            AviFrameIndxField = [AviFrameIndxField; k];
        end
    end
    clear i
    tsROIField  = tsROIField - tsStart;

    tFrames2BpodmsField = nan(length(tsROIField), 1);
    for i = 1:length(tsROIField)
        [~, id] = min(abs(tsROIField(i) - tsROI));
        tFrames2BpodmsField(i) = tFrames2Bpodms(id) + tsROIField(i) - tsROI(id);
    end
    clear i id

    FrameInfoField              =   struct();
    FrameInfoField.tframe       =   tsROIField;
    FrameInfoField.AviFileIndx  =   AviFileIndxField;
    FrameInfoField.AviFrameIndx =   AviFrameIndxField;

    SaveNameFrameInfoField      =   fullfile(VideoFolder, "FrameInfoField.mat");
    save(SaveNameFrameInfoField, "FrameInfoField");

    MyVidFilesField = cell(length(AviFileIndxField), 1);
    for i = 1:length(MyVidFilesField)
        MyVidFilesField{i} = vidFilesField{AviFileIndxField(i)};
    end
    clear i

    FrameTableField = table(tsROIField, tFrames2BpodmsField, AviFileIndxField, AviFrameIndxField, MyVidFilesField, ...
        'VariableNames', {'tsROI', 'tFrames2Bpodms', 'AviFileIndx', 'AviFrameIndx', 'MyVidFiles'});
    [~, ia] = unique(FrameTableField.tsROI);
    FrameTableField = FrameTableField(ia, :);

    SaveNameFrameTableField = fullfile(VideoFolder, "FrameInfoField.csv");
    writetable(FrameTableField, SaveNameFrameTableField);
end

%%
if view_init
    tsROIInit        = [];
    AviFileIndxInit  = [];
    AviFrameIndxInit = [];

    FrameTs           = [];
    FrameIdx          = [];
    for i = 1:length(tsFilesInit)

        fileID          =   fopen(fullfile(VideoFolderInit, tsFilesInit{i}), 'r');
        formatSpec      =   '%f' ;
        NumOuts         =   fscanf(fileID, formatSpec); % this contains frame time (in ms) and frame index
        fclose(fileID);

        ind_brk         =   find(NumOuts==0);
        FrameTs         =   [FrameTs; NumOuts(1:ind_brk-1)];   % frame timestamps
        FrameIdx        =   [FrameIdx; NumOuts(ind_brk+1:end)]; % frame index
    end

    j = 0;
    for i = 1:length(vidFilesInit)
        
        filename = fullfile(VideoFolderInit, vidFilesInit{i});
        vidObj   = VideoReader(filename);
        
        for k = 1:vidObj.NumFrames
            j = j + 1;
            tsROIInit        = [tsROIInit; FrameTs(j)];
            AviFileIndxInit  = [AviFileIndxInit; i];
            AviFrameIndxInit = [AviFrameIndxInit; k];
        end
    end
    clear i
    tsROIInit  = tsROIInit - tsStart;

    tFrames2BpodmsInit_org = MapVidFrameTime2Bpod(tLEDon, tMarkingEvent(IndOut)*1000, zeros(length(IndOut), 1), tsROI);
    tFrames2BpodmsInit = nan(length(tsROIInit), 1);
    for i = 1:length(tsROIInit)
        [~, id] = min(abs(tsROIInit(i) - tsROI));
        tFrames2BpodmsInit(i) = tFrames2BpodmsInit_org(id) + tsROIInit(i) - tsROI(id);
    end
    clear i id

    FrameInfoInit              =   struct();
    FrameInfoInit.tframe       =   tsROIInit;
    FrameInfoInit.AviFileIndx  =   AviFileIndxInit;
    FrameInfoInit.AviFrameIndx =   AviFrameIndxInit;

    SaveNameFrameInfoInit      =   fullfile(VideoFolder, "FrameInfoInit.mat");
    save(SaveNameFrameInfoInit, "FrameInfoInit");

    MyVidFilesInit = cell(length(AviFileIndxInit), 1);
    for i = 1:length(MyVidFilesInit)
        MyVidFilesInit{i} = vidFilesInit{AviFileIndxInit(i)};
    end
    clear i

    FrameTableInit = table(tsROIInit, tFrames2BpodmsInit, AviFileIndxInit, AviFrameIndxInit, MyVidFilesInit, ...
        'VariableNames', {'tsROI', 'tFrames2Bpodms', 'AviFileIndx', 'AviFrameIndx', 'MyVidFiles'});
    [~, ia] = unique(FrameTableInit.tsROI);
    FrameTableInit = FrameTableInit(ia, :);

    SaveNameFrameTableInit = fullfile(VideoFolder, "FrameInfoInit.csv");
    writetable(FrameTableInit, SaveNameFrameTableInit);
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
                switch view_init
                    case 1
                        save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'FrameTableFront', 'FrameTableField', 'FrameTableInit', 'SessionInfo', 'ClipInfo', 'ClipInfoField', 'ClipInfoInit');
                    case 0
                        save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'FrameTableFront', 'FrameTableField', 'SessionInfo', 'ClipInfo', 'ClipInfoField');
                end
            case 0
                save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'FrameTableFront', 'SessionInfo', 'ClipInfo');

        end
    case 0
        save(SaveNameUsedData, 'BehTable', 'IntTable', 'FrameTable', 'SessionInfo', 'ClipInfo');
end

%%
clearvars -except SaveNameUsedData ScnScl;
load(SaveNameUsedData);

%%
fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");

ExportVideoClipFromAvi(BehTable, IntTable, FrameTable, SessionInfo, ClipInfo, ScnScl, 0, "Top", 0);

% ExportVideoClipFromAvi(BehTable, IntTable, FrameTableFront, SessionInfo, ClipInfo, ScnScale, 0, "Front", 1);
% 
% ExportVideoClipFromAvi_Field(BehTable, FrameTableField, SessionInfo, ClipInfoField, ScnScl, 0, 0);
% ExportVideoClipFromAvi_Init(BehTable, FrameTableInit, SessionInfo, ClipInfoInit, ScnScl, 0, 0);

fprintf("\n----------------------------------------");
fprintf("\n----------------------------------------");
fprintf("\nDone.\n\n");

