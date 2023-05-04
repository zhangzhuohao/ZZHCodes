%%  File information
% Use this func to find video files and associated timestamp files.
% go to folder that includes video files, use 'ShowAviFiles' to get vid
% file names and folder name

% run ShowAviFiles.m, select any vid file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Fill in the following information %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

fprintf("\nGPS video clipping ...\n");

VideoFolder     =   uigetdir("E:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolder
    return;
end
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
        end
    end
end
clear i entry

if ~view_top
    fprintf("\nFind no top-view vedios.\n");
    return
else
    fprintf("\nTop-view vedios folder: \n%s\n", VideoFolderTop);
    [vidFilesTop, tsFilesTop, ~] = ShowAviFiles(VideoFolderTop);
    fprintf(" Videos\t\tTimestamps\n");
    cellfun(@(vid, ts) fprintf(" %s\t%s\n", vid, ts), vidFilesTop, tsFilesTop);
end

if ~view_front
    fprintf("\nFind no front-view vedios.\n");
else
    fprintf("\nFront-view vedios folder: \n%s\n", VideoFolderFront);
    [vidFilesFront, tsFilesFront, ~] = ShowAviFiles(VideoFolderFront);
    fprintf(" Videos\t\tTimestamps\n");
    cellfun(@(vid, ts) fprintf(" %s\t%s\n", vid, ts), vidFilesFront, tsFilesFront);
end

%%
vidInfo         =   split(VideoFolder, filesep);
ANM             =   string(vidInfo{end-1});
Session         =   string(vidInfo{end});

ANMInfoFile     =   "E:\YuLab\Work\GPS\Data\ANMInfo.xlsx";
ANMInfo         =   readtable(ANMInfoFile, "Sheet", ANM, "TextType", "string");
ANMInfo.Session =   string(ANMInfo.Session);

if isempty(ANMInfo.Session==Session)
    fprintf("\nCannot find session information.\n")
    return
elseif ANMInfo(ANMInfo.Session==Session, :).Task == "None"
    fprintf("\nThe session class has not been generated.\n")
    return
else
    SessionInfo = ANMInfo(ANMInfo.Session==Session, :);
    disp(SessionInfo(:, 1:6));
end

CSVFile         =   dir(SessionInfo.SessionFolder+"\*.csv");
BehTableName    =   CSVFile(1).name;
BehTableInfo    =   split(BehTableName, '_');

%%
vidfolder           =           pwd;
CSVFile             =           dir('*.csv');
BehTableName        =           CSVFile(1).name;
BehTableInfo        =           split(BehTableName, '_');
ANM                 =           BehTableInfo{1};
Session             =           BehTableInfo{5}(1:end-4);
BehaviorType        =           BehTableInfo{4};
thisTable           =           fullfile(pwd, BehTableName);
MarkingEvent        =           'PortCenterPokeTime'; % use this time to align bpod and video ts
CheckingEvent       =           'CenterLightTime'; % use this time to align bpod and video ts
VideoEvent          =           'PortCenterOutTime'; % make video clips based on this event
Pre                 =           1000; % ms pre event time for video clips
Post                =           2500; % ms post event time for video clips

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if isempty(thisTable)
    [thisTable, path] = uigetfile('.csv', 'Select a behavior table');
    thisTable = fullfile(path, thisTable);
end

thisTable = readtable(thisTable);
tMarkingEvent = eval(['thisTable.' MarkingEvent]); % in seconds
tCheckingEvent = eval(['thisTable.' CheckingEvent]); % in seconds

%%  ROI extraction
% Define mask 
mask = ExtractMask(fullfile(vidfolder, vidFiles{1}), [9000 10000]);

%% based on "mask", extract pixel intensity in ROI from all frames
tsROI = [];
SummedROI = [];
AviFrameIndx = [];
AviFileIndx = [];
tic

clc
sprintf('Extracting ......')


for i=1:length(vidFiles)

    fileID          =       fopen(fullfile(vidfolder, tsFiles{i}), 'r');
    formatSpec      =       '%f' ;
    NumOuts         =       fscanf(fileID, formatSpec); % this contains frame time (in ms) and frame index
    fclose(fileID);

    ind_brk = find(NumOuts ==0);
    FrameTs = NumOuts(1:ind_brk-1);  % frame times

    FrameIdx = NumOuts(ind_brk+1:end);                     % frame idx
    filename = fullfile(vidfolder, vidFiles{i});

    vidObj = VideoReader(filename);
    parfor k =1:  vidObj.NumberOfFrames
        thisFrame = read(vidObj, [k k]);
        thisFrame = thisFrame(:, :, 1);
        roi_k = sum(thisFrame(mask));

        tsROI           =       [tsROI FrameTs(k)];
        SummedROI       =       [SummedROI roi_k];
        AviFileIndx     =       [AviFileIndx i];
        AviFrameIndx    =       [AviFrameIndx k];

    end
    toc
end


tsROI = tsROI - tsROI(1); % onset normalized to 0
FraeInfo                          =   [];
FrameInfo.tframe                   =     tsROI;
FrameInfo.mask                      =    mask;
FrameInfo.ROI                         =    SummedROI;
FrameInfo.AviFile                   =    vidFiles;
FrameInfo.AviFileIndx           =     AviFileIndx;
FrameInfo.AviFrameIndx      =      AviFrameIndx;

% Save for now because it takes a long time to get tsROI
save FrameInfo FrameInfo

tsROI                   =       tsROI';
SummedROI       =       SummedROI';
AviFileIndx          =       AviFileIndx';
AviFrameIdx        =       AviFrameIndx';

MyVidFiles          =       cell(length(AviFileIndx), 1);
for i =1:length(MyVidFiles)
    MyVidFiles{i} = vidFiles{AviFileIndx(i)};
end
% 
% aTable = table(tsROI, SummedROI, AviFileIndx, AviFrameIdx, MyVidFiles);
% 
% aGoodTableName =  ['FrameInfo' '.csv'];
% writetable(aTable, aGoodTableName);

%% Alignment: tsROI and MarkingEvent 
% use threshold method to find LED onset
tLEDon = FindLEDonGSP(tsROI, SummedROI);
Indout = findseqmatchrev(tMarkingEvent*1000, tLEDon, 0, 1);
% these LEDon times are the ones that cannot be matched to trigger. It must be a false positive signal that was picked up by mistake in "tLEDon = FindLEDon(tsROI, SummedROI);"
tLEDon(isnan(Indout)) = []; % remove them
Indout(isnan(Indout)) = []; % at this point, each LEDon time can be mapped to a trigger time in b (Indout)

FrameInfo.tLEDon = tLEDon;
FrameInfo.Indout = Indout;
FrameInfo.tLEDonBpod = tMarkingEvent(Indout); % LED timing in Bpod's world

%% Now, let's redefine the frame time. Each frame time should be re-mapped to the timespace in bpod.
% all frame times are here: FrameInfo.tframe
tFrames2Bpodms = MapVidFrameTime2Bpod(tLEDon,  tMarkingEvent(Indout)*1000, tsROI); 
FrameTable = table(tsROI, tFrames2Bpodms, SummedROI, AviFileIndx, AviFrameIdx, MyVidFiles);
FrameTableName =  ['FrameInfo' '.csv'];
writetable(FrameTable, FrameTableName);
FrameInfo.tFramesInBpod = tFrames2Bpodms;

%% Check if another event is also aligned
% Previously, we have defined checking event:
% CheckingEvent = 'CenterLightTime'; % use this time to align bpod and video ts
% tCheckingEvent = eval(['thisTable.' CheckingEvent]); % in seconds
CheckAnotherEvent = 1;
if CheckAnotherEvent% still working on this
end

%% Make video clips
% Extrct           =       "ExportVideoFiles(b, FrameInfo, 'Event', 'Press', 'TimeRange', [2000 3000], 'SessionName',strrep(fileinfo.MED(1:10), '-', ''), 'RatName', MEDfile(27:strfind(MEDfile, '.')-1), 'Remake', 0)";
% MakeSheet   =      "MakeSpreadSheet(b, FrameInfo, 'Event', 'Press', 'TimeRange', [2000 3000], 'SessionName', strrep(fileinfo.MED(1:10), '-', ''), 'RatName', MEDfile(27:strfind(MEDfile, '.')-1))";

% make video clip based on two tables: a behavioral data table and a frame
% info table
% 
thisTable               =           readtable(BehTableName);
FrameTable          =           readtable('FrameInfo.csv');

save('UsedData', 'thisTable', 'FrameTable', 'VideoEvent', 'ANM', 'Pre', 'Post', 'BehaviorType', 'Session');
clear
load('UsedData.mat');
ExportVideoClipFromAvi(thisTable, FrameTable, 'Event', VideoEvent,'ANM', ANM, 'Pre', Pre, 'Post', Post, 'BehaviorType', BehaviorType, 'Session', Session, 'Remake', 0)