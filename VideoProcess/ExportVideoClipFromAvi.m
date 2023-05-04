function ExportVideoClipFromAvi(BehTable, IntTable, FrameTable, SessionInfo, ClipInfo, remake)

% ExportVideoClipFromAvi(thisTable, FrameTable, 'Event', VideoEvent,'ANM', ANM, ...
%     'Pre', Pre, 'Post', Post, 'BehaviorType', BehaviorType, 'Session', Session, 'Remake', 1)
% Jianing Yu

% 5/1/2021
% 4/16/2022

event       = ClipInfo.VideoEvent;
tPre        = ClipInfo.Pre;
tPost       = ClipInfo.Post;
anm         = SessionInfo.ANM;
beh_type    = SessionInfo.Task;
session     = SessionInfo.Session;

if nargin<6
    remake = 0;
end

%  in mili-seconds, timing of selected events.
tBehEvent = 1000*BehTable.(event);
% in mili-seconds, timing of each frame in Bpod's world
tFramesBpod = FrameTable.tFrames2Bpodms;

% set up video clip storage folder
thisFolder = fullfile(ClipInfo.VideoFolder, 'Clips');
if ~exist(thisFolder, 'dir')
    mkdir(thisFolder);
end

% setting up metadata
VidsMeta = struct('Session', [], 'Event', [], 'EventIndex', [], 'Performance', [], 'EventTime', [], 'FrameTimesB', [], 'VideoOrg', [], 'FrameIndx', [], 'Code', [], 'CreatedOn', []);

%% Start making videos
clc
video_acc = 0;

for i =1:length(tBehEvent) % i is also the trial number

    itEvent = tBehEvent(i);

    IndThisClip         =   find(tFramesBpod>=itEvent - tPre & tFramesBpod<= itEvent + tPost);
    [~, IndThisFrame]   =   min(abs(tFramesBpod - itEvent));

    % check if a video has been created and check if we want to
    % re-create the same video
    switch event
        case 'PortSamplePokeTime'
            ClipName = sprintf('%s_%s_SamplePokeTrial%03d', anm, session, i);
            ievent = 'SamplePoke';

        case 'PortCenterPokeTime'
            ClipName = sprintf('%s_%s_CenterPokeTrial%03d', anm, session, i);
            ievent = 'CenterPoke';
    end

    VidClipFileName = fullfile(thisFolder, [ClipName '.avi']);
    check_this_file = dir(VidClipFileName);

    if ~isempty(check_this_file)  && ~remake % found a video clip with the same name, and we don't want to remake the video clip
        continue % move on
    end


    iFrameTimesBpod = tFramesBpod(IndThisClip);
    % make sure the videoclip can be constructed from a single video file
    if itEvent -iFrameTimesBpod(1) < tPre-50
        continue
    elseif  iFrameTimesBpod(end)-itEvent < tPost-50
        continue
    elseif ~strcmp(FrameTable.MyVidFiles{(IndThisClip(1))}, FrameTable.MyVidFiles{(IndThisFrame)}) || ~strcmp(FrameTable.MyVidFiles{(IndThisClip(end))}, FrameTable.MyVidFiles{(IndThisFrame)})
        % same video file
        continue
    end

    % a few important events
    SamplePokeTime = BehTable.PortSamplePokeTime(i)*1000-iFrameTimesBpod(1);
    CenterPokeTime = BehTable.PortCenterPokeTime(i)*1000-iFrameTimesBpod(1);
    ChoicePokeTime = BehTable.ChoicePortTime(i)*1000-iFrameTimesBpod(1);

    % build video clips, frame by frame
    F= struct('cdata', [], 'colormap', []);

    VidMeta.Session               =          session;
    VidMeta.Event                   =          event;
    VidMeta.EventIndex          =         i;
    VidMeta.EventTime           =        itEvent/1000;       % Event time in sec (Bpod)
    VidMeta.FrameTimesB     =          tFramesBpod(IndThisClip);                        % frame time in ms in behavior time
    VidMeta.FrameIndx          =           IndThisClip;                   % frame index in original video
    VidMeta.Code                   =             mfilename('fullpath');
    VidMeta.CreatedOn          =            date;                                % today's date

    video_acc = video_acc+1;
    if video_acc ==1
        VidsMeta = VidMeta;
    else
        VidsMeta(video_acc) = VidMeta;
    end

    % Extract frames
    VidFrameIdx_thisfile = FrameTable.AviFrameIdx(IndThisClip);  % these are the frame index in this video
    this_video = FrameTable.MyVidFiles{IndThisFrame};
    vidObj = VideoReader(this_video);
    img_extracted = [];
    frames_ifile = read(vidObj, [VidFrameIdx_thisfile(1) VidFrameIdx_thisfile(end)]);
    for ii =1:size(frames_ifile, 4)
        img_extracted = cat(3, img_extracted, rgb2gray(frames_ifile(:, :, :, ii)));
    end

    [height, width, nframe] = size(img_extracted); %

    height = height - 200;

    % Make videos

    for k =1:nframe

        hf25 = figure(25); clf
        set(hf25, 'name', 'side view', 'units', 'centimeters', 'position', [ 3 5 15 3+15*(height-100)/width], 'PaperPositionMode', 'auto', 'color', 'w')

        ha= axes;
        set(ha, 'units', 'centimeters', 'position', [0 3 15 15*(height-100)/width], 'nextplot', 'add', 'xlim',[0 width], 'ylim', [100 height], 'ydir','reverse')
        axis off

        % plot this frame:
        imagesc(img_extracted(:, :, k), [0 250]);
        colormap('gray')

        % plot some behavior data
        tthis_frame =  round(iFrameTimesBpod(k) - iFrameTimesBpod(1));
        time_of_frame = sprintf('%3.0f',tthis_frame);

        text(10, height-150,  sprintf('%s %s',anm, session), 'color', [255 255 255]/255, 'fontsize',  8, 'fontweight', 'bold')
        text(10, height-110,  sprintf('%ss',beh_type), 'color', [255 255 255]/255, 'fontsize',  8, 'fontweight', 'bold')
        text(10, height-70,  sprintf('%s %03d', 'Trial#', i), 'color', [255 255 255]/255, 'fontsize',  8, 'fontweight', 'bold')
        text(10, height-30,  sprintf('PortCorrect:%2.0d; PortChosen:%2.0d', BehTable.CorrectPort(i), BehTable.PortChosen(i) ), 'color', [255 255 255]/255, 'fontsize',  8, 'fontweight', 'bold')
        text(10, 130, [time_of_frame ' ms'], 'color', [255 255 255]/255, 'fontsize', 12,'fontweight', 'bold')
        % plot some important behavioral events
        ha2= axes;
        set(ha2, 'units', 'centimeters', 'position', [1.5 0.75 13 2], ...
            'nextplot', 'add', 'xtick', [0:1000:tPost+tPre], 'xlim', [0 iFrameTimesBpod(end)-iFrameTimesBpod(1)], 'ytick', [0.5 1.5 2.5], ...
            'yticklabel', {'Sample', 'Center', 'Choice'},'ylim', [0 3], 'tickdir', 'out') %#ok<NBRAK>

        line([SamplePokeTime SamplePokeTime], [0 1], 'linewidth', 4, 'color', 'k')
        line([CenterPokeTime CenterPokeTime], [1 2], 'linewidth', 4, 'color', 'c')
        if BehTable.CorrectPort(i)== BehTable.PortChosen(i)
            line([ChoicePokeTime ChoicePokeTime], [2 3], 'linewidth', 4, 'color', 'g')
        else
            line([ChoicePokeTime ChoicePokeTime], [2 3], 'linewidth', 4, 'color', 'r')
        end

        % plot current time
        line([tthis_frame tthis_frame], [0 4], 'linestyle', '-', 'linewidth', 1)

        % plot or update data in this plot
        F(k) = getframe(hf25) ;

    end
    % make a video clip and save it to the correct location

    writerObj = VideoWriter([ClipName '.avi']);
    writerObj.FrameRate = median(1000./diff(iFrameTimesBpod)); % this is 1 x
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for ifrm=1:length(F)
        % convert the image to a frame
        frame = F(ifrm) ;
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    movefile( [ClipName '.avi'], thisFolder)

    MetaFileName = fullfile(thisFolder, [ClipName, '.mat']);
    save(MetaFileName, 'VidMeta');
end
toc
