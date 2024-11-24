function ExportVideoClipFromAvi_Field(BehTable, FrameTable, SessionInfo, ClipInfo, scn_scale, remake, x_rev)

% ExportVideoClipFromAvi(thisTable, FrameTable, 'Event', VideoEvent,'ANM', ANM, ...
%     'Pre', Pre, 'Post', Post, 'BehaviorType', BehaviorType, 'Session', Session, 'Remake', 1)
% Jianing Yu

% 5/1/2021
% 4/16/2022

% revised by ZZH, 5/5/2023

scale_ratio     =   1 / scn_scale;

color           =   GPSColor();

anm             =   SessionInfo.ANM;
beh_type        =   SessionInfo.Task;
session         =   SessionInfo.Session;

event           =   ClipInfo.VideoEvent;
tPreMax         =   ClipInfo.Pre;
tPost           =   ClipInfo.Post;

if nargin<6
    remake = 0;
elseif nargin < 7
    x_rev  = 0;
end

%  in mili-seconds, timing of selected events.
tBehEvent = 1000*(BehTable.(event) + BehTable.TrialStartTime);
% in mili-seconds, timing of each frame in Bpod's world
tFramesBpod = FrameTable.tFrames2Bpodms;

% set up video clip storage folder
thisView   = "Field";
viewFolder = ClipInfo.VideoFolderField;
thisFolder = fullfile(ClipInfo.VideoFolderField, 'Clips');

if ~isfolder(thisFolder)
    mkdir(thisFolder);
end

% setting up metadata
VidsMeta = struct('Session', [], 'Event', [], 'EventIndex', [], 'Performance', [], 'EventTime', [], 'FrameTimesB', [], 'VideoOrg', [], 'FrameIndx', [], 'Code', [], 'CreatedOn', []);

%% Start making videos
fprintf("\nStart video clipping ...\n")
video_accum = 0;

wait_bar = waitbar(0, sprintf('1 / %d', length(tBehEvent)), 'Name', sprintf('Clipping_%s_%s', anm, session));
for i = 1:length(tBehEvent) % i is also the trial number
    
    if ~isvalid(wait_bar)
        fprintf("\n****** Interrupted ******\n");
        fprintf("%d / %d clips have been generated\n", i-1, length(tBehEvent));
        return
    end
    waitbar(i/length(tBehEvent), wait_bar, sprintf('%d / %d', i, length(tBehEvent)));

    itEvent = tBehEvent(i);
    iShuttleTime = BehTable.ST(i); % Shuttle time of this trial (in sec)
    if (500 + iShuttleTime*1000) < tPreMax
        tPre = iShuttleTime*1000;
    else
        tPre = tPreMax;
    end

    IndThisClip         =   find(tFramesBpod>=itEvent-tPre & tFramesBpod<=itEvent+tPost);
    if isempty(IndThisClip)
        continue
    end
    [~, IndThisFrame] = min(abs(tFramesBpod - itEvent));

    % check if a video has been created and check if we want to
    % re-create the same video
    ClipName = sprintf('%s_%s_Trial%03d_FieldView', anm, session, i);

    VidClipFileName = fullfile(thisFolder, [ClipName '.avi']);
    check_this_file = dir(VidClipFileName);

    if ~isempty(check_this_file) && ~remake % found a video clip with the same name, and we don't want to remake the video clip
        continue % move on
    end

    iFrameTimesBpod = tFramesBpod(IndThisClip);
    % make sure the videoclip can be constructed from a single video file
    if itEvent-iFrameTimesBpod(1) < tPre-50
        continue
    elseif iFrameTimesBpod(end)-itEvent < tPost-50
        continue
    elseif ~strcmp(FrameTable.MyVidFiles{(IndThisClip(1))}, FrameTable.MyVidFiles{(IndThisFrame)}) || ~strcmp(FrameTable.MyVidFiles{(IndThisClip(end))}, FrameTable.MyVidFiles{(IndThisFrame)})
        % same video file
        continue
    end

    % build video clips, frame by frame
    F = struct('cdata', [], 'colormap', []);

    VidMeta.Session     =   session;
    VidMeta.Event       =   event;
    VidMeta.EventIndex  =   i;
    VidMeta.EventTime   =   itEvent/1000; % Event time in sec (Bpod)
    VidMeta.FrameTimesB =   tFramesBpod(IndThisClip); % frame time in ms in behavior time
    VidMeta.FrameIndx   =   IndThisClip; % frame index in original video
    VidMeta.Code        =   mfilename('fullpath');
    VidMeta.CreatedOn   =   date; % today's date

    video_accum = video_accum + 1;
    if video_accum==1
        VidsMeta = VidMeta;
    else
        VidsMeta(video_accum) = VidMeta;
    end

    % Extract frames
    VidFrameIndx_thisfile = FrameTable.AviFrameIndx(IndThisClip);  % these are the frame index in this video
    this_video = fullfile(viewFolder, FrameTable.MyVidFiles{IndThisFrame});
    vidObj = VideoReader(this_video);
    img_extracted = [];
    for ii = 1:length(VidFrameIndx_thisfile)
        img_extracted = cat(3, img_extracted, rgb2gray(read(vidObj, VidFrameIndx_thisfile(ii))));
    end
    clear frames_ifile vidObj
% 
%     if size(img_extracted, 1)==1240
%         img_extracted = img_extracted(101:1124, 201:1480, :);
%     end
    [H, W, nframe] = size(img_extracted); %

    % height = height - 200;

    %% Make videos

    k = 1;
    hf25 = figure(25); clf
    set(hf25, 'name', thisView, 'units', 'pixels', 'position', [5 50 scale_ratio*W scale_ratio*H], ...
        'PaperPositionMode', 'auto', 'color', 'w', 'renderer', 'opengl', 'toolbar', 'none', 'resize', 'off', 'Visible', 'on');

    ha = axes;
    set(ha, 'units', 'pixels', 'position', [0 0 scale_ratio*W scale_ratio*H], 'nextplot', 'add', 'xlim', [.5 W+.5], 'ylim', [.5 H+.5], 'ydir', 'reverse')
    if x_rev==1
        set(ha, 'xdir', 'reverse');
    end
    axis off

    % plot this frame:
    img = imagesc(ha, img_extracted(:, :, k), [0 250]);
    colormap('gray');

    % plot some behavior data
    tthis_frame   = round(iFrameTimesBpod(k) - iFrameTimesBpod(1) - tPre);
    time_of_frame = sprintf('%3.0f', tthis_frame);

    text(W-20, 40,  sprintf('%s %s', anm, session), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 90,  beh_type, 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 140,  sprintf('Trial %03d', BehTable.Trials(i)), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    time_text = text(20, 40, [time_of_frame ' ms'], 'color', [255 215 0]/255, 'FontSize', 22,'fontweight', 'bold');
    % plot some important behavioral events

    F(k) = getframe(hf25);
    % plot or update data in this plot
    for k = 2:nframe

        tthis_frame = round(iFrameTimesBpod(k) - iFrameTimesBpod(1) - tPre);
        time_of_frame = sprintf('%3.0f', tthis_frame);
        time_text.String = [time_of_frame ' ms'];

        time_line.Value = tthis_frame;

        img.CData = img_extracted(:, :, k);

%         drawnow;
        % plot or update data in this plot
        F(k) = getframe(hf25);
    end
    % make a video clip and save it to the correct location
    close(hf25);
    clear img_extracted

    writerObj = VideoWriter(VidClipFileName);
    writerObj.FrameRate = 0.5 * median(1000./diff(iFrameTimesBpod));
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for ifrm = 1:length(F)
        % convert the image to a frame
        frame = F(ifrm);
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    clear writerObj F IndThisClip IndThisFrame

    MetaFileName = fullfile(thisFolder, [ClipName, '.mat']);
    save(MetaFileName, 'VidMeta');

end

close(wait_bar);

