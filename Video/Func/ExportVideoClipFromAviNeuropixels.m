function ExportVideoClipFromAviNeuropixels(r, FrameTable, SessionInfo, ClipInfo, scn_scale, remake, view, x_rev)

% ExportVideoClipFromAvi(thisTable, FrameTable, 'Event', VideoEvent,'ANM', ANM, ...
%     'Pre', Pre, 'Post', Post, 'BehaviorType', BehaviorType, 'Session', Session, 'Remake', 1)
% Jianing Yu

% 5/1/2021
% 4/16/2022

% revised by ZZH, 5/5/2023

img_ratio = .5;

BehClass = r.BehaviorClass;
rb = r.Behavior;

scale_ratio     =   1 / scn_scale;

color           =   GPSColor();

anm             =   BehClass.Subject;
beh_type        =   SessionInfo.Task;
session         =   SessionInfo.Session;

event           =   ClipInfo.VideoEvent;
tPre            =   ClipInfo.Pre;
tPost           =   ClipInfo.Post;
time_elapsed    =   -tPre:0.1:tPost;

if nargin<6
    remake = 0;
    view   = "Top";
elseif nargin < 7
    view   = "Top";
end

%  in mili-seconds, timing of selected events.
lb_beh = find(strcmp(rb.Labels, event));
tBehEvent = rb.EventTimings(rb.EventMarkers==lb_beh); % in milliseconds
% in mili-seconds, timing of each frame in Bpod's world
tFramesPixel = FrameTable.tFrames2Pixelms;

% set up video clip storage folder
switch view
    case "Top"
        thisView   = "Top";
        viewFolder = ClipInfo.VideoFolderTop;
        thisFolder = fullfile(ClipInfo.VideoFolderTop, 'Clips');
    case "Front"
        thisView   = "Front";
        viewFolder = ClipInfo.VideoFolderFront;
        thisFolder = fullfile(ClipInfo.VideoFolderFront, 'Clips');
    case "Field"
        thisView   = "Field";
        viewFolder = ClipInfo.VideoFolderField;
        thisFolder = fullfile(ClipInfo.VideoFolderField, 'Clips');
end
if ~isfolder(thisFolder)
    mkdir(thisFolder);
end

% setting up metadata
VidsMeta = struct('Session', [], 'Event', [], 'EventIndex', [], 'Performance', [], 'EventTime', [], 'FrameTimesE', [], 'VideoOrg', [], 'FrameIndx', [], 'Code', [], 'CreatedOn', []);

%% colors
% events
c = colororder;
cCentIn = c(2,:);
cTrigger = c(3,:);
cCentOut = c(4,:);
cChoice = c(5,:);

% spikes
nUnits = length(r.Units.SpikeTimes);
cUnits = 1-varycolor(nUnits+1);

%% extract behavior events
lb_centin = find(strcmp(rb.Labels, "PokeCentIn"));
tCentIn   = rb.EventTimings(rb.EventMarkers==lb_centin); % in milliseconds

lb_centout = find(strcmp(rb.Labels, "PokeCentOut"));
tCentOut   = rb.EventTimings(rb.EventMarkers==lb_centout); % in milliseconds

lb_choice = find(strcmp(rb.Labels, "PokeChoiceIn"));
tChoice   = rb.EventTimings(rb.EventMarkers==lb_choice); % in milliseconds

lb_trigger = find(strcmp(rb.Labels, "Trigger"));
tTrigger   = rb.EventTimings(rb.EventMarkers==lb_trigger); % in milliseconds

lb_initin = find(strcmp(rb.Labels, "PokeInitIn"));
tInitIn   = rb.EventTimings(rb.EventMarkers==lb_initin); % in milliseconds

lb_initout = find(strcmp(rb.Labels, "PokeInitOut"));
tInitOut   = rb.EventTimings(rb.EventMarkers==lb_initout); % in milliseconds

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

    IndThisClip = find(tFramesPixel>=itEvent-tPre & tFramesPixel<=itEvent+tPost);
    if isempty(IndThisClip)
        continue
    end
    [~, IndThisFrame] = min(abs(tFramesPixel - itEvent));

    % check if a video has been created and check if we want to
    % re-create the same video
    switch event
        case 'PortSamplePokeTime'
            ClipName = sprintf('%s_%s_SamplePoke%03d', anm, session, i);

        case 'PortCenterPokeTime'
            ClipName = sprintf('%s_%s_CenterPoke%03d', anm, session, i);

        case {'CentInTime', 'PokeCentIn'}
            ClipName = sprintf('%s_%s_Hold%03d_%sView', anm, session, i, view);
    end

    VidClipFileName = fullfile(thisFolder, [ClipName '.avi']);
    check_this_file = dir(VidClipFileName);

    if ~isempty(check_this_file) && ~remake % found a video clip with the same name, and we don't want to remake the video clip
        continue % move on
    end

    iFrameTimesPixel = tFramesPixel(IndThisClip);
    % make sure the videoclip can be constructed from a single video file 
    if itEvent-iFrameTimesPixel(1) < tPre-50
        continue
    elseif iFrameTimesPixel(end)-itEvent < tPost-50
        continue
    elseif ~strcmp(FrameTable.MyVidFiles{(IndThisClip(1))}, FrameTable.MyVidFiles{(IndThisFrame)}) || ~strcmp(FrameTable.MyVidFiles{(IndThisClip(end))}, FrameTable.MyVidFiles{(IndThisFrame)})
        % same video file
        continue
    end
    tPre_this = tFramesPixel(IndThisFrame) - tFramesPixel(IndThisClip(1));

    % poke events
    %     SamplePokeTime  =   BehTable.PortSamplePokeTime(i)*1000-iFrameTimesBpod(1);
    [~, this_CentIn] = min(abs(tCentIn - itEvent));
    itCentIn = tCentIn(this_CentIn);

    this_InitIn = find(tInitIn<itCentIn, 1, 'last');
    if ~isempty(this_InitIn)
        itInitIn = tInitIn(this_InitIn);
    else
        itInitIn = nan;
    end
    this_InitOut = find(tInitOut<itCentIn, 1, 'last');
    if ~isempty(this_InitOut)
        itInitOut = tInitOut(this_InitOut);
    else
        itInitOut = nan;
    end
    next_InitIn = find(tInitIn>itCentIn, 1, 'first');
    if ~isempty(next_InitIn)
        nextInitIn = tInitIn(next_InitIn);
    else
        nextInitIn = nan;
    end
    next_InitOut = find(tInitOut>itCentIn, 1, 'first');
    if ~isempty(next_InitOut)
        nextInitOut = tInitOut(next_InitOut);
    else
        nextInitOut = nan;
    end

    this_CentOut = find(tCentOut>itCentIn, 1, 'first');
    if ~isempty(this_CentOut)
        itCentOut = tCentOut(this_CentOut);
    else
        itCentOut = nan;
    end

    this_Trigger = find(tTrigger>itCentIn, 1, 'first');
    if ~isempty(this_Trigger)
        itTrigger = tTrigger(this_Trigger);
        if isnan(nextInitIn)
            if itTrigger > nextInitIn
                itTrigger = nan;
            end
        end
    else
        itTrigger = nan;
    end

    this_Choice = find(tChoice>itCentIn, 1, 'first');
    if ~isempty(this_Choice)
        itChoice = tChoice(this_Choice);
        if isnan(nextInitIn)
            if itChoice > nextInitIn
                itChoice = nan;
            end
        end
    else
        itChoice = nan;
    end

    % build video clips, frame by frame
    F = struct('cdata', [], 'colormap', []);

    VidMeta.Session     =   session;
    VidMeta.Event       =   event;
    VidMeta.EventIndex  =   i;
    VidMeta.EventTime   =   itEvent/1000; % Event time in sec (Bpod)
    VidMeta.FrameTimesE =   tFramesPixel(IndThisClip); % frame time in ms in behavior time
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
    if x_rev==1
        for ii = 1:length(VidFrameIndx_thisfile)
            img_this = rgb2gray(read(vidObj, VidFrameIndx_thisfile(ii)));
            img_extracted = cat(3, img_extracted, img_this(:, end:-1:1));
        end
    else
        for ii = 1:length(VidFrameIndx_thisfile)
            img_extracted = cat(3, img_extracted, rgb2gray(read(vidObj, VidFrameIndx_thisfile(ii))));
        end
    end
    clear frames_ifile vidObj

    [H, W, nframe] = size(img_extracted); %
    H = H * img_ratio;
    W = W * img_ratio;

    % height = height - 200;

    % Make videos
    k = 1;
    hf25 = figure(25); clf
    set(hf25, 'name', thisView, 'units', 'pixels', 'position', [5 50 scale_ratio*W 1.8*scale_ratio*H], ...
        'PaperPositionMode', 'auto', 'color', 'k', 'renderer', 'opengl', 'toolbar', 'none', 'resize', 'off', 'Visible', 'on');

    ha = axes;
    set(ha, 'units', 'pixels', 'position', [0 .8*scale_ratio*H scale_ratio*W scale_ratio*H], 'nextplot', 'add', 'xlim', [.5 W+.5], 'ylim', [.5 H+.5], 'ydir', 'reverse', 'Color', 'none')
    axis off

    % plot this frame:
    img_this = imresize(img_extracted(:, :, k), img_ratio);
    img = imagesc(ha, img_this, [0 250]);

    colormap('gray');

    % plot some behavior data
    tthis_frame   = round(iFrameTimesPixel(k) - iFrameTimesPixel(1) - tPre_this);

%     text(W-20, 40,  sprintf('%s %s', anm, session), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
%     text(W-20, 90,  beh_type, 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
%     text(W-20, 160,  sprintf('Trial %03d', BehTable.Trials(i)), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
%     text(W-20, 190,  sprintf('FP: %d ms', BehTable.FP(i)*1000), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
%     text(W-20, 240,  sprintf('RT: %d ms', round(1000*BehTable.RT(i))), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
%     text(W-20, 290,  thisOutcome, 'color', color.(thisOutcome), 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    time_text = text(10, H-20, sprintf('Time: %0.2f s', tthis_frame ./ 1000), 'color', [255 215 0]/255, 'FontSize', 11,'fontweight', 'bold');
    % plot some important behavioral events

    % raster
    ha_raster = axes;
    set(ha_raster, 'units', 'pixels', 'position', [0 0.13*scale_ratio*H scale_ratio*W 0.63*scale_ratio*H], 'color', 'none', ...
        'nextplot', 'add', 'xtick', [-tPre:500:tPost], 'xlim', [-tPre tPost], 'xcolor', 'none', ...
        'ycolor', 'none', 'ydir', 'reverse', 'ylim', [.5 nUnits+.5], 'tickdir', 'out', 'FontSize', 11) %#ok<NBRAK>

    for u = 1:nUnits
        c_unit = cUnits(u, :);

        spike_timing = r.Units.SpikeTimes(u).timings - itEvent;
        spike_timing = spike_timing(spike_timing>=-tPre & spike_timing<=tPost);

        for spx = 1:length(spike_timing)
            line(ha_raster, spike_timing(spx)*[1 1], u+[-.5 .5], 'Color', c_unit, 'LineWidth', 1);
        end
    end

    % event lines
    ha_eline = axes; 
    set(ha_eline, 'units', 'pixels', 'position', [0 0.125*scale_ratio*H scale_ratio*W 0.67*scale_ratio*H], 'color', 'none', ...
        'nextplot', 'add', 'xtick', -tPre:500:tPost, 'xlim', [-tPre tPost], 'xcolor', 'none', ...
        'ycolor', 'none', 'ylim', [0 1], 'FontSize', 10) % #ok<NBRAK>

    xline(ha_eline, itCentIn - itCentIn, 'Color', cCentIn, 'LineStyle', '-', 'LineWidth', 1.5, 'Alpha', 1);
    xline(ha_eline, itCentOut - itCentIn, 'Color', cCentOut, 'LineStyle', '-', 'LineWidth', 1.5, 'Alpha', 1);
    xline(ha_eline, itTrigger - itCentIn, 'Color', cTrigger, 'LineStyle', '-', 'LineWidth', 1.5, 'Alpha', 1);
    xline(ha_eline, itChoice - itCentIn, 'Color', cChoice, 'LineStyle', '-', 'LineWidth', 1.5, 'Alpha', 1);
   
    text(ha_eline, itCentIn-itCentIn-20, 1, 'Cent-In', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(ha_eline, itTrigger-itCentIn-20, 1, 'Trigger', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(ha_eline, itCentOut-itCentIn+20, 1, 'Cent-Out', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(ha_eline, itChoice-itCentIn+20, 1, 'Choice-In', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    ha_tline = axes;
    set(ha_tline, 'units', 'pixels', 'position', [0 0.01*scale_ratio*H scale_ratio*W 0.785*scale_ratio*H], 'color', 'none', ...
        'nextplot', 'add', 'xtick', -tPre:500:tPost, 'xlim', [-tPre tPost], 'xcolor', 'none', ...
        'ycolor', 'none', 'ylim', [0 1]) % #ok<NBRAK>
    time_line = xline(ha_tline, tthis_frame, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5, 'Alpha', 1);

    % x axis
    x_present = [-800 2800];
    x_ratio = range(x_present) ./ (tPost+tPre);
    x_width = x_ratio * scale_ratio * W;
    x_height = 0.01*scale_ratio*H;
    x_pos_x = scale_ratio * W * (tPre+x_present(1)) ./ (tPost+tPre);
    x_pos_y = 0.12*scale_ratio*H;

    ha_x = axes;
    set(ha_x, 'units', 'pixels', 'position', [x_pos_x x_pos_y x_width x_height], 'color', 'none', ...
        'nextplot', 'add', 'xtick', (-tPre:500:tPost)./1000, 'xlim', x_present./1000, 'xcolor', 'w', ...
        'ycolor', 'none', 'ydir', 'reverse', 'ylim', [0 1], 'tickdir', 'out', 'FontSize', 10, 'LineWidth', 1) % #ok<NBRAK>
    ha_x.XLabel.String = 'Time from Cent-In (s)';
    ha_x.XLabel.FontWeight = 'bold';
    ha_x.XLabel.FontSize = 10;

%
    F(k) = getframe(hf25);
    % plot or update data in this plot
    for k = 2:nframe

        tthis_frame = round(iFrameTimesPixel(k) - iFrameTimesPixel(1) - tPre_this);
        time_text.String = sprintf('Time: %0.2f ms', tthis_frame ./ 1000);

        time_line.Value = tthis_frame;

        img_this = imresize(img_extracted(:, :, k), img_ratio);
        img.CData = img_this;

%         drawnow;
        % plot or update data in this plot
        F(k) = getframe(hf25);
    end
    % make a video clip and save it to the correct location
    close(hf25);
    clear img_extracted

    writerObj = VideoWriter(VidClipFileName);
    writerObj.FrameRate = 0.4 * median(1000./diff(iFrameTimesPixel));
%     writerObj.Quality = 100;
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

