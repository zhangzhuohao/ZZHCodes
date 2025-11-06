function ExportVideoClipFromAvi(BehTable, IntTable, FrameTable, SessionInfo, ClipInfo, scn_scale, remake, view, x_rev)

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
tBehEvent = 1000*(BehTable.(event) + BehTable.TrialStartTime);
% in mili-seconds, timing of each frame in Bpod's world
tFramesBpod = FrameTable.tFrames2Bpodms;

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
VidsMeta = struct('Session', [], 'Event', [], 'EventIndex', [], 'Performance', [], 'EventTime', [], 'FrameTimesB', [], 'VideoOrg', [], 'FrameIndx', [], 'Code', [], 'CreatedOn', []);

%% Start making videos
fprintf("\nStart video clipping ...\n")
video_accum = 0;

wait_bar = waitbar(0, sprintf('1 / %d', length(tBehEvent)), 'Name', sprintf('Clipping_%s_%s', anm, session));
for i = 1:length(tBehEvent) % i is also the trial number

    i_trial = BehTable.Trials(i);

    if ~isvalid(wait_bar)
        fprintf("\n****** Interrupted ******\n");
        fprintf("%d / %d clips have been generated\n", i-1, length(tBehEvent));
        return
    end
    waitbar(i/length(tBehEvent), wait_bar, sprintf('%d / %d', i, length(tBehEvent)));

    itEvent = tBehEvent(i);

    IndThisClip         =   find(tFramesBpod>=itEvent-tPre & tFramesBpod<=itEvent+tPost);
    if isempty(IndThisClip)
        continue
    end
    [~, IndThisFrame]   =   min(abs(tFramesBpod - itEvent));

    % check if a video has been created and check if we want to
    % re-create the same video
    switch event
        case 'PortSamplePokeTime'
            ClipName = sprintf('%s_%s_SamplePokeTrial%03d', anm, session, i_trial);

        case 'PortCenterPokeTime'
            ClipName = sprintf('%s_%s_CenterPokeTrial%03d', anm, session, i_trial);

        case 'CentInTime'
            ClipName = sprintf('%s_%s_HoldTrial%03d_%sView', anm, session, i_trial, view);

        case 'CentOutTime'
            ClipName = sprintf('%s_%s_ChoiceTrial%03d_%sView', anm, session, i_trial, view);
    end

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
    tPre_this = tFramesBpod(IndThisFrame) - tFramesBpod(IndThisClip(1));

    % poke events
%     SamplePokeTime  =   BehTable.PortSamplePokeTime(i)*1000-iFrameTimesBpod(1);
    CentInTime      =   (BehTable.TrialStartTime(i) + BehTable.CentInTime(i))*1000     - iFrameTimesBpod(1) - tPre_this;
    CentOutTime     =   (BehTable.TrialStartTime(i) + BehTable.CentOutTime(i))*1000    - iFrameTimesBpod(1) - tPre_this;
    ChoicePokeTime  =   (BehTable.TrialStartTime(i) + BehTable.ChoicePokeTime(i))*1000 - iFrameTimesBpod(1) - tPre_this;

    if ~isempty(IntTable)
        IntOnTime       =   CentInTime + 1000*IntTable.On(IntTable.Trials==i_trial);
        IntOffTime      =   CentInTime + IntOnTime + 1000*IntTable.Dur(IntTable.Trials==i_trial);
    end

    poke_state      =   .5*ones(1, length(time_elapsed));
    poke_state(time_elapsed>=CentInTime & time_elapsed<CentOutTime) = 0.2;
    if ~isempty(IntTable)
        for int = 1:length(IntOnTime)
            poke_state(time_elapsed>=IntOnTime(int) & time_elapsed<IntOffTime(int)) = .35;
        end
    end

    thisFP          =   CentInTime + BehTable.FP(i)*1000;
    switch BehTable.Outcome{i}
        case {'Premature', 'Pre'}
            thisOutcome = "Premature";
        case {'Correct', 'Cor'}
            thisOutcome = "Correct";
        case {'Late', 'LateCorrect', 'LateWrong', 'LateMiss'}
            thisOutcome = "Late";
        case {'Wrong', 'Wro'}
            thisOutcome = "Wrong";
        case {'Probe'}
            thisOutcome = "Probe";
    end

    % cue events
    try
        ChoiceCueTime = (BehTable.TrialStartTime(i) + BehTable.ChoiceCueTime(i,:))*1000 - iFrameTimesBpod(1) - tPre_this;
    catch
        ChoiceCueTime = (BehTable.TrialStartTime(i) + [BehTable.ChoiceCueTime_1(i) BehTable.ChoiceCueTime_2(i)])*1000 - iFrameTimesBpod(1) - tPre_this;
    end
    TriggerCueTime  =   [0 250] + (BehTable.TrialStartTime(i) + BehTable.TriggerCueTime(i))*1000 - iFrameTimesBpod(1) - tPre_this;

    if strcmp(beh_type, "KornblumSRT")
        if BehTable.Cued(i)==0
            TriggerCueTime = nan(1,2);
        end
    end

    % build video clips, frame by frame
    F = struct('cdata', [], 'colormap', []);

    VidMeta.Session     =   session;
    VidMeta.Event       =   event;
    VidMeta.EventIndex  =   i_trial;
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

    % height = height - 200;

    %% Make videos

    k = 1;
    hf25 = figure(25); clf
    set(hf25, 'name', thisView, 'units', 'pixels', 'position', [5 50 scale_ratio*W 1.3*scale_ratio*H], ...
        'PaperPositionMode', 'auto', 'color', 'w', 'renderer', 'opengl', 'toolbar', 'none', 'resize', 'off', 'Visible', 'on');

    ha = axes;
    set(ha, 'units', 'pixels', 'position', [0 .3*scale_ratio*H + 1 scale_ratio*W scale_ratio*H], 'nextplot', 'add', 'xlim', [.5 W+.5], 'ylim', [.5 H+.5], 'ydir', 'reverse')
    axis off

    % plot this frame:
    img = imagesc(ha, img_extracted(:, :, k), [0 250]);

    colormap('gray');

    % plot some behavior data
    tthis_frame   = round(iFrameTimesBpod(k) - iFrameTimesBpod(1) - tPre_this);
    time_of_frame = sprintf('%3.0f', tthis_frame);

    text(W-20, 40,  sprintf('%s %s', anm, session), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 90,  beh_type, 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 140,  sprintf('Trial %03d', i_trial), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 190,  sprintf('FP: %d ms', BehTable.FP(i)*1000), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 240,  sprintf('RT: %d ms', round(1000*BehTable.RT(i))), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    text(W-20, 290,  thisOutcome, 'color', color.(thisOutcome), 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
    time_text = text(20, 40, [time_of_frame ' ms'], 'color', [255 215 0]/255, 'FontSize', 22,'fontweight', 'bold');
    % plot some important behavioral events

    ha2 = axes;
    set(ha2, 'units', 'pixels', 'position', [0.05*scale_ratio*W 0.11*scale_ratio*H 0.9*scale_ratio*W 0.18*scale_ratio*H], ...
        'nextplot', 'add', 'xtick', [-tPre:500:tPost], 'xlim', [-tPre tPost], ...
        'ycolor', 'none', 'ylim', [0 1.25], 'tickdir', 'out', 'FontSize', 20) %#ok<NBRAK>
    ha2.XLabel.String = 'Time (ms)';
    ha2.XLabel.FontWeight = 'bold';
    ha2.XLabel.FontSize = 20;

    time_line = xline(ha2, tthis_frame, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'Alpha', 0.6);

    xline(ha2, thisFP, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);

    stairs(ha2, time_elapsed, poke_state, 'Color', 'k', 'LineWidth', 2.5);
    text(ha2, -tPre+5, 0.35, "Center poke", 'Color', 'k', 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');

    patch(ha2, 'XData', [ChoicePokeTime ChoicePokeTime ChoicePokeTime ChoicePokeTime] + [0 40 40 0], ...
        'YData', [.2 .2 .45 .45], ...
        'FaceColor', color.(thisOutcome), 'EdgeColor', 'none');

    text(ha2, -tPre+5, 0.8, "Choice cue", 'Color', color.Cue, 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    patch(ha2, 'XData', [ChoiceCueTime flip(ChoiceCueTime)], ...
        'YData', [.7 .7 .85 .85], ...
        'FaceColor', color.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

    text(ha2, -tPre+5, 1.1, "Trigger cue", 'Color', [30 144 255] / 255, 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    patch(ha2, 'XData', [TriggerCueTime flip(TriggerCueTime)], ...
        'YData', [1.0 1.0 1.15 1.15], ...
        'FaceColor', [30 144 255] / 255, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

    F(k) = getframe(hf25);
    % plot or update data in this plot
    for k = 2:nframe

        tthis_frame = round(iFrameTimesBpod(k) - iFrameTimesBpod(1) - tPre_this);
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

