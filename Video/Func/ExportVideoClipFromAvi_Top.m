function ExportVideoClipFromAvi_Top(BehClass, FrameInfo, ClipInfo, scn_scale, remake, x_rev)

% ExportVideoClipFromAvi(thisTable, FrameTable, 'Event', VideoEvent,'ANM', ANM, ...
%     'Pre', Pre, 'Post', Post, 'BehaviorType', BehaviorType, 'Session', Session, 'Remake', 1)
% Jianing Yu

% 5/1/2021
% 4/16/2022

% revised by ZZH, 5/5/2023

if nargin < 5
    remake = 0;
    x_rev  = 0;
elseif nargin < 6
    x_rev  = 0;
end

scale_ratio     =   1 / scn_scale;
color           =   GPSColor();

%% get information
% get behavior data
BehTable = BehClass.BehavTable;
anm      = BehClass.Subject;
beh_type = BehClass.Task;
session  = BehClass.Session;

% get clip parameters
event        = ClipInfo.VideoEvent;
tPre         = ClipInfo.Pre;
tPost        = ClipInfo.Post;
time_elapsed = -tPre:0.1:tPost;

%  in mili-seconds, timing of selected events.
tEventBpod  = FrameInfo.tEventBpod*1000;
Trials      = FrameInfo.Trials;
% in mili-seconds, timing of each frame in Bpod's world
tFramesBpod  = FrameInfo.tFramesInBpod;

% set up video clip storage folder
thisView   = "Top";
viewFolder = fullfile(pwd, thisView);
% viewFolder = ClipInfo.VideoFolderTop;
clipFolder = fullfile(viewFolder, 'Clips');
if ~isfolder(clipFolder)
    mkdir(clipFolder);
end

% % setting up metadata
% VidsMeta = struct('Session', [], 'Event', [], 'EventIndex', [], 'Performance', [], 'EventTime', [], 'FrameTimesB', [], 'VideoOrg', [], 'FrameIndx', [], 'Code', [], 'CreatedOn', []);

%% Start making videos
fprintf("\nStart video clipping ...\n")
video_accum = 0;

wait_bar = waitbar(0, sprintf('1 / %d', length(tEventBpod)), 'Name', sprintf('Clipping_%s_%s', anm, session));
for i = 1:length(tEventBpod) % i is also the trial number

    i_trial = Trials(i);

    if ~isvalid(wait_bar)
        fprintf("\n****** Interrupted ******\n");
        fprintf("%d / %d clips have been generated\n", i-1, length(tEventBpod));
        return
    end
    waitbar(i/length(tEventBpod), wait_bar, sprintf('%d / %d', i, length(tEventBpod)));

    itEvent = tEventBpod(i);

    IndThisClip = find(tFramesBpod>=itEvent-tPre & tFramesBpod<=itEvent+tPost);
    if isempty(IndThisClip)
        continue
    end
    [~, IndThisFrame] = min(abs(tFramesBpod - itEvent));

    % check if a video has been created and check if we want to
    % re-create the same video
    switch event
        case 'PortSamplePokeTime'
            ClipName = sprintf('%s_%s_SamplePokeTrial%03d', anm, session, round(itEvent));

        case 'PortCenterPokeTime'
            ClipName = sprintf('%s_%s_CenterPokeTrial%03d', anm, session, round(itEvent));

        case 'CentInTime'
            ClipName = sprintf('%s_%s_CentIn%03d_%sView', anm, session, round(itEvent), thisView);

        case 'CentOutTime'
            ClipName = sprintf('%s_%s_ChoiceTrial%03d_%sView', anm, session, round(itEvent), thisView);
    end

    VidClipFileName = fullfile(clipFolder, [ClipName '.mp4']);
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
    elseif ~strcmp(FrameInfo.MyVidFiles{(IndThisClip(1))}, FrameInfo.MyVidFiles{(IndThisFrame)}) || ~strcmp(FrameInfo.MyVidFiles{(IndThisClip(end))}, FrameInfo.MyVidFiles{(IndThisFrame)})
        % same video file
        continue
    end
    EventFrame   = IndThisFrame - IndThisClip(1) + 1;
    NumFrames    = length(IndThisClip);
    NumFramePre  = EventFrame - 1;
    NumFramePost = NumFrames - EventFrame;

    tPre_this = tFramesBpod(IndThisFrame) - tFramesBpod(IndThisClip(1));

    % poke events
%     SamplePokeTime  =   BehTable.PortSamplePokeTime(i)*1000-iFrameTimesBpod(1);
    CentInTime      =   (BehTable.TrialStartTime(i) + BehTable.CentInTime(i))*1000     - iFrameTimesBpod(1) - tPre_this;
    CentOutTime     =   (BehTable.TrialStartTime(i) + BehTable.CentOutTime(i))*1000    - iFrameTimesBpod(1) - tPre_this;
    ChoicePokeTime  =   (BehTable.TrialStartTime(i) + BehTable.ChoicePokeTime(i))*1000 - iFrameTimesBpod(1) - tPre_this;

    poke_state      =   .5*ones(1, length(time_elapsed));
    poke_state(time_elapsed>=CentInTime & time_elapsed<CentOutTime) = 0.2;

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

    VidMeta.Subject      = anm;
    VidMeta.Session      = session;
    VidMeta.Event        = event;
    VidMeta.EventIndex   = i_trial;
%     VidMeta.EventTimeE   = itEvent/1000; % Event time in sec (Ephys)
    VidMeta.EventTimeB   = itEvent/1000; % Event time in sec (Bpod)
%     VidMeta.FrameTimesE  = tFramesEphys(IndThisClip); % frame time in ms in behavior time
    VidMeta.FrameTimesB  = tFramesBpod(IndThisClip); % frame time in ms in behavior time
    VidMeta.FrameIndx    = IndThisClip; % frame index in original video
    VidMeta.NumFrames    = NumFrames;
    VidMeta.EventFrame   = EventFrame; % frame index of event onset
    VidMeta.NumFramePre  = NumFramePre;
    VidMeta.NumFramePost = NumFramePost;
    VidMeta.VideoName    = VidClipFileName;
    VidMeta.Code         = mfilename('fullpath');
    VidMeta.CreatedOn    = date; % today's date

    % Extract frames
    VidFrameIndx_thisfile = FrameInfo.AviFrameIndx(IndThisClip);  % these are the frame index in this video
    this_video = fullfile(viewFolder, FrameInfo.MyVidFiles{IndThisFrame});
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
    H_scl = scale_ratio*H;
    W_scl = scale_ratio*W;
% % 
% %     H_beh = ceil(.3*H_scl);
% %     if mod(H_beh, 2)
% %         H_beh = H_beh+1;
% %     end

    VidMeta.SizeVideo = [H W];

    video_accum = video_accum + 1;
    if video_accum==1
        VidsMeta = VidMeta;
    else
        VidsMeta(video_accum) = VidMeta;
    end

    % height = height - 200;

    %% Make videos

    k = 1;
    hf25 = figure(25); clf
    set(hf25, 'name', thisView, 'units', 'pixels', 'position', [5 50 W_scl H_scl], ...
        'PaperPositionMode', 'auto', 'color', 'w', 'renderer', 'opengl', 'toolbar', 'none', 'resize', 'off', 'Visible', 'on');

    ha = axes;
    set(ha, 'units', 'pixels', 'position', [1 1 W_scl H_scl], 'nextplot', 'add', 'xlim', [.5 W+.5], 'ylim', [.5 H+.5], 'ydir', 'reverse')
    axis off

    % plot this frame:
    img = imagesc(ha, img_extracted(:, :, k), [0 250]);
    colormap('gray');

    % plot some behavior data
% %     tthis_frame   = round(iFrameTimesBpod(k) - iFrameTimesBpod(1) - tPre_this);
% %     time_of_frame = sprintf('%3.0f', tthis_frame);
% % 
% %     text(W-20, 40,  sprintf('%s %s', anm, session), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
% %     text(W-20, 90,  beh_type, 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
% %     text(W-20, 140,  sprintf('Trial %03d', i_trial), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
% %     text(W-20, 190,  sprintf('FP: %d ms', BehTable.FP(i)*1000), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
% %     text(W-20, 240,  sprintf('RT: %d ms', round(1000*BehTable.RT(i))), 'color', [255 255 255]/255, 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
% %     text(W-20, 290,  thisOutcome, 'color', color.(thisOutcome), 'FontSize', 20, 'fontweight', 'bold', 'HorizontalAlignment', 'right')
% %     time_text = text(20, 40, [time_of_frame ' ms'], 'color', [255 215 0]/255, 'FontSize', 22,'fontweight', 'bold');
% %     % plot some important behavioral events
% % 
% %     ha2 = axes;
% %     set(ha2, 'units', 'pixels', 'position', [0.05*scale_ratio*W 0.11*scale_ratio*H 0.9*scale_ratio*W 0.18*scale_ratio*H], ...
% %         'nextplot', 'add', 'xtick', [-tPre:500:tPost], 'xlim', [-tPre tPost], ...
% %         'ycolor', 'none', 'ylim', [0 1.25], 'tickdir', 'out', 'FontSize', 20) %#ok<NBRAK>
% %     ha2.XLabel.String = 'Time (ms)';
% %     ha2.XLabel.FontWeight = 'bold';
% %     ha2.XLabel.FontSize = 20;
% % 
% %     time_line = xline(ha2, tthis_frame, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'Alpha', 0.6);
% % 
% %     xline(ha2, thisFP, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
% % 
% %     stairs(ha2, time_elapsed, poke_state, 'Color', 'k', 'LineWidth', 2.5);
% %     text(ha2, -tPre+5, 0.35, "Center poke", 'Color', 'k', 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
% % 
% %     patch(ha2, 'XData', [ChoicePokeTime ChoicePokeTime ChoicePokeTime ChoicePokeTime] + [0 40 40 0], ...
% %         'YData', [.2 .2 .45 .45], ...
% %         'FaceColor', color.(thisOutcome), 'EdgeColor', 'none');
% % 
% %     text(ha2, -tPre+5, 0.8, "Choice cue", 'Color', color.Cue, 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
% %     patch(ha2, 'XData', [ChoiceCueTime flip(ChoiceCueTime)], ...
% %         'YData', [.7 .7 .85 .85], ...
% %         'FaceColor', color.Cue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
% % 
% %     text(ha2, -tPre+5, 1.1, "Trigger cue", 'Color', [30 144 255] / 255, 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
% %     patch(ha2, 'XData', [TriggerCueTime flip(TriggerCueTime)], ...
% %         'YData', [1.0 1.0 1.15 1.15], ...
% %         'FaceColor', [30 144 255] / 255, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

    F(k) = getframe(hf25);
    % plot or update data in this plot
    for k = 2:nframe

% %         tthis_frame = round(iFrameTimesBpod(k) - iFrameTimesBpod(1) - tPre_this);
% %         time_of_frame = sprintf('%3.0f', tthis_frame);
% %         time_text.String = [time_of_frame ' ms'];
% % 
% %         time_line.Value = tthis_frame;

        img.CData = img_extracted(:, :, k);

%         drawnow;
        % plot or update data in this plot
        F(k) = getframe(hf25);
    end
    % make a video clip and save it to the correct location
    close(hf25);
    clear img_extracted

    warning('off', 'MATLAB:audiovideo:VideoWriter:mp4FramePadded');
    writerObj = VideoWriter(VidClipFileName, 'MPEG-4');
    Fs = median(1000./diff(iFrameTimesBpod));
    Fs = roundn(Fs, 1);
    writerObj.FrameRate = 0.4 * Fs;
    writerObj.Quality   = 100;
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

    MetaFileName = fullfile(clipFolder, [ClipName, '.mat']);
    save(MetaFileName, 'VidMeta');

end

close(wait_bar);

