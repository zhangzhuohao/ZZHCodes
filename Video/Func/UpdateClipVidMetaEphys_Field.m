function UpdateClipVidMetaEphys_Field(r, FrameInfo, ClipInfo, scn_scale, remake, x_rev)

% ExportVideoClipFromAvi(thisTable, FrameInfo, 'Event', VideoEvent,'ANM', ANM, ...
%     'Pre', Pre, 'Post', Post, 'BehaviorType', BehaviorType, 'Session', Session, 'Remake', 1)
% Jianing Yu

% 5/1/2021
% 4/16/2022

% revised by ZZH, 5/5/2023
% 
if nargin < 5
    remake = 0;
    x_rev  = 0;
elseif nargin < 6
    x_rev  = 0;
end

scale_ratio     =   1 / scn_scale;
color           =   GPSColor();

%% get information
% load Ephys data
rb = r.Behavior;

% get behavior data
BehClass = r.BehaviorClass;
BehTable = BehClass.BehavTable;
anm      = BehClass.Subject;
beh_type = BehClass.Task;
session  = BehClass.Session;

% get clip parameters
event        = ClipInfo.VideoEvent;
tPreMax      = ClipInfo.Pre;
tPost        = ClipInfo.Post;

% in mili-seconds, timing of selected events.
tEventEphys = FrameInfo.tEventEphys;
tEventBpod  = FrameInfo.tEventBpod;
Trials      = FrameInfo.Trials;
% in mili-seconds, timing of each frame in Bpod's world
tFramesBpod  = FrameInfo.tFramesInBpod;
% in mili-seconds, timing of each frame in Ephys' world
tFramesEphys = FrameInfo.tFramesInEphys;

% set up video clip storage folder
thisView   = "Field";
viewFolder = ClipInfo.VideoFolderField;
clipFolder = fullfile(ClipInfo.VideoFolderField, 'Clips');
if ~isfolder(clipFolder)
    mkdir(clipFolder);
end

% Get Ephys event timings
id_centin = find(strcmp(rb.Labels, 'PokeCentIn'));
t_centin  = rb.EventTimings(rb.EventMarkers==id_centin);
id_centout = find(strcmp(rb.Labels, 'PokeCentOut'));
t_centout  = rb.EventTimings(rb.EventMarkers==id_centout);
id_choicein = find(strcmp(rb.Labels, 'PokeChoiceIn'));
t_choicein  = rb.EventTimings(rb.EventMarkers==id_choicein);
id_trigger = find(strcmp(rb.Labels, 'Trigger'));
t_trigger  = rb.EventTimings(rb.EventMarkers==id_trigger);
id_initout = find(strcmp(rb.Labels, 'PokeInitOut'));
t_initout  = rb.EventTimings(rb.EventMarkers==id_initout);

% setting up metadata
VidsMeta = struct('Subject', [], 'Session', [], 'Event', [], 'EventIndex', [], 'Performance', [], 'EventTime', [], 'FrameTimesE', [], 'FrameTimesB', [], 'VideoOrg', [], 'FrameIndx', [], 'Code', [], 'CreatedOn', []);

%% Start making videos
fprintf("\nStart video clipping ...\n")
video_accum = 0;

wait_bar = waitbar(0, sprintf('1 / %d', length(tEventEphys)), 'Name', sprintf('Clipping_%s_%s', anm, session));
for i = 1:length(tEventEphys) % i is also the trial number
    
    i_trial = Trials(i);
    ind_bpod = find(BehClass.Trials==i_trial);
    ind_ephys = find(rb.TrialID==i_trial);

    if ~isvalid(wait_bar)
        fprintf("\n****** Interrupted ******\n");
        fprintf("%d / %d clips have been generated\n", i-1, length(tEventEphys));
        return
    end
    waitbar(i/length(tEventEphys), wait_bar, sprintf('%d / %d', i, length(tEventEphys)));

    itEvent = tEventEphys(i);
    iShuttleTime = t_centin(ind_ephys) - t_initout(ind_ephys); % Shuttle time of this trial (in ms)
    if (400 + iShuttleTime) < tPreMax
        tPre = iShuttleTime + 400;
    else
        tPre = tPreMax;
    end
    time_elapsed = -tPre:0.1:tPost;

    IndThisClip = find(tFramesEphys>=itEvent-tPre & tFramesEphys<=itEvent+tPost);
    if isempty(IndThisClip)
        continue
    end
    [~, IndThisFrame] = min(abs(tFramesEphys - itEvent));

    % check if a video has been created and check if we want to
    % re-create the same video
    ClipName = sprintf('%s_%s_Trial%03d_FieldView', anm, session, i);

    VidClipFileName = fullfile(clipFolder, [ClipName '.avi']);
    check_this_file = dir(VidClipFileName);

    if ~isempty(check_this_file) && ~remake % found a video clip with the same name, and we don't want to remake the video clip
        continue % move on
    end

    iFrameTimesEphys = tFramesEphys(IndThisClip);
    % make sure the videoclip can be constructed from a single video file
    if itEvent-iFrameTimesEphys(1) < tPre-50
        continue
    elseif iFrameTimesEphys(end)-itEvent < tPost-50
        continue
    elseif ~strcmp(FrameInfo.MyVidFiles{(IndThisClip(1))}, FrameInfo.MyVidFiles{(IndThisClip(end))})
        % same video file
        continue
    end

    EventFrame   = IndThisFrame - IndThisClip(1) + 1;
    NumFrames    = length(IndThisClip);
    NumFramePre  = EventFrame - 1;
    NumFramePost = NumFrames - EventFrame;

    tPre_this = tFramesEphys(IndThisFrame) - tFramesEphys(IndThisClip(1));
    
    % poke events
    InitOutTime = t_initout(ind_ephys) - tFramesEphys(IndThisFrame);
    CentInTime  = t_centin(ind_ephys) - tFramesEphys(IndThisFrame);
    CentOutTime = t_centout(ind_ephys) - tFramesEphys(IndThisFrame);
    poke_state  = .5*ones(1, length(time_elapsed));
    poke_state(time_elapsed>=CentInTime & time_elapsed<CentOutTime) = 0.2;

    ChoicePokeTime = t_choicein - tFramesEphys(IndThisFrame);
    ChoicePokeTime = ChoicePokeTime(ChoicePokeTime<=tPost & ChoicePokeTime>=-tPre);
    if isempty(ChoicePokeTime)
        ChoicePokeTime = nan;
    end

    thisFP = CentInTime + rb.Foreperiods(ind_ephys)*1000;
    switch rb.Outcome(ind_ephys)
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
    ChoiceCueTime  = (BehTable.ChoiceCueTime(ind_bpod,:) - BehTable.CentInTime(ind_bpod))*1000 + t_centin(ind_ephys) - tFramesEphys(IndThisFrame);

    TriggerCueTime = t_trigger - tFramesEphys(IndThisFrame);
    TriggerCueTime = TriggerCueTime(TriggerCueTime<=tPost & TriggerCueTime>=-tPre);
    if isempty(TriggerCueTime)
        TriggerCueTime = [nan nan];
    else
        TriggerCueTime = [0 250] + TriggerCueTime;
    end

    if strcmp(beh_type, "KornblumSRT")
        if rb.CueIndex(ind_ephys)==0
            TriggerCueTime = nan(1,2);
        end
    end

    VidMeta.Subject      = anm;
    VidMeta.Session      = session;
    VidMeta.Event        = event;
    VidMeta.EventIndex   = i_trial;
    VidMeta.EventTimeE   = itEvent/1000; % Event time in sec (Ephys)
    VidMeta.EventTimeB   = tEventBpod(i); % Event time in sec (Bpod)
    VidMeta.FrameTimesE  = tFramesEphys(IndThisClip); % frame time in ms in behavior time
    VidMeta.FrameTimesB  = tFramesBpod(IndThisClip); % frame time in ms in behavior time
    VidMeta.FrameIndx    = IndThisClip; % frame index in original video
    VidMeta.NumFrames    = NumFrames;
    VidMeta.EventFrame   = EventFrame; % frame index of event onset
    VidMeta.NumFramePre  = NumFramePre;
    VidMeta.NumFramePost = NumFramePost;
    VidMeta.VideoName    = VidClipFileName;
    VidMeta.Code         = mfilename('fullpath');
    VidMeta.CreatedOn    = date; % today's date

    video_accum = video_accum + 1;
    if video_accum==1
        VidsMeta = VidMeta;
    else
        VidsMeta(video_accum) = VidMeta;
    end

    MetaFileName = fullfile(clipFolder, [ClipName, '.mat']);
    save(MetaFileName, 'VidMeta');

end

close(wait_bar);

