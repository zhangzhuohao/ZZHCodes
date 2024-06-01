% clear;

DLCCrop = [0 1280 0 1024];

BodyParts = {
    'ear_base_left', ...
    'ear_base_right', ...
    'body_1'
    };

Ports = {
    'port_left', ...
    'port_right', ...
    };

%%
opts = delimitedTextImportOptions("NumVariables", 37);

opts.DataLines = [4, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["frame", ...
    "ear_tip_left_x",   "ear_tip_left_y",   "ear_tip_left_lh", ...
    "ear_base_left_x",  "ear_base_left_y",  "ear_base_left_lh", ...
    "ear_tip_right_x",  "ear_tip_right_y",  "ear_tip_right_lh", ...
    "ear_base_right_x", "ear_base_right_y", "ear_base_right_lh", ...
    "body_1_x",         "body_1_y",         "body_1_lh", ...
    "body_2_x",         "body_2_y",         "body_2_lh", ...
    "body_3_x",         "body_3_y",         "body_3_lh", ...
    "tail_left_x",      "tail_left_y",      "tail_left_lh", ...
    "tail_right_x",     "tail_right_y",     "tail_right_lh", ...
    "port_left_x",      "port_left_y",      "port_left_lh", ...
    "port_right_x",     "port_right_y",     "port_right_lh", ...
    "port_center_x",    "port_center_y",    "port_center_lh"];

opts.VariableTypes = repmat("double", 1, 37);

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%%
% ClipFolder  = uigetdir('D:\YuLab\Work\GPS\Video\Leopold\GPS_12_ThreeFPHoldSRT\');
% if ~ClipFolder
%     return
% end
% [ViewFolder, dir_name] = fileparts(ClipFolder);
% [DateFolder, ~] = fileparts(ViewFolder);
% [TaskFolder, ~] = fileparts(DateFolder);
% if ~strcmp(dir_name, 'Clips')
%     fprintf("\nPlease select a 'Clips' folder.\n");
%     return
% end
%%
TaskFolder = uigetdir('D:\YuLab\Work\GPS\Video\');
if ~TaskFolder
    return
end

% TaskFolder = VideoFolderParent;
SessionFolders = get_folders(TaskFolder, "FolderType", 'Session');

[~, SessionsAll] = arrayfun(@(x) fileparts(x), SessionFolders);
Sessions = unique(SessionsAll);

[SessionInd, tf] = listdlg("ListString", Sessions, "ListSize", [200, 200]);
if ~tf
    return
end

SessionFolders = SessionFolders(ismember(SessionsAll, Sessions(SessionInd)));

for s = 1:length(SessionFolders)
% ClipFolder  = 'D:\YuLab\Work\GPS\Video\Kennard\GPS_13_ThreeFPHoldSRT\20240313\Top\Clips';
ClipFolder = fullfile(SessionFolders(s), "Top", "Clips");
if ~isfolder(ClipFolder)
    fprintf("no clip folder in %s", fullfile(SessionFolders(s), "Top"));
end

%%
ClipInfo = split(ClipFolder, '\');
view     = ClipInfo{end-1};
session  = ClipInfo{end-2};
anm      = ClipInfo{end-4};

Drives = string(char('A':'Z')');
for i = 1:length(Drives)
    if isfolder(Drives(i)+":\OneDrive")
        DataFolder  = Drives(i) + ":\OneDrive\YuLab\Work\GPS\Data\";
        break;
    end
end

ANMInfoFile     =   fullfile(DataFolder, "ANMInfo.xlsx");
ANMInfo         =   readtable(ANMInfoFile, "Sheet", anm, "TextType", "string");
ANMInfo.Session =   string(ANMInfo.Session);

SessionInfo = ANMInfo(ANMInfo.Session==session, :);
disp(SessionInfo);

SessionFolderPart = extractAfter(SessionInfo.SessionFolder, "Data\");
SessionFolder = fullfile(DataFolder, SessionFolderPart);

BehTableFile = dir(fullfile(SessionFolder, '*SessionTable*.csv'));
BehTable     = readtable(fullfile(SessionFolder, BehTableFile.name));

%%
OutFile = fullfile(ClipFolder, 'DLCTrackingOutAuto.mat');

%%
ClipFiles = dir(fullfile(ClipFolder, '*.avi'));
NumClips  = length(ClipFiles);

%%
wait_bar = waitbar(0, sprintf('1 / %d', NumClips), 'Name', sprintf('Extracting_%s_%s', anm, session));

fprintf("\n");
for i = 1:NumClips

    if ~isvalid(wait_bar)
        fprintf("\n****** Interrupted ******\n");
        fprintf("%d / %d trials have been extracted\n", i-1, NumClips);
        return
    end
    waitbar(i/NumClips, wait_bar, sprintf('%d / %d', i, NumClips));

    clip_filename = ClipFiles(i).name;
    meta_filename = [clip_filename(1:end-4) '.mat'];
    dlcTable_name = dir(fullfile(ClipFolder, [clip_filename(1:end-4) '*00.csv']));
    dlcTable_name = dlcTable_name.name;

    load(fullfile(ClipFolder, meta_filename), "VidMeta");

    if isfield(VidMeta, 'FrameTimesB')
        t_frames = VidMeta.FrameTimesB - VidMeta.EventTime*1000;  % in ms
    elseif isfield(VidMeta, 'FrameTimesMED')
        t_frames = VidMeta.FrameTimesMED*1000 - VidMeta.EventTime*1000;  % in ms
    else
        error('Check VidMeta. Cannot find FrameTimesB or FrameTimesMED')
    end

    trial_info = BehTable(BehTable.Trials==VidMeta.EventIndex, :);

    FrameBeg = find(t_frames>=-120, 1);
    switch trial_info.Outcome{1}
        case {'Correct', 'Wrong'}
            FrameEnd = find(t_frames<=(trial_info.ChoicePokeTime-trial_info.(VidMeta.Event))*1000, 1, 'last');
        case {'Premature'}
            FrameEnd = find(t_frames<=(trial_info.CentOutTime-trial_info.(VidMeta.Event))*1000+200, 1, 'last');
        case {'Late', 'LateMiss', 'LateCorrect', 'LateWrong'}
            FrameEnd = find(t_frames<=(trial_info.CentOutTime-trial_info.(VidMeta.Event))*1000+200, 1, 'last');
            trial_info.Outcome{1} = 'Late';
    end
    t_frames = t_frames(FrameBeg:FrameEnd);

    % start to update video file
    D = readtable(fullfile(ClipFolder, dlcTable_name), opts);
    % This is the list of body parts tracked

    port_loc.L = [D.port_left_x D.port_left_y];
    port_loc.R = [D.port_right_x D.port_right_y];

    if ~any(D.port_left_lh>0.99) || ~any(D.port_right_lh>0.99)
        fprintf("\nDrop %d for invisible left/right port\n", VidMeta.EventIndex);
        continue
    else
        port_loc.L = mean(port_loc.L(D.port_left_lh>0.99 , :));
        port_loc.R = mean(port_loc.R(D.port_right_lh>0.99, :));
    end

    DropOut = 0;
    if size(D, 1) < 198
        DropOut = 1;
        fprintf("\nDrop %d for frame loss\n", VidMeta.EventIndex);
        continue;
    end

    for j = 1:length(BodyParts)

        body_part =     BodyParts{j};

        nframes   =     D.frame(FrameBeg:FrameEnd, 1) + 1;
        x_pos     =     D.([body_part '_x'])(FrameBeg:FrameEnd) + DLCCrop(1);
        y_pos     =     D.([body_part '_y'])(FrameBeg:FrameEnd) + DLCCrop(3);
        lh        =     D.([body_part '_lh'])(FrameBeg:FrameEnd);

        if j <= 2 % check ear_base_left and ear_base_right
            bad_label = find(lh < 0.8);
            if ~isempty(bad_label)
                DropOut = 1;

                this_clip = VideoReader(fullfile(ClipFolder, clip_filename));

                for k = 1:length(bad_label)

                    this_frame_bad = read(this_clip, nframes(bad_label(k)));

                    this_frame_filename = sprintf('%s_%03d.jpg', clip_filename(1:end-4), nframes(bad_label(k)));
                    
                    bad_label_folder = fullfile(ClipFolder, 'BadLabels');
                    if ~isfolder(bad_label_folder)
                        mkdir(bad_label_folder);
                    end

                    this_frame_filepath = fullfile(bad_label_folder, this_frame_filename);

                    imwrite(this_frame_bad, this_frame_filepath);

                end
            end
        end
    end
    if DropOut
        fprintf("\nDrop %d for bad labelling\n", VidMeta.EventIndex);
        continue;
    end

    for j = 1:length(BodyParts)

        body_part =     BodyParts{j};

        nframes   =     D.frame(FrameBeg:FrameEnd, 1) + 1;
        x_pos     =     D.([body_part '_x'])(FrameBeg:FrameEnd) + DLCCrop(1);
        y_pos     =     D.([body_part '_y'])(FrameBeg:FrameEnd) + DLCCrop(3);
        lh        =     D.([body_part '_lh'])(FrameBeg:FrameEnd);

        if j <=2
            bad_label = find(lh < 0.8);
        end
        if ~isempty(bad_label)
            fh    =     figure(31); clf(fh);
            set(fh, 'name', body_part, 'units', 'centimeter', 'position', [5 2 20 1.3*20*1024/1280]);

            ax = axes;
            set(ax, 'units', 'norm', 'position', [.05 .05 .9 .9], 'YDir', 'Reverse', 'nextplot', 'add', 'xcolor', 'none', 'ycolor', 'none');

            this_clip = VideoReader(fullfile(ClipFolder, clip_filename));
            imagesc(ax, read(this_clip, nframes(bad_label(1))));
            drawnow();

            for k = 1:length(bad_label)
                clf(fh);
                ax = axes;
                set(ax, 'units', 'norm', 'position', [.05 .05 .9 .9], 'YDir', 'Reverse', 'nextplot', 'add', 'xcolor', 'none', 'ycolor', 'none');

                imagesc(ax, read(this_clip, nframes(bad_label(k))));
                scatter(ax, x_pos(bad_label(k)), y_pos(bad_label(k)), 'or');

                title(body_part, 'Interpreter', 'none');
                axis equal

                [xi, yi] = getpts(ax);

                x_pos(bad_label(k)) = xi(end-1);
                y_pos(bad_label(k)) = yi(end-1);

            end

            close(fh);

        end
        pos_data = [x_pos, y_pos, lh, t_frames, nframes];

        if isempty(dir(OutFile))

            DLCTrackingOut = [];
            DLCTrackingOut.Session = VidMeta.Session;
            DLCTrackingOut.Event = VidMeta.Event;
            DLCTrackingOut.BodyParts{1} = body_part;
            ind_body = 1;

            DLCTrackingOut.PoseTracking(1).PosInfo          = {'x', 'y', 'likelihood', 'time(ms)'};
            DLCTrackingOut.PoseTracking(1).PosData{1}       = pos_data;
            DLCTrackingOut.PoseTracking(1).Performance{1}   = trial_info.Outcome{1};

            if isfield(VidMeta, 'BpodTrialNum')
                DLCTrackingOut.PoseTracking(1).BpodEventIndex(:, 1) = [VidMeta.BpodTrialNum; VidMeta.EventTime];
                % Trials tracked
                app.TrialTrackedLabel.Text = sprintf('tracked: %2.0d', size(DLCTrackingOut.PoseTracking(1).BpodEventIndex, 2));
            elseif isfield(VidMeta, 'PressTrialNum')
                DLCTrackingOut.PoseTracking(1).MEDPressIndex(:, 1) = [VidMeta.PressTrialNum; VidMeta.EventTime];
                % Trials tracked
                app.TrialTrackedLabel.Text = sprintf('tracked: %2.0d', size(DLCTrackingOut.PoseTracking(1).MEDPressIndex, 2));
            elseif isfield(VidMeta, 'EventIndex')
                DLCTrackingOut.PoseTracking(1).BpodEventIndex(:, 1) = [VidMeta.EventIndex; VidMeta.EventTime];
                % Trials tracked
                app.TrialTrackedLabel.Text = sprintf('tracked: %2.0d', size(DLCTrackingOut.PoseTracking(1).BpodEventIndex, 2));
            else
                error('Cannot find BpodTrialNum or PressTrialNum or EventIndex. Check index')
            end

            DLCTrackingOut.PortLoc(1) = port_loc;

            save(OutFile, 'DLCTrackingOut');
        
        else

            load(OutFile, "DLCTrackingOut");
 
            % load DLCTracingOut
            % find index of body parts
            ind_body = find(strcmp(DLCTrackingOut.BodyParts, body_part));

            if isempty(ind_body) % no such body parts
                DLCTrackingOut.BodyParts = [DLCTrackingOut.BodyParts body_part];
                ind_body = length(DLCTrackingOut.BodyParts);
                DLCTrackingOut.PoseTracking(ind_body) = struct('PosInfo', [], 'PosData', [], 'BpodEventIndex', [], 'Performance', []);
                DLCTrackingOut.PoseTracking(ind_body).PosInfo = {'x', 'y', 'likelihood', 'time(ms)' };
            end

            % new index
            if isfield(VidMeta, 'BpodTrialNum')
                if isempty(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex) || isempty(find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.BpodTrialNum, 1))
                    DLCTrackingOut.PoseTracking(ind_body).PosData{end+1}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, end+1)  = [VidMeta.BpodTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{end+1}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(end+1) = port_loc;
                    end
                else
                    ind_old = find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.BpodTrialNum);
                    DLCTrackingOut.PoseTracking(ind_body).PosData{ind_old}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, ind_old)  = [VidMeta.BpodTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{ind_old}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(ind_old) = port_loc;
                    end
                end

            elseif isfield(VidMeta, 'EventIndex')
                if isempty(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex) || isempty(find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.EventIndex, 1))
                    DLCTrackingOut.PoseTracking(ind_body).PosData{end+1}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, end+1)  = [VidMeta.EventIndex; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{end+1}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(end+1) = port_loc;
                    end                
                else
                    ind_old = find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.EventIndex);
                    DLCTrackingOut.PoseTracking(ind_body).PosData{ind_old}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, ind_old)  = [VidMeta.EventIndex; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{ind_old}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(ind_old) = port_loc;
                    end
                end

            elseif isfield(VidMeta, 'PressTrialNum')
                if isempty(DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex) || isempty(find(DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(1, :)==VidMeta.PressTrialNum, 1))
                    DLCTrackingOut.PoseTracking(ind_body).PosData{end+1}        = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(:, end+1)  = [VidMeta.PressTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{end+1}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(end+1) = port_loc;
                    end                
                else
                    ind_old = find(DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(1, :)==VidMeta.PressTrialNum);
                    DLCTrackingOut.PoseTracking(ind_body).PosData{ind_old}        = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(:, ind_old)  = [VidMeta.PressTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{ind_old}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(ind_old) = port_loc;
                    end
                end

            else
                error('Check VidMeta')
            end

            save(OutFile, 'DLCTrackingOut');
        end
    end

end

close(wait_bar);
end
