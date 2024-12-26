clear;

DLCCrop = [0 1280 300 900];

BodyParts = {
    'snout', ...
    'ear_left', ...
    'ear_right', ...
    'body_1'
    };

Corridor = {
    'corridor_l_u', ...
    'corridor_l_d', ...
    'corridor_r_u', ...
    'corridor_r_d', ...
    };

%%
opts = delimitedTextImportOptions("NumVariables", 34);

opts.DataLines = [4, Inf];
opts.Delimiter = ",";

opts.VariableNames = [
    "frame", ...
    "corridor_l_u_x", "corridor_l_u_y", "corridor_l_u_lh", ...
    "corridor_l_d_x", "corridor_l_d_y", "corridor_l_d_lh", ...
    "corridor_r_u_x", "corridor_r_u_y", "corridor_r_u_lh", ...
    "corridor_r_d_x", "corridor_r_d_y", "corridor_r_d_lh", ...
    "snout_x",        "snout_y",        "snout_lh", ...
    "ear_left_x",     "ear_left_y",     "ear_left_lh", ...
    "ear_right_x",    "ear_right_y",    "ear_right_lh", ...
    "body_1_x",       "body_1_y",       "body_1_lh", ...
    "body_2_x",       "body_2_y",       "body_2_lh", ...
    "body_3_x",       "body_3_y",       "body_3_lh", ...
    "tail_x",         "tail_y",         "tail_lh"
    ];

opts.VariableTypes = repmat("double", 1, 34); 

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%%
TaskFolder = uigetdir('T:\YuLab\Work\GPS\Video\');
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

%%
for s = 1:length(SessionFolders)
% ClipFolder  = 'D:\YuLab\Work\GPS\Video\Kennard\GPS_13_ThreeFPHoldSRT\20240313\Field\Clips';
ClipFolder = fullfile(SessionFolders(s), "Field", "Clips");
if ~isfolder(ClipFolder)
    fprintf("no clip folder in %s", fullfile(SessionFolders(s), "Field"));
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
    st = trial_info.ST * 1000;

    D = readtable(fullfile(ClipFolder, dlcTable_name), opts);
    lh_l = D.ear_left_lh;
    lh_r = D.ear_right_lh;

    valid_id = find(lh_l>.8 & lh_r>.8);
    if isempty(valid_id)
        fprintf("\nDrop %d for no good labels\n", VidMeta.EventIndex);
        continue
    end
    df = diff(valid_id);
    f_seg = find(df>2);
    if ~isempty(f_seg)
        f_seg = [0; f_seg; length(valid_id)];
        valid_seg = cell(1, length(f_seg)-1);
        for j = 1:length(f_seg)-1
            valid_seg{j} = valid_id(f_seg(j)+1:f_seg(j+1));
        end
        len_seg = cellfun(@length, valid_seg);
        [~, id] = max(len_seg);
        valid_id = valid_seg{id};
    end

    FrameBeg = valid_id(1);
    FrameEnd = valid_id(end);

    t_frames = t_frames(FrameBeg:FrameEnd);

    % start to update video file
    
    % This is the list of body parts tracked
    corridor_loc.L_U = [D.corridor_l_u_x D.corridor_l_u_y] + DLCCrop([1 3]);
    corridor_loc.L_D = [D.corridor_l_d_x D.corridor_l_d_y] + DLCCrop([1 3]);
    corridor_loc.R_U = [D.corridor_r_u_x D.corridor_r_u_y] + DLCCrop([1 3]);
    corridor_loc.R_D = [D.corridor_r_d_x D.corridor_r_d_y] + DLCCrop([1 3]);

    if ~any(D.corridor_l_u_lh>0.95) || ~any(D.corridor_l_d_lh>0.95) || ~any(D.corridor_r_u_lh>0.95) || ~any(D.corridor_r_d_lh>0.95)
        fprintf("\nDrop %d for invisible corridor corner\n", VidMeta.EventIndex);
        continue
    else
        corridor_loc.L_U = mean(corridor_loc.L_U(D.corridor_l_u_lh>0.95, :), 1);
        corridor_loc.L_D = mean(corridor_loc.L_D(D.corridor_l_d_lh>0.95, :), 1);
        corridor_loc.R_U = mean(corridor_loc.R_U(D.corridor_r_u_lh>0.95, :), 1);
        corridor_loc.R_D = mean(corridor_loc.R_D(D.corridor_r_d_lh>0.95, :), 1);
    end

    DropOut = 0;

    dt = diff(t_frames);
    if any(dt>=3*median(dt))
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

        if ismember(body_part, {'ear_left', 'ear_right'}) % check ear_base_left and ear_base_right
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
% % % % % % % % % % % % % %         continue;
    end

    for j = 1:length(BodyParts)

        body_part =     BodyParts{j};

        nframes   =     D.frame(FrameBeg:FrameEnd, 1) + 1;
        x_pos     =     D.([body_part '_x'])(FrameBeg:FrameEnd) + DLCCrop(1);
        y_pos     =     D.([body_part '_y'])(FrameBeg:FrameEnd) + DLCCrop(3);
        lh        =     D.([body_part '_lh'])(FrameBeg:FrameEnd);

        if ismember(body_part, {'ear_left', 'ear_right'}) % check ear_base_left and ear_base_right
            bad_label = find(lh < 0.8);
            x_pos(bad_label) = [];
            y_pos(bad_label) = [];
            t_frames_nan = t_frames;
            t_frames_nan(bad_label) = [];
            x_pos = interp1(t_frames_nan, x_pos, t_frames, 'linear');
            y_pos = interp1(t_frames_nan, y_pos, t_frames, 'linear');
        end
%         if ~isempty(bad_label)
%             fh    =     figure(31); clf(fh);
%             set(fh, 'name', body_part, 'units', 'centimeter', 'position', [5 2 20 1.3*20*1024/1280]);
% 
%             ax = axes;
%             set(ax, 'units', 'norm', 'position', [.05 .05 .9 .9], 'YDir', 'Reverse', 'nextplot', 'add', 'xcolor', 'none', 'ycolor', 'none');
% 
%             this_clip = VideoReader(fullfile(ClipFolder, clip_filename));
%             imagesc(ax, read(this_clip, nframes(bad_label(1))));
%             drawnow();
% 
%             for k = 1:length(bad_label)
%                 clf(fh);
%                 ax = axes;
%                 set(ax, 'units', 'norm', 'position', [.05 .05 .9 .9], 'YDir', 'Reverse', 'nextplot', 'add', 'xcolor', 'none', 'ycolor', 'none');
% 
%                 imagesc(ax, read(this_clip, nframes(bad_label(k))));
%                 scatter(ax, x_pos(bad_label(k)), y_pos(bad_label(k)), 'or');
% 
%                 title(body_part, 'Interpreter', 'none');
%                 axis equal
% 
%                 [xi, yi] = getpts(ax);
% 
%                 x_pos(bad_label(k)) = xi(end-1);
%                 y_pos(bad_label(k)) = yi(end-1);
% 
%             end
%             close(fh);
%         end

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

            DLCTrackingOut.PortLoc(1) = corridor_loc;

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
                        DLCTrackingOut.PortLoc(end+1) = corridor_loc;
                    end
                else
                    ind_old = find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.BpodTrialNum);
                    DLCTrackingOut.PoseTracking(ind_body).PosData{ind_old}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, ind_old)  = [VidMeta.BpodTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{ind_old}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(ind_old) = corridor_loc;
                    end
                end

            elseif isfield(VidMeta, 'EventIndex')
                if isempty(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex) || isempty(find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.EventIndex, 1))
                    DLCTrackingOut.PoseTracking(ind_body).PosData{end+1}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, end+1)  = [VidMeta.EventIndex; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{end+1}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(end+1) = corridor_loc;
                    end                
                else
                    ind_old = find(DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(1, :)==VidMeta.EventIndex);
                    DLCTrackingOut.PoseTracking(ind_body).PosData{ind_old}            = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).BpodEventIndex(:, ind_old)  = [VidMeta.EventIndex; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{ind_old}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(ind_old) = corridor_loc;
                    end
                end

            elseif isfield(VidMeta, 'PressTrialNum')
                if isempty(DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex) || isempty(find(DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(1, :)==VidMeta.PressTrialNum, 1))
                    DLCTrackingOut.PoseTracking(ind_body).PosData{end+1}        = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(:, end+1)  = [VidMeta.PressTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{end+1}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(end+1) = corridor_loc;
                    end                
                else
                    ind_old = find(DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(1, :)==VidMeta.PressTrialNum);
                    DLCTrackingOut.PoseTracking(ind_body).PosData{ind_old}        = pos_data;
                    DLCTrackingOut.PoseTracking(ind_body).MEDPressIndex(:, ind_old)  = [VidMeta.PressTrialNum; VidMeta.EventTime];
                    DLCTrackingOut.PoseTracking(ind_body).Performance{ind_old}        = trial_info.Outcome{1};
                    if j == 1
                        DLCTrackingOut.PortLoc(ind_old) = corridor_loc;
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
