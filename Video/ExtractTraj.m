clear;

DLCCrop = [120 1180 0 1000];

BodyParts = {
    'ear_base_left', ...
    'ear_base_right', ...
    'pattern_left', ...
    'pattern_right', ...
    };

Ports = {
    'port_left', ...
    'port_right', ...
    };

%%
opts = delimitedTextImportOptions("NumVariables", 28);

opts.DataLines = [4, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["frame", "ear_tip_left_x", "ear_tip_left_y", "ear_tip_left_lh", "ear_base_left_x", "ear_base_left_y", "ear_base_left_lh", "ear_tip_right_x", "ear_tip_right_y", "ear_tip_right_lh", "ear_base_right_x", "ear_base_right_y", "ear_base_right_lh", "pattern_left_x", "pattern_left_y", "pattern_left_lh", "pattern_right_x", "pattern_right_y", "pattern_right_lh", "port_left_x", "port_left_y", "port_left_lh", "port_right_x", "port_right_y", "port_right_lh", "port_center_x", "port_center_y", "port_center_lh"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%%
ClipFolder  = uigetdir('E:\YuLab\Work\GPS\Video\');
[ViewFolder, dir_name] = fileparts(ClipFolder);
if ~strcmp(dir_name, 'Clips')
    fprintf("\nPlease select a 'Clips' folder.\n");
    return
end

%%
ClipInfo = split(ClipFolder, '\');
view     = ClipInfo{end-1};
session  = ClipInfo{end-2};
anm      = ClipInfo{end-3};

ANMInfoFile = 'E:\YuLab\Work\GPS\Data\ANMInfo.xlsx';
ANMInfo     = readtable(ANMInfoFile, 'Sheet', anm);
SessionInfo = ANMInfo(ANMInfo.Session==str2double(session), :);

BehTableFile = dir(fullfile(SessionInfo.SessionFolder{1}, '*SessionTable*.csv'));
BehTable     = readtable(fullfile(SessionInfo.SessionFolder{1}, BehTableFile.name));

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
        continue
    else
        port_loc.L = mean(port_loc.L(D.port_left_lh>0.99 , :));
        port_loc.R = mean(port_loc.R(D.port_right_lh>0.99, :));
    end

    DropOut = 0;
    for j = 1:length(BodyParts)

        body_part =     BodyParts{j};

        nframes   =     D.frame(FrameBeg:FrameEnd, 1) + 1;
        x_pos     =     D.([body_part '_x'])(FrameBeg:FrameEnd) + DLCCrop(1);
        y_pos     =     D.([body_part '_y'])(FrameBeg:FrameEnd) + DLCCrop(3);
        lh        =     D.([body_part '_lh'])(FrameBeg:FrameEnd);

        bad_label =     find(lh < 0.8);
        if ~isempty(bad_label)
            DropOut = 1;
        end
    end
    if DropOut
        continue;
    end

    for j = 1:length(BodyParts)

        body_part =     BodyParts{j};

        nframes   =     D.frame(FrameBeg:FrameEnd, 1) + 1;
        x_pos     =     D.([body_part '_x'])(FrameBeg:FrameEnd) + DLCCrop(1);
        y_pos     =     D.([body_part '_y'])(FrameBeg:FrameEnd) + DLCCrop(3);
        lh        =     D.([body_part '_lh'])(FrameBeg:FrameEnd);

        bad_label =     find(lh < 0.8);
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

