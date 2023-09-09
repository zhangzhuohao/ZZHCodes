%%
clear; close all

ScnScale   = 1.5;
cRange     = [-60 60];
lineLength = 60;
fade_ratio = 0.85;

labelL = 'ear_base_left';
labelR = 'ear_base_right';

mycolormap = customcolormap_preset("red-white-blue");

%%
ClipFolder = uigetdir('E:\YuLab\Work\GPS\Video\');
[ViewFolder, dir_name] = fileparts(ClipFolder);
if ~strcmp(dir_name, 'Clips')
    fprintf("\nPlease select a 'Clips' folder.\n");
    return
end

ClassFile = dir(fullfile(ClipFolder, 'GPSTrajectoryClass*.mat'));
load(fullfile(ClipFolder, ClassFile.name))

indL = find(strcmp(obj.DLCTracking.BodyParts, labelL));
indR = find(strcmp(obj.DLCTracking.BodyParts, labelR));

%%
clipFiles = dir(fullfile(ClipFolder, '*.avi'));
NumClips  = length(clipFiles);
TrialID   = arrayfun(@(x) str2double(x.name(end-6:end-4)), clipFiles);

LabeledVidFolder = fullfile(ClipFolder, "LabeledVideo");
if ~isfolder(LabeledVidFolder)
    mkdir(LabeledVidFolder);
end

for m = 1:NumClips
    %%

    VidLabelFileName = fullfile(LabeledVidFolder, strcat(clipFiles(m).name(1:end-4), '_Labeled.avi'));
    if ~isempty(dir(VidLabelFileName))
        continue;
    end
    
    trial_id = TrialID(m);
    this_id = find(obj.DLCTracking.PoseTracking(indL).BpodEventIndex(1,:)==trial_id);
    if isempty(this_id)
        continue
    end

    traceL = obj.DLCTracking.PoseTracking(indL).PosData{this_id};
    traceR = obj.DLCTracking.PoseTracking(indR).PosData{this_id};
    NumFrames = height(traceL);

    vid = VideoReader(fullfile(ClipFolder, clipFiles(m).name));

    %%
    fig = figure(36); clf(36);
    set(fig, 'name', "TrackedVideo", 'units', 'pixels', 'position', [5 50 vid.Width/ScnScale vid.Height/ScnScale], ...
        'PaperPositionMode', 'auto', 'color', 'w', 'renderer', 'opengl', 'toolbar', 'none', 'resize', 'off');

    ha1 = axes;
    set(ha1, 'units', 'pixels', 'position', [0 0 vid.Width/ScnScale vid.Height/ScnScale], 'nextplot', 'add', ...
        'xlim', [0 vid.Width], 'ylim', [0 vid.Height], 'ydir', 'reverse', 'XColor', 'none', 'YColor', 'none', 'Color', 'none');

    ha2 = axes;
    set(ha2, 'units', 'pixels', 'position', [0 0 vid.Width/ScnScale vid.Height/ScnScale], 'nextplot', 'add', ...
        'xlim', [0 vid.Width], 'ylim', [0 vid.Height], 'ydir', 'reverse', 'XColor', 'none', 'YColor', 'none', 'Color', 'none', ...
        'Colormap', mycolormap);

    %%
    F = struct('cdata', [], 'colormap', []);

    img = imagesc(ha1, read(vid, traceL(1, 5)));

    p = plot(ha2, [traceL(1, 1) traceR(1, 1)], [traceL(1, 2) traceR(1, 2)], '-k', 'LineWidth', 1.5);

    k = -1 * diff(p.XData) / diff(p.YData);

    if ~isinf(k)
        dx = lineLength / sqrt(1 + k^2);

        x1 = [mean(p.XData) mean(p.XData)] + [dx -dx];
        y1 = [mean(p.YData) mean(p.YData)] + [dx -dx]*k;
    else
        x1 = [mean(p.XData) mean(p.XData)];
        y1 = [mean(p.YData) mean(p.YData)] + [100 -100];
    end

    angle_level = (obj.AngleHead{this_id}(1)-min(cRange)) / range(cRange);
    angle_level = min([1 angle_level]);
    angle_level = max([0 angle_level]);

    color_index = round(angle_level * size(mycolormap, 1));

    this_color = mycolormap(color_index, :);

    d(1) = patch(ha2, 'XData', [x1 nan], 'YData', [y1 nan], 'LineWidth', 8, 'EdgeColor', this_color, 'EdgeAlpha', 0.8);

    F(1) = getframe(fig);

    for i = 2:NumFrames

        img.CData = read(vid, traceL(i, 5));

        p.XData = [traceL(i, 1) traceR(i, 1)];
        p.YData = [traceL(i, 2) traceR(i, 2)];

        k = -1 * diff(p.XData) / diff(p.YData);

        if ~isinf(k)
            dx = lineLength / sqrt(1 + k^2);

            x1 = [mean(p.XData) mean(p.XData)] + [dx -dx];
            y1 = [mean(p.YData) mean(p.YData)] + [dx -dx]*k;
        else
            x1 = [mean(p.XData) mean(p.XData)];
            y1 = [mean(p.YData) mean(p.YData)] + [100 -100];
        end

        angle_level = (obj.AngleHead{this_id}(i)-min(cRange)) / range(cRange);
        angle_level = min([1 angle_level]);
        angle_level = max([0 angle_level]);

        color_index = round(angle_level * size(mycolormap, 1));
        if color_index == 0
            color_index = 1;
        end

        this_color = mycolormap(color_index, :);

        arrayfun(@(x) set(x, 'EdgeAlpha', x.EdgeAlpha*fade_ratio), d);
        d(i) = patch(ha2, 'XData', [x1 nan], 'YData', [y1 nan], 'LineWidth', 2, 'EdgeColor', this_color, 'EdgeAlpha', 0.8);

        F(i) = getframe(fig);
    end

    %
    arrayfun(@(x) set(x, 'EdgeAlpha', 0.5), d); 

    F(end+1) = getframe(fig);
    for i = 1:10
        F(end+1) = F(end);
    end
    
    %%

    close(fig);

    writerObj = VideoWriter(VidLabelFileName);
    writerObj.FrameRate = 20;
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
    clear writerObj F vid d traceL traceR

end
