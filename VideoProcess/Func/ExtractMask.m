function [maskout, fig_mask] = ExtractMask(vidFile, frmindex)

% Jianing Yu 
% 4/27/2021
% Draw mask from selected frame in vidFile
% frmindex gives the beginning and end of frame selection

maskout = [];

%% extract ROI that matters
% read 100 frames:

filename = vidFile;
% Construct a multimedia reader object associated with file
vidObj = VideoReader(filename);

% read one min of data
% framebeg = 1000;
% frameend = framebeg + 60*50;
list_of_frames = read(vidObj, frmindex);
list_of_frames = squeeze(list_of_frames(:, :, 1, :));
maxproj_frames = max(list_of_frames, [], 3);

% plot the max projection of these frames:
fig_mask = figure(13); clf
set(gcf, 'name', 'ROI selection', 'units', 'centimeters', 'position', [15 5 20 20])
imagesc(maxproj_frames, [0 400]);
colormap('gray');
axis equal;
axis off;

fprintf('\nPlease select region of interests\n');

roi_selected = drawfreehand();
maskout = createMask(roi_selected); % this mask determines what pixels are included. this is the mask to use in the future.
