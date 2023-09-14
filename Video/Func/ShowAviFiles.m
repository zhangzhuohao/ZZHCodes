function [vidFiles, tsFiles, cfolder]=ShowAviFiles(cfolder, view)


if nargin<1
    [~, path, ~] = uigetfile({'*.avi';'*.txt'},...
        'Select any file');

    cfolder = path;
end

allVidFiles = dir(fullfile(cfolder, 'Cam*.avi'));
if ~isempty(allVidFiles)
    vidFiles = arrayfun(@(x)x.name, allVidFiles, 'UniformOutput',0);
    vidFileTime = cell2mat(arrayfun(@(x)x.datenum, allVidFiles, 'UniformOutput',0));
    [~, indsort] = sort(vidFileTime);
    vidFiles = vidFiles(indsort);
end

alltsFiles = dir(fullfile(cfolder, 'Cam*.txt'));

if ~isempty(alltsFiles)
    tsFiles = arrayfun(@(x)x.name, alltsFiles, 'UniformOutput',0);
    tsFileTime = cell2mat(arrayfun(@(x)x.datenum, alltsFiles, 'UniformOutput',0));
    [~, indsort] = sort(tsFileTime);
    tsFiles = tsFiles(indsort);
end

% Sometimes the recording program interrupts early, causing the last video file to have no corresponding timestamp file
switch view
    case {'Top', 'Front'}
        if length(vidFiles)>length(tsFiles)
            fprintf("\nVideo files do not match to timestamp files, ")
            vidFiles(end) = [];
        end
    case {'Field', 'Init'}
        if length(vidFiles)>length(tsFiles)
            fprintf("\nVideo files do not match to timestamp files, ")
        end
end
sprintf('%s',cfolder);