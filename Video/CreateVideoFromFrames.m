% Will use all images in the current folder to make videos

imagefiles = dir('*.jpg');
nfiles = length(imagefiles);    % Number of files found

VidClipName = 'RefineVideo.avi';

writerObj = VideoWriter(VidClipName);
writerObj.FrameRate = 25; % this is 2 x slower

% set the seconds per image
% open the video writer

open(writerObj);
for ii=1:nfiles
    currentfilename = imagefiles(ii).name;
    currentimage = imread(currentfilename);
    writeVideo(writerObj, currentimage);
end
close(writerObj);