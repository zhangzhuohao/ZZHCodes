% Will use all images in the current folder to make videos
%%
View = "Field";

video_files = dir("./*.avi");


%% Create videos from video

VidClipName = fullfile(GatherFolder, "RefineVideo.avi");

writerObj = VideoWriter(VidClipName);
writerObj.FrameRate = 25; % this is 2 x slower
writerObj.Quality = 100;

% set the seconds per image
% open the video writer

open(writerObj);
for ii=1:n_bad_labels
    currentfilename = fullfile(GatherFolder, bad_label_frames(ii).name);
    currentimage = imread(currentfilename);
    writeVideo(writerObj, currentimage);
end
close(writerObj);

fprintf("\nRefine video created.\n");