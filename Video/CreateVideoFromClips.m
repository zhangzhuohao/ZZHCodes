% Will use all images in the current folder to make videos
%%
View = "Field";

clip_files = dir("./*.avi");
n_clips = length(clip_files);

%% Create videos from video

VidClipName = "./GatheredVideo.avi";

writerObj = VideoWriter(VidClipName);
writerObj.FrameRate = 25; % this is 2 x slower
writerObj.Quality = 100;


% set the seconds per image
% open the video writer

open(writerObj);
for ii=1:n_clips
    disp(ii);
    clip = VideoReader(clip_files(ii).name);
    frames = read(clip);
    writeVideo(writerObj, frames);
end
close(writerObj);

fprintf("\nGathered video created.\n");