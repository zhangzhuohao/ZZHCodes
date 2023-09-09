%%
VideoFolder = uigetdir("E:\YuLab\Work\GPS\Video\", "Choose target vedio folder");
if ~VideoFolder
    return;
end

matFiles = dir(VideoFolder + "\*.mat");

%%
opts = delimitedTextImportOptions("NumVariables", 28);

opts.DataLines = [4, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["frame", "ear_tip_left_x", "ear_tip_left_y", "ear_tip_left_p", "ear_base_left_x", "ear_base_left_y", "ear_base_left_p", "ear_tip_right_x", "ear_tip_right_y", "ear_tip_right_p", "ear_base_right_x", "ear_base_right_y", "ear_base_right_p", "pattern_left_x", "pattern_left_y", "pattern_left_p", "pattern_right_x", "pattern_right_y", "pattern_right_p", "port_left_x", "port_left_y", "port_left_p", "port_right_x", "port_right_y", "port_right_p", "port_center_x", "port_center_y", "port_center_p"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%%
NumTrials = length(matFiles);

AngleHead  = cell(NumTrials, 1);
AngleHeadS = cell(NumTrials, 1);
% Import the tracking data

for i = 1:NumTrials
    
    csvFile = dir(strcat(VideoFolder, filesep, matFiles(i).name(1:end-4), "*00.csv"));

    thisTrack = readtable(fullfile(VideoFolder, csvFile.name), opts);

    PortC = [mean(thisTrack.port_center_x(thisTrack.port_center_p>.999)) mean(thisTrack.port_center_y(thisTrack.port_center_p>.999))];
    PortL = [mean(thisTrack.port_left_x(thisTrack.port_left_p>.999))     mean(thisTrack.port_left_y(thisTrack.port_left_p>.999))];
    PortR = [mean(thisTrack.port_right_x(thisTrack.port_right_p>.999))   mean(thisTrack.port_right_y(thisTrack.port_right_p>.999))];
    
    vecPort = PortR - PortL;

    nframes = height(thisTrack);
    angle_head = zeros(1, nframes);
    for t = 1:nframes
        if thisTrack.ear_base_left_p(t) > .99 && thisTrack.ear_base_right_p(t) > .99

            ear_left  = [thisTrack.ear_base_left_x(t)  thisTrack.ear_base_left_y(t)];
            ear_right = [thisTrack.ear_base_right_x(t) thisTrack.ear_base_right_y(t)];
            vec_ear   = ear_right - ear_left;

            if vec_ear(2)-vecPort(2) <= 0
                angle_sign = 1; % face to right port
            else
                angle_sign = -1; % face to left port
            end

            angle_head(t) = acos(dot(vec_ear, vecPort) / (norm(vec_ear) * norm(vecPort))) * 180 / pi;
            angle_head(t) = angle_head(t) * angle_sign;
        else
            angle_head(t) = nan;
        end
    end
    angle_head_s = smoothdata(angle_head, "gaussian", 5, "omitnan");

    AngleHead{i}  = angle_head;
    AngleHeadS{i} = angle_head_s;
end

%%
minFrame = min(cellfun(@(x) length(x), AngleHead));
AngleHeadMat  = cell2mat(cellfun(@(x) x(1:minFrame), AngleHead,  "UniformOutput", false));
AngleHeadMat(isnan(AngleHeadMat)) = 0;
AngleHeadSMat = cell2mat(cellfun(@(x) x(1:minFrame), AngleHeadS, "UniformOutput", false));
AngleHeadSMat(isnan(AngleHeadSMat)) = 0;

%%
[~, sortid] = sort(AngleHeadMat(:,end));
mycolormap = customcolormap_preset("red-white-blue");
figure(); colormap(mycolormap);
imagesc(AngleHeadSMat(sortid, :));
caxis([-80 80]);
cb = colorbar(); cb.Limits = [-80 80];

%%
plot(mean(AngleHeadSMat(AngleHeadSMat(:, end)>0, :)))
hold on;
plot(mean(AngleHeadSMat(AngleHeadSMat(:, end)<0, :)))

