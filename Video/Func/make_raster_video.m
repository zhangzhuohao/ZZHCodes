function make_raster_video(r, units, trials, range, video_dir, gain, t_query)
% Kin.make_raster_video(r, units, trials, range_press, video_dir)

if nargin < 7 || isempty(t_query)
    t_query = 0;
end

t_query = t_query(:)';   % ensure row vector
n_query = length(t_query);
did_save_frame = false(1, n_query);

if nargin < 6
    gain = 2.5;
end


% --- find index of units ---
n_units = length(units);
spike_times = cell(1, n_units);
for i =1:length(units)
    unit = units{i};
    out = parse_unit_name(unit);
    ind = find(r.Units.SpikeNotes(:, 1)==out.ch & r.Units.SpikeNotes(:, 2)==out.unit);
    spike_times{i} = r.Units.SpikeTimes(ind).timings;
end


out_file = fullfile(pwd, sprintf('raster_video_%s.mp4', out.anm_session));

vw = VideoWriter(out_file, 'MPEG-4');
vw.FrameRate = 20;
open(vw);

event_tab = r.EventTable;

% --- find video ---
n_trials = length(trials);
vid = repmat(struct('top', [], 'side', [], 'ind_plot', []), 1, n_trials);
meta = repmat(struct('top', [], 'side', []), 1, n_trials);

FP = zeros(1, n_trials);
rt = zeros(1, n_trials);
poke = zeros(1, n_trials);

for i =1:length(trials)

    i_trial = trials(i);
    beh_row = event_tab(find(abs(event_tab.t_press-i_trial)<10),:);
    FP(i) = beh_row.FP;
    rt(i) = beh_row.rt;
    poke(i) = beh_row.t_poke-beh_row.t_press;

    top_vid = fullfile(video_dir.top, sprintf('Top_Press_%d.mp4', i_trial));
    vid(i).top = VideoReader(top_vid);

    top_vid_meta = fullfile(video_dir.top, sprintf('Top_Press_%d.mat', i_trial));
    meta(i).top = load(top_vid_meta);

    side_vid = fullfile(video_dir.side, sprintf('Side_Press_%d.mp4', i_trial));
    vid(i).side = VideoReader(side_vid);

    side_vid_meta = fullfile(video_dir.side, sprintf('Side_Press_%d.mat', i_trial));
    meta(i).side = load(side_vid_meta);

    time_top_rel = meta(i).top.VideoInfo.EphysTimeStamps - meta(i).top.VideoInfo.Time;
    time_side_rel = meta(i).side.VideoInfo.EphysTimeStamps - meta(i).side.VideoInfo.Time;

    vid(i).time_top_rel = time_top_rel;
    vid(i).time_side_rel = time_side_rel;

    % use TOP as master
    vid(i).ind_plot = find(time_top_rel >= range(1) & time_top_rel <= range(2));
end

% --- initiate two plots ---
fig_w = 23.5;
fig_h = 13.5;
vid_w = 5;

hf = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 fig_w fig_h], 'Visible', 'on');

width_side = vid_w;
height_side = vid(1).side.Height*width_side/vid(1).side.Width;

% --- show the first frame --
xp = 1.5;
yp = 4;
h_gap = 0.25;

ax1 = gobjects(n_trials,1);
himg_side = gobjects(n_trials,1);
htext_side = gobjects(n_trials,1);

for i =1:n_trials
    ax1(i) = axes('Parent', hf,...
        'Units', 'centimeters',...
        'XLim',[0 vid(i).side.Width], 'YLim',[0 vid(i).side.Height],...
        'Position', [xp+(i-1)*(h_gap+width_side), yp, width_side, height_side],...
        'NextPlot','add');

    side_frame = read(vid(i).side, 1);
    side_frame = uint8(min(double(side_frame) * gain, 255));
    himg_side(i) = imshow(side_frame);

end

width_top = width_side;
height_top = vid(1).top.Height*width_top/vid(1).top.Width;

% --- show the first frame ---
yp = yp + height_side + .1;
ax2 = gobjects(n_trials,1);
himg_top = gobjects(n_trials,1);

for i =1:n_trials
    ax2(i) = axes('Parent', hf,...
        'Units', 'centimeters',...
        'XLim',[0 vid(i).top.Width], 'YLim',[0 vid(i).top.Height],...
        'Position', [xp+(i-1)*(h_gap+width_top), yp, width_top, height_top],...
        'NextPlot','add');

    title(sprintf('Press_%d', trials(i)), 'Interpreter','none', 'FontSize',10)

    top_frame = read(vid(i).top, 1);
    top_frame = uint8(min(double(top_frame) * gain, 255));
    himg_top(i) = imshow(top_frame);

    htext_top(i) = text(0, vid(i).top.Height-40, ...
        sprintf('%d ms', round(meta(i).top.VideoInfo.EphysTimeStamps(1)-meta(i).top.VideoInfo.Time)), ...
        'Color', 'w', 'FontSize',10);

end

height_raster = 2;
yp = 1.5;

ax3 = gobjects(n_trials,1);

for i =1:n_trials
    ax3(i) = axes('Parent', hf,...
        'Units', 'centimeters',...
        'XLim', range, 'YLim',[0.5 n_units+.5],...
        'YTick',[.5:1:n_units],...
        'YGrid','on',...
        'FontSize',7,...
        'YTickLabel', num2cell(1:n_units),...
        'Position', [xp+(i-1)*(h_gap+width_top), yp, width_top, height_raster],...
        'NextPlot','add', 'YDir','reverse');

    if i == 1
        xlabel('ms');
        ylabel('units')
    else
        ax3(i).YTickLabel = [];
    end

    plotshaded([0 FP(i)], [0 0; n_units+1 n_units+1],[255, 200, 30]/255)
    xline(ax3(i), rt(i)+FP(i), 'Color','#F13E93', 'LineWidth',1);
    xline(ax3(i), poke(i), 'Color','#468432', 'LineWidth',1);
end

n_frames = min(arrayfun(@(v) length(v.ind_plot), vid));

h_line = gobjects(n_trials,1);
h_raster = gobjects(n_trials,1);

t_start = zeros(1,n_trials);
t_end = zeros(1,n_trials);

for i = 1:n_trials
    t_start(i) = meta(i).side.VideoInfo.Time + range(1);
    t_end(i) = meta(i).side.VideoInfo.Time + range(2);
    % --- raster ---

    time_sofar = [t_start(i), t_end(i)];

    X = [];
    Y = [];

    for j = 1:n_units

        spk = spike_times{j};
        spk = spk(spk >= time_sofar(1) & spk < time_sofar(2));

        if isempty(spk), continue; end

        % OPTIONAL: align to top instead of side
        % spk = spk - meta(i).top.VideoInfo.Time;

        spk = spk - meta(i).side.VideoInfo.Time;

        n = length(spk);

        xx = [spk; spk; nan(1,n)];
        yy = [(j-0.5)*ones(1,n);
            (j+0.5)*ones(1,n);
            nan(1,n)];

        X = [X, xx(:)'];
        Y = [Y, yy(:)'];
    end

    plot(ax3(i), X, Y, 'color', [.75 .75 .75], 'LineWidth',0.5);
end

for k = 1:n_frames
    for i = 1:n_trials

        ind_plot = vid(i).ind_plot;
        idx = ind_plot(k);

        % --- time ---
        top_time = meta(i).top.VideoInfo.EphysTimeStamps(idx);
        top_time_rel = top_time - meta(i).top.VideoInfo.Time;
        htext_top(i).String = sprintf('%d ms', round(top_time_rel));

        % --- vertical line ---
        if ~isgraphics(h_line(i))
            h_line(i) = xline(ax3(i), top_time_rel, 'Color',[.1 .1 .8], 'LineWidth',1);
        else
            h_line(i).Value = top_time_rel;
        end

        % --- side video ---
        [dt, ind_side] = min(abs(meta(i).side.VideoInfo.EphysTimeStamps - top_time));

        if ind_side <= vid(i).side.NumFrames && meta(i).side.VideoInfo.EphysTimeStamps(end) > top_time && dt < 5
            side_frame = read(vid(i).side, ind_side);
            side_frame = uint8(min(double(side_frame) * gain, 255));
            himg_side(i).CData = side_frame;

            side_time_rel = meta(i).side.VideoInfo.EphysTimeStamps(ind_side) - meta(i).side.VideoInfo.Time;

            % OPTIONAL: multiply by 1000 if seconds
            % htext_side(i).String = sprintf('%d ms', round(side_time_rel*1000));

        end
  
        % --- top video ---
        [~, ind_top] = min(abs(meta(i).top.VideoInfo.EphysTimeStamps - top_time));
        top_frame = read(vid(i).top, ind_top);
        top_frame = uint8(min(double(top_frame) * gain, 255));
        himg_top(i).CData = top_frame;

        % --- raster ---

        time_sofar = [t_start(i), top_time];

        X = [];
        Y = [];

        for j = 1:n_units

            spk = spike_times{j};
            spk = spk(spk >= time_sofar(1) & spk < time_sofar(2));

            if isempty(spk), continue; end

            % OPTIONAL: align to top instead of side
            % spk = spk - meta(i).top.VideoInfo.Time;

            spk = spk - meta(i).side.VideoInfo.Time;

            n = length(spk);

            xx = [spk; spk; nan(1,n)];
            yy = [(j-0.5)*ones(1,n);
                  (j+0.5)*ones(1,n);
                  nan(1,n)];

            X = [X, xx(:)'];
            Y = [Y, yy(:)'];
        end

        if ~isempty(X)
            if ~isgraphics(h_raster(i))
                h_raster(i) = plot(ax3(i), X, Y, 'k', 'LineWidth',1);
            else
                h_raster(i).XData = X;
                h_raster(i).YData = Y;
            end
        end

        % --- still frame ---
        for q = 1:n_query

            if ~did_save_frame(q) && abs(top_time_rel - t_query(q)) < 9

                frame = getframe(hf);

                fname = sprintf('raster_still_%dms.png', round(t_query(q)));
                imwrite(frame.cdata, fullfile(pwd, fname));

                did_save_frame(q) = true;

            end

        end

    end

    drawnow
    frame = getframe(hf);
    writeVideo(vw, frame);
end

close(vw);