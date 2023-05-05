function [tLEDon, FigLEDon] = FindLEDonGPS(tsROI, SummedROI)

% Jianing Yu
% 4/27/2021
% 4/15/2022
% Track the change in ROI intensity to judge when LED is on
% Revised version applies to GPS system

led = SummedROI;

FigLEDon = figure(14); clf(14)
set(gcf, 'name', 'ROI dynamics', 'units', 'centimeters', 'position', [15 5 25 15])
led = detrend(led);
led = led - median(led);
ha1 = subplot(2, 2, 1);
set(ha1, 'nextplot', 'add', 'ylim', [-2*10^4  5*10^4])
plot(ha1, tsROI, led, 'k');

ha2 = subplot(2, 2, 3);
set(ha2, 'nextplot', 'add')
histogram(ha2, led, 'BinWidth', 0.01);
axis 'auto y'

% select threshold
fprintf('\nSelect threshold, end selection by right click\n')
[x_thrh, ~] = getpts(ha2);

roi_th = min(x_thrh); % this is the threshold to extract LED_on times
xline(ha2, roi_th, 'color', 'g', 'linestyle', ':', 'linewidth', 1)

yline(ha1, roi_th, 'color', 'g', 'linestyle', ':', 'linewidth', 1)

% find those ROIs that are above threshold
above_th = find(led > roi_th);

% find begs
above_th_begs = above_th([1; 1+find(diff(above_th)>1)]);
above_th_ends = above_th([find(diff(above_th)>1); end]);
plot(ha1, tsROI(above_th_begs), led(above_th_begs), 'go');
plot(ha1, tsROI(above_th_ends), led(above_th_ends), 'bo');

% check the duration of these LED on periods
above_th_dur = tsROI(above_th_ends) - tsROI(above_th_begs);

% check the amplitudes
amp_above = zeros(length(above_th_dur), 1);
for m = 1:length(amp_above)
    m_led = led(above_th_begs(m):above_th_ends(m));
    amp_above(m) = median(m_led);
end

ha3 = subplot(2, 2, 2);
above_th_dur_log = log(above_th_dur);
set(ha3, 'nextplot', 'add', 'xscale', 'linear')
plot(ha3, above_th_dur_log, amp_above, 'ko')

title('Select multiple points for decision boundary')
xlabel('log(dur)')
ylabel('Amp')

% remove "bad" ROIs
disp('Select mulitple points to define a decision boundary, end selection by right click')
[x_thrh, y_thrh] = getpts(ha3);

plot(ha3, x_thrh, y_thrh, 'bx')

tbl = table(x_thrh(1:end-1), y_thrh(1:end-1));
lm = fitlm(tbl, 'linear'); % y = ax + b

ypredict = predict(lm, above_th_dur_log);
[above_sort, indsort] = sort(above_th_dur_log);

plot(above_sort, ypredict(indsort), 'Color','m');

% find those that are above this line
ind_good = find(amp_above > ypredict);

falseindex_LEDon = find(amp_above < ypredict);
plot(ha3, above_th_dur_log(ind_good), amp_above(ind_good), 'g+', 'linewidth', 2)
axis tight

plot(ha1, tsROI(above_th_begs(falseindex_LEDon)), led(above_th_begs(falseindex_LEDon)), 'r*', 'markersize', 8);
plot(ha1, tsROI(above_th_ends(falseindex_LEDon)), led(above_th_ends(falseindex_LEDon)), 'r*', 'markersize', 8);

if above_th_begs(1) == 1
    falseindex_LEDon = [1; falseindex_LEDon];
end

if above_th_ends(end) == length(tsROI)
    falseindex_LEDon = [falseindex_LEDon; length(above_th_ends)];
end

above_th_begs(falseindex_LEDon) = [];
above_th_ends(falseindex_LEDon) = [];
above_th_dur(falseindex_LEDon)  = [];

plot(ha1, tsROI(above_th_begs), led(above_th_begs), 'g+', 'linewidth', 2);

% empirically, push the onset time by one frame
% above_th_begs = above_th_begs-1;

ha4=subplot(2, 2, 4);
set(ha4, 'nextplot', 'add', 'xlim', [-1 5]);

abv_seg =[];

for j = 1:length(above_th_begs)
    abv_seg{j} = led(above_th_begs(j):above_th_ends(j));
    plot(ha4, abv_seg{j}, 'k')
end

abv_seg =[];

cla(ha4);
for j = 1:length(above_th_begs)
    abv_seg{j} = led(above_th_begs(j)-3:above_th_ends(j));
    tseg = (0:length(abv_seg{j})-1)-3;
    plot(ha4, tseg, abv_seg{j}, 'k')
    plot(ha1, tsROI(above_th_begs(j)-3:above_th_ends(j)), abv_seg{j}, 'g')
end

xlabel('Frames')
ylabel('ROI')

tLEDon = tsROI(above_th_begs);
