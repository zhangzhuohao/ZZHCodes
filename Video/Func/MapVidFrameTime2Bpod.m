function tframes_in_b = MapVidFrameTime2Bpod(tLEDon, tbeh_trigger, tbeh_lag, tsROI)

% Jianing Yu
% 5/1/2021
% given LED time, indout, and trigger time in b, conver the time of each
% frame to a time in behavior domain

tframes_in_b = zeros(length(tsROI), 1); % this is the frame time in behavior time domain

nLEDon = length(tLEDon);
if nLEDon==1
    tframes_in_b = tsROI - tLEDon;
else
    for i = 1:nLEDon
        if i==1
            frames_sofar = find(tsROI<=tLEDon(i+1));
        elseif i==length(tLEDon)
            frames_sofar = find(tsROI>tLEDon(i)-tbeh_lag(i)-500);
        else
            frames_sofar = find(tsROI>tLEDon(i)-tbeh_lag(i)-500 & tsROI<=tLEDon(i+1));
        end
        tframes_in_b(frames_sofar) = tsROI(frames_sofar) - tLEDon(i) + tbeh_trigger(i); % convert the frame time to time in the behavior domain
    end
end