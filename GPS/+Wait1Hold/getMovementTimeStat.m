function MovementTimeStat = getMovementTimeStat(obj)

datain = obj.MovementTime;
datain(isnan(datain))=[];
[data2575] = prctile(datain, [25, 75]);
interq = data2575(2) - data2575(1);
c = 5;

% Only include correct trials
MT_Port1 = obj.MovementTime(obj.PortCorrect==1 & strcmp(obj.Outcome, 'Correct'));
MT_Port1(MT_Port1>data2575(2)+interq*c | MT_Port1<data2575(1)-interq*c) = [];
iMTOut = calRT(MT_Port1*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);
MedianMT = iMTOut.median*0.001;
MedianMT_ksdensity = iMTOut.median_ksdensity*0.001;
N = length(MT_Port1);
Port = 1;

MT_Port2 = obj.MovementTime(obj.PortCorrect ==2 & strcmp(obj.Outcome, 'Correct'));
MT_Port2(MT_Port2>data2575(2)+interq*c | MT_Port2<data2575(1)-interq*c) = [];
iMTOut = calRT(MT_Port2*1000, [], "Remove100ms", 0, "RemoveOutliers", 0, 'CalSE', 0);
MedianMT = [MedianMT; iMTOut.median*0.001];
MedianMT_ksdensity = [MedianMT_ksdensity; iMTOut.median_ksdensity*0.001];
N = [N; length(MT_Port2)];
Port = [Port; 2];

MovementTimeStat = table(Port, N, MedianMT, MedianMT_ksdensity);

end
