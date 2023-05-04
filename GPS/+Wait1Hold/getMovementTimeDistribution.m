function MovementTimeDistribution = getMovementTimeDistribution(obj)

datain = obj.MovementTime;
datain(isnan(datain))=[];
[data2575] = prctile(datain, [25, 75]);
interq = data2575(2) - data2575(1);
c = 5;

binEdges = 0:0.05:3;
binCenters = (binEdges(1:end-1)+binEdges(2:end))/2;

% Only include correct trials
MT_Port1 = obj.MovementTime(obj.PortCorrect == 1 & strcmp(obj.Outcome, 'Correct'));
MT_Port1(MT_Port1>data2575(2)+interq*c | MT_Port1<data2575(1)-interq*c) = [];
MT_Port2 = obj.MovementTime(obj.PortCorrect == 2 & strcmp(obj.Outcome, 'Correct'));
MT_Port2(MT_Port2>data2575(2)+interq*c | MT_Port2<data2575(1)-interq*c) = [];

Ncounts = histcounts(MT_Port1,binEdges);
Counts_Port1 = Ncounts';
Ncounts = histcounts(MT_Port2,binEdges);
Counts_Port2 = Ncounts';
MT_Centers = binCenters';

MovementTimeDistribution = table(MT_Centers, Counts_Port1, Counts_Port2);

end

