function MTSorted = getMovementTimeSorted(obj)

MTCollected = cell(1, length(obj.Ports));
for i = 1:length(obj.Ports)
    ind_thisFP = find(obj.PortCorrect==obj.Ports(i) & strcmp(obj.Outcome, 'Correct'));
    iMTs = obj.MovementTime(ind_thisFP);
    MTCollected{i} = iMTs;
end
MTSorted = MTCollected;

end

