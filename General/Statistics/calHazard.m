function hazard_ratio = calHazard(hold_dur, fp, bin_edges)

hazard_ratio = zeros(length(bin_edges), 1);
for i = 1:length(bin_edges)
    bin_i = bin_edges(i);
    if i==1
        n_pre = sum(hold_dur<bin_i);
        n_tot = length(hold_dur);
    else
        n_pre = sum(hold_dur<bin_i & hold_dur>=bin_edges(i-1) & (fp>=bin_i | fp==-1));
        n_tot = sum(hold_dur>=bin_edges(i-1) & (fp>=bin_i | fp==-1));
    end
    hazard_ratio(i) = n_pre / n_tot;
end
end % calHazard
