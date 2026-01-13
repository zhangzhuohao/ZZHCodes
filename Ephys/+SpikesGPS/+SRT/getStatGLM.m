function StatGLM = getStatGLM(r)

rs = struct();
rs.Units = r.Units;
rs.EphysTable = r.EphysTable;

% Ramping
num_units = size(rs.Units.SpikeNotes, 1);
for ku = 1:num_units
    stat_k = SpikesGPS.SRT.testRampingGLM(rs, ku);
    if ku==1
        stat_ramping = stat_k;
    else
        stat_ramping = [stat_ramping; stat_k];
    end
end

StatGLM.ramping = stat_ramping;

save_name = sprintf('StatGLM_%s_%s.mat', r.BehaviorClass.Subject, r.BehaviorClass.Session);
save(save_name, 'StatGLM');

end