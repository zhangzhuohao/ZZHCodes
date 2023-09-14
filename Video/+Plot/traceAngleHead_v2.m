function traceAngleHead_v2(ax, obj, ForePeriod, AlignTo, opts)

arguments
    ax  (1, 1) handle
    obj (1, 1) 
    ForePeriod (1, 1) {mustBeMember(ForePeriod, ["Short", "Med", "Long"])} = "Med"
    AlignTo (1, 1) {mustBeMember(AlignTo, ["In", "Out"])} = "In"
    opts.color (1, 1) GPSColor = GPSColor()
    opts.PortChosen (1, 1) {mustBeMember(opts.PortChosen, ["L", "R", "Both"])} = "Both"
    opts.Performance (1, 1) {mustBeMember(opts.Performance, ["Correct", "Wrong", "Late", "Premature", "Valid"])} = "Valid"
    opts.SortBy (1, 1) {mustBeMember(opts.SortBy, ["HD", "MT", "FP"])} = "HD"
    opts.Label (1, 1) {mustBeMember(opts.Label, ["All", "None", "Saline", "Control", "Chemo"])} = "All"
end

ind = obj.Ind;
info = obj.TrialInfo;

%% extract by performance
switch opts.Performance
    case {'Premature', 'Late'}
        switch  ForePeriod
            case {'All'}
                these_trials = find(ind.(opts.Performance));
            otherwise
                these_trials = find(ind.(opts.Performance) & ind.( ForePeriod));
        end
        ax.Title.String = opts.Performance;
    case {'Valid'}
        switch  ForePeriod
            case {'All'}
                switch opts.PortChosen
                    case {'Both'}
                        these_trials = find(ind.Wrong | ind.Correct);
                    otherwise
                        these_trials = find(ind.("Choose"+opts.PortChosen) & (ind.Wrong | ind.Correct));
                end
            otherwise
                switch opts.PortChosen
                    case {'Both'}
                        these_trials = find((ind.Wrong | ind.Correct) & ind.( ForePeriod));
                    otherwise
                        these_trials = find(ind.("Choose"+opts.PortChosen) & (ind.Wrong | ind.Correct) & ind.( ForePeriod));
                end
        end
        yline(ax, sum(ind.("Choose"+opts.PortChosen) & ind.Correct) + 0.5, ':', 'LineWidth', 1);
%         ax.Title.String = "Chose Port"+opts.PortChosen;
%         ax.Title.Color  = opts.color.("Port"+opts.PortChosen);
    otherwise
        switch  ForePeriod
            case {'All'}
                these_trials = find(ind.("Choose"+opts.PortChosen) & ind.(opts.Performance));
            otherwise
                these_trials = find(ind.("Choose"+opts.PortChosen) & ind.(opts.Performance) & ind.( ForePeriod));
        end
%         ax.Title.String = "Chose Port"+opts.PortChosen;
%         ax.Title.Color  = opts.color.("Port"+opts.PortChosen);
end







