function IndSort = VisualizePopPSTH(Pop)
% Pop = PopOut;
% 5/14/2023 revised JY

% IndSort = PSTHOut.IndSort;
% IndSignificant = PSTHOut.IndSignificant;
% nsig = sum(IndSignificant);

% visualize PSTH (z score)
zrange     = [-4 4];
pval_pop   = zeros(1, size(Pop.Units, 1));
tpeaks_pop = zeros(1, size(Pop.Units, 1));
n_unit     = size(Pop.Units, 1);

[nFP, nPort] = size(Pop.CentIn);
MixedFPs     = Pop.FPs * 1000;

press_col   = [5 191 219]/255;
trigger_col = [242 182 250]/255;
release_col = [87, 108, 188]/255;
reward_col  = [164, 208, 164]/255;
FP_cols     = [255, 217, 90; 192, 127, 0; 76, 61, 61]/255;
%
IndSort = Pop.IndSort;
nsig = n_unit - length(Pop.IndUnmodulated);

CentInPSTHs  = cell(size(Pop.CentIn));
CentInPSTHZs = cell(size(Pop.CentIn));
tCentInPSTHs = cell(size(Pop.CentIn));
CentOutPSTHZs = cell(size(Pop.CentOut));
CentOutPSTHs  = cell(size(Pop.CentOut));
tCentOutPSTHs = cell(size(Pop.CentOut));
TriggerPSTHs  = cell(size(Pop.CentOut));
TriggerPSTHZs = cell(size(Pop.CentOut));
tTriggerPSTHs = cell(size(Pop.CentOut));
RewardPSTHs  = cell(size(Pop.Reward));
RewardPSTHZs = cell(size(Pop.Reward));
tRewardPSTHs = cell(size(Pop.Reward));

%%
for ifp = 1:nFP
    for jport = 1:nPort
        thisPSTH                 = Pop.CentInZ{ifp, jport}(2:end, :);
        CentInPSTHZs{ifp, jport} = thisPSTH(IndSort{jport}, :);
        thisPSTH                 = Pop.CentIn{ifp, jport}(2:end, :);
        CentInPSTHs{ifp, jport}  = thisPSTH(IndSort{jport}, :);
        tCentInPSTHs{ifp, jport} = Pop.CentInZ{ifp, jport}(1, :);

        thisPSTH                  = Pop.CentOutZ{ifp, jport}(2:end, :);
        CentOutPSTHZs{ifp, jport} = thisPSTH(IndSort{jport}, :);
        thisPSTH                  = Pop.CentOut{ifp, jport}(2:end, :);
        CentOutPSTHs{ifp, jport}  = thisPSTH(IndSort{jport}, :);
        tCentOutPSTHs{ifp, jport} = Pop.CentOutZ{ifp, jport}(1, :);

        thisPSTH                  = Pop.TriggerZ{ifp, jport}(2:end, :);
        TriggerPSTHZs{ifp, jport} = thisPSTH(IndSort{jport}, :);
        thisPSTH                  = Pop.Trigger{ifp, jport}(2:end, :);
        TriggerPSTH{ifp, jport}   = thisPSTH(IndSort{jport}, :);
        tTriggerPSTHs{ifp, jport} = Pop.TriggerZ{ifp, jport}(1, :);

        thisPSTH                       = Pop.RewardZ{ifp, jport}(2:end, :);
        RewardPSTHZs{ifp, jport} = thisPSTH(IndSort{jport}, :);
        thisPSTH                       = Pop.Reward{ifp, jport}(2:end, :);
        RewardPSTHs{ifp, jport}  = thisPSTH(IndSort{jport}, :);
        tRewardPSTHs{ifp, jport} = Pop.Reward{ifp, jport}(1, :);
    end
end

%% Concatenate raw PSTH for computing normalized PSTHs
PSTH_Concatenate = cell(nFP, nPort);
PSTH_ConNorm = cell(nFP, nPort);
IndCentIn =  cell(nFP, nPort);
IndCentOut = cell(nFP, nPort);
IndTrigger = cell(nFP, nPort);
IndReward = cell(nFP, nPort);
norm_range = [0 1];

for ifp = 1:nFP
    for jport = 1:nPort
        PSTH_Concatenate{ifp, jport} = [CentInPSTHs{ifp, jport} CentOutPSTHs{ifp, jport} RewardPSTHs{ifp, jport}];
        PSTH_ConNorm{ifp, jport}     = normalize(PSTH_Concatenate{ifp, jport}' , 'range', norm_range);
        PSTH_ConNorm{ifp, jport}     = PSTH_ConNorm{ifp, jport}';
        IndCentIn{ifp, jport}        = 1:size(CentInPSTHs{ifp, jport}, 2);
        IndCentOut{ifp, jport}       = IndCentIn{ifp, jport}(end) + (1:size(CentOutPSTHs{ifp, jport}, 2));
        IndReward{ifp, jport}        = IndCentOut{ifp, jport}(end) + (1:size(RewardPSTHs{ifp, jport}, 2));
    end
end

%% Plot
set_matlab_default;
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
% visualize PSTH
hf=48;
figure(hf); clf(hf)
set(gcf, 'unit', 'centimeters', 'position', [2 2 21 20], 'paperpositionmode', 'auto' ,'color', 'w')
size_factor = 1;
space = 0.1;

% A. Plot range normalize
xlevel_start = 1.25;
xlevel_now   = xlevel_start;
for jport = 1:nPort
    ylevel_start = 2;
    map_height = 3;

    % CentIn
    tRange = [tCentInPSTHs{ifp, jport}(1) 2000];
    Width  = size_factor*diff(tRange)/2000;
    for ifp = 1:nFP
        ha_centin(ifp) = axes('unit', 'centimeters', 'position',...
            [xlevel_now ylevel_start+(map_height+0.5)*(ifp-1) Width map_height], 'nextplot', 'add',...
            'xlim', tRange, 'xtick', -3500:500:2000,...
            'ylim', [0.5 n_unit+0.5], 'ydir', 'reverse', 'ticklength', [0.025 0.01], 'XTickLabelRotation', 90);
        if ifp > 1
            set(ha_centin(ifp), 'xticklabel', [])
        else
            ylabel('Units')
        end
        if jport > 1
            yticklabels([]);
            ylabel([]);
        end
        CentInPSTH_thisFP = PSTH_ConNorm{ifp, jport}(:, IndCentIn{ifp, jport});
        himage1 = imagesc(tCentInPSTHs{ifp, jport}, 1:n_unit, CentInPSTH_thisFP, norm_range);
        colormap('Parula')
        yrange = [0.5 n_unit+0.5];
        line([0 0], yrange, 'color', press_col, 'linestyle', ':', 'linewidth', 1.5);
        line([0 0]+MixedFPs(ifp), yrange, 'color', 'm', 'linestyle', ':', 'linewidth', 1.5)
        %         plotshaded([get(gca, 'xlim')], [nsig+.5 nsig+.5; n_unit+0.5 n_unit+0.5], 'w', 0.75);
        fill([get(gca, 'xlim') flip(get(gca, 'xlim'))], [nsig+.5 nsig+.5 n_unit+0.5 n_unit+0.5], 'w', 'FaceColor', 'w', 'EdgeColor', 'none', 'FaceAlpha', .75);
        if jport==1
            title(sprintf('FP=%2.0dms, Ntrials=%2.0d', MixedFPs(ifp), Pop.Trials(ifp, jport)));
        else
            title(sprintf('Ntrials=%2.0d', Pop.Trials(ifp, jport)));
        end
    end
    ylevel_now = ylevel_start+(map_height+0.5)*(nFP-1)+map_height+0.5;
    if jport==1
        uicontrol('Style','text','Units','centimeters','Position',[xlevel_start-1 ylevel_now  6 0.7],...
            'string', 'A. Normalized activity ([0-1])', ...
            'FontName','Dejavu Sans',  'fontweight', 'bold','fontsize', 9,'BackgroundColor',[1 1 1],...
            'HorizontalAlignment','Left');
    end
    xlevel_now = xlevel_now + Width + space;

    % Release
    ylevel_start = 2;
    map_height = 3;
    tRange = [-500 tCentOutPSTHs{1}(end)];
    Width  = size_factor*diff(tRange)/2000;
    for ifp =1:nFP
        ha_release(ifp) = axes('unit', 'centimeters', 'position',...
            [xlevel_now ylevel_start+(map_height+0.5)*(ifp-1)  Width map_height],...
            'nextplot', 'add', 'xlim', tRange, 'xtick', -500:500:2000,'ytick', [], 'yticklabel', [],...
            'ylim', [0.5 n_unit+0.5], 'ydir', 'reverse', 'ticklength', [0.025 0.01], 'XTickLabelRotation', 90);
        if ifp > 1
            set(ha_release(ifp), 'xticklabel', []);
        elseif jport==1
            xlabel('From CentIn / CentOut / Reward (ms)');
        end
        if jport > 1
            yticklabels([]);
            ylabel([]);
        end
        ReleasePSTH_thisFP = PSTH_ConNorm{ifp, jport}(:, IndCentOut{ifp, jport});
        imagesc(tCentOutPSTHs{ifp, jport}, 1:n_unit, ReleasePSTH_thisFP, norm_range);
        colormap('Parula')
        yrange = [0.5 n_unit+0.5];
        line([0 0], yrange, 'color', 'w', 'linestyle', ':', 'linewidth', 1.5);
        %         plotshaded([get(gca, 'xlim')], [nsig+.5 nsig+.5; n_unit+0.5 n_unit+0.5], 'w', 0.75);
        fill([get(gca, 'xlim') flip(get(gca, 'xlim'))], [nsig+.5 nsig+.5 n_unit+0.5 n_unit+0.5], 'w', 'FaceColor', 'w', 'EdgeColor', 'none', 'FaceAlpha', .75);
    end
    xlevel_now = xlevel_now + Width + space;

    % Reward
    ylevel_start = 2;
    map_height = 3;
    tRange = [-500 tRewardPSTHs{1}(end)];
    Width  = size_factor*diff(tRange)/2000;
    for ifp =1:nFP
        ha_reward(ifp) = axes('unit', 'centimeters', 'position',...
            [xlevel_now ylevel_start+(map_height+0.5)*(ifp-1)  Width map_height],...
            'nextplot', 'add', 'xlim', tRange, 'xtick', -500:500:2000,'ytick', [], 'yticklabel', [],...
            'ylim', [0.5 n_unit+0.5], 'ydir', 'reverse', 'ticklength', [0.025 0.01], 'XTickLabelRotation', 90);
        if ifp >1
            set(ha_reward(ifp), 'xticklabel', [])
        else
%             xlabel('Reward (ms)')
        end
        if jport > 1
            yticklabels([]);
            ylabel([]);
        end
        RewardPSTH_thisFP = PSTH_ConNorm{ifp, jport}(:, IndReward{ifp, jport});
        imagesc(tRewardPSTHs{ifp, jport}, 1:n_unit, RewardPSTH_thisFP, norm_range);
        colormap('Parula')
        yrange = [0.5 n_unit+0.5];
        line([0 0], yrange, 'color', 'w', 'linestyle', ':', 'linewidth', 1.5);
        fill([get(gca, 'xlim') flip(get(gca, 'xlim'))], [nsig+.5 nsig+.5 n_unit+0.5 n_unit+0.5], 'w', 'FaceColor', 'w', 'EdgeColor', 'none', 'FaceAlpha', .75);
    end

    xlevel_now = xlevel_now + Width + 2*space;
end
% ylevel_now = ylevel_start+(map_height+0.5)*(nFP-1)+map_height+0.5;
xloc = get(gca, 'position');
xboundbar= sum(xloc([1, 3]))+space;
% add color bar
hcbar = colorbar('location', 'eastoutside');
set(hcbar, 'units', 'centimeters', 'position', [xboundbar, ylevel_start, 0.25 3], 'TickDirection', 'out', 'ticklength', 0.025)
hcbar.Label.String = 'normalized';
hcbar.Label.FontSize = 9;

%% B. Plot z scores
% Press
xlevel_start = xboundbar+2.5;
xlevel_now = xlevel_start;
ylevel_start = 2;
map_height = 3;
for jport = 1:nPort
    tRange = [tCentInPSTHs{1}(1) 2000];
    Width  = size_factor*diff(tRange)/2000;
    for ifp =1:nFP
        ha_centin(ifp) = axes('unit', 'centimeters', 'position',...
            [xlevel_now ylevel_start+(map_height+0.5)*(ifp-1)  Width map_height], 'nextplot', 'add',...
            'xlim', tRange, 'xtick', -3500:500:2000,...
            'ylim', [0.5 n_unit+0.5], 'ydir', 'reverse', 'ticklength', [0.025 0.01], 'XTickLabelRotation', 90);
        colormap(ha_centin(ifp), mycolormap);
        if ifp >1
            set(ha_centin(ifp), 'xticklabel', [])
        else
            ylabel('Units')
        end
        if jport > 1
            yticklabels([]);
            ylabel([]);
        end
        imagesc(tCentInPSTHs{ifp, jport}, 1:n_unit, CentInPSTHZs{ifp, jport}, 'AlphaData', ~isnan(CentInPSTHZs{ifp, jport})); clim(zrange);
        yrange = [0.5 n_unit+0.5];
        line([0 0], yrange, 'color', press_col, 'linestyle', ':', 'linewidth', 1.5);
        line([0 0]+MixedFPs(ifp), yrange, 'color', 'm', 'linestyle', ':', 'linewidth', 1.5)
        fill([get(gca, 'xlim') flip(get(gca, 'xlim'))], [nsig+.5 nsig+.5 n_unit+0.5 n_unit+0.5], 'w', 'FaceColor', 'w', 'EdgeColor', 'none', 'FaceAlpha', .75);
        if jport==1
            title(sprintf('FP=%2.0dms, Ntrials=%2.0d', MixedFPs(ifp), Pop.Trials(ifp, jport)));
        else
            title(sprintf('Ntrials=%2.0d', Pop.Trials(ifp, jport)));
        end
    end

    ylevel_now = ylevel_start+(map_height+0.5)*(nFP-1)+map_height+0.5;
    uicontrol('Style','text','Units','centimeters','Position',[xlevel_start-1 ylevel_now  6 0.7],...
        'string', ['B. z-scored activity [' num2str(zrange(1)) '-' num2str(zrange(2)) ']'], ...
        'FontName','Dejavu Sans',  'fontweight', 'bold','fontsize', 9,'BackgroundColor',[1 1 1],...
        'HorizontalAlignment','Left');
    xlevel_now = xlevel_now + Width+space;

    % Release
    ylevel_start = 2;
    map_height = 3;
    tRange = [-600 tCentOutPSTHs{1}(end)];
    Width  = size_factor*diff(tRange)/2000;
    for ifp =1:nFP
        ha_release(ifp) = axes('unit', 'centimeters', 'position',...
            [xlevel_now ylevel_start+(map_height+0.5)*(ifp-1)  Width map_height],...
            'nextplot', 'add', 'xlim', tRange, 'xtick', -500:500:2000,'ytick', [], 'yticklabel', [],...
            'ylim', [0.5 n_unit+0.5], 'ydir', 'reverse', 'ticklength', [0.025 0.01], 'XTickLabelRotation', 90);
        colormap(ha_release(ifp), mycolormap);
        if ifp >1
            set(ha_release(ifp), 'xticklabel', [])
        elseif jport==1
            xlabel('From CentIn / CentOut / Reward (ms)');
        end
        if jport > 1
            yticklabels([]);
            ylabel([]);
        end
        imagesc(tCentOutPSTHs{ifp, jport} , 1:n_unit,  CentOutPSTHZs{ifp, jport}, 'AlphaData', ~isnan(CentOutPSTHZs{ifp, jport})); clim(zrange);
        yrange = [0.5 n_unit+0.5];
        line([0 0], yrange, 'color', 'k', 'linestyle', ':', 'linewidth', 1.5);
        fill([get(gca, 'xlim') flip(get(gca, 'xlim'))], [nsig+.5 nsig+.5 n_unit+0.5 n_unit+0.5], 'w', 'FaceColor', 'w', 'EdgeColor', 'none', 'FaceAlpha', .75);
    end
    xlevel_now = xlevel_now + Width+space;

    % Reward
    ylevel_start = 2;
    map_height = 3;
    tRange = [-1000 tRewardPSTHs{1}(end)];
    Width  = size_factor*diff(tRange)/2000;
    for ifp =1:nFP
        ha_reward(ifp) = axes('unit', 'centimeters', 'position',...
            [xlevel_now ylevel_start+(map_height+0.5)*(ifp-1)  Width map_height],...
            'nextplot', 'add', 'xlim', tRange, 'xtick', -500:500:2000,'ytick', [], 'yticklabel', [],...
            'ylim', [0.5 n_unit+0.5], 'ydir', 'reverse', 'ticklength', [0.025 0.01], 'XTickLabelRotation', 90);
        colormap(ha_reward(ifp), mycolormap);
        if ifp >1
            set(ha_reward(ifp), 'xticklabel', [])
        else
%             xlabel('Reward (ms)')
        end
        if jport > 1
            yticklabels([]);
            ylabel([]);
        end
        imagesc(tRewardPSTHs{ifp, jport}, 1:n_unit, RewardPSTHZs{ifp, jport}, 'AlphaData', ~isnan(RewardPSTHZs{ifp, jport})); clim(zrange);
        yrange = [0.5 n_unit+0.5];
        line([0 0], yrange, 'color', 'w', 'linestyle', ':', 'linewidth', 1.5);
        fill([get(gca, 'xlim') flip(get(gca, 'xlim'))], [nsig+.5 nsig+.5 n_unit+0.5 n_unit+0.5], 'w', 'FaceColor', 'w', 'EdgeColor', 'none', 'FaceAlpha', .75);
    end

    xlevel_now = xlevel_now + Width + 2 * space;
end
xloc = get(gca, 'position');
xboundbar= sum(xloc([1, 3]))+0.1;

% add color bar
hcbar = colorbar('location', 'eastoutside');
set(hcbar, 'units', 'centimeters', 'position', [xboundbar, ylevel_start, 0.25 3], 'TickDirection', 'out', 'ticklength', 0.025)
hcbar.Label.String = 'z score';
hcbar.Label.FontSize = 9;

ylevel_now = ylevel_start+(map_height+0.5)*(nFP-1)+map_height+2;

%% add unit table
% % 
% % % write the info to a csv table.
% % % number of sorted units in this experiment-this is also the number of
% % % rows
% % n_unit = size(Pop.Units, 1);
% % % anm name
% % Name = repmat(Pop.Name, n_unit, 1);
% % % session
% % Session = repmat(Pop.Date, n_unit, 1);
% % 
% % % Sorting order
% % Chs              = zeros(n_unit, 1);
% % Ch_Units         = zeros(n_unit, 1);
% % Unit_Quality_Num = zeros(n_unit, 1);
% % Unit_Quality     = cell(n_unit, 1);
% % % IndUnmodulated
% % SignificantMod = ones(n_unit, nPort);
% % for j = 1:nPort
% %     SignificantMod(Pop.IndUnmodulated{j}, j) = 0;
% % end
% % Unit_Sorted = cell2mat(Pop.IndSort);
% % 
% % for i = 1:size(Unit_Sorted, 1)
% %     iCode = Pop.Units(Pop.IndSort(i), :);
% %     Chs(i)      = iCode(1);
% %     Ch_Units(i) = iCode(2);
% %     if iCode(3) == 1
% %         Unit_Quality{i} = 's';
% %     else
% %         Unit_Quality{i} = 'm';
% %     end
% %     Unit_Quality_Num(i) = iCode(3);
% % end
% % 
% % tab = table(Name, Session, Unit_Sorted, Chs, Ch_Units, Unit_Quality_Num, SignificantMod);
% % % aGoodName = ['PopOut', Pop.Name, '_' Pop.Session '.csv'];
% % % writetable(tab, aGoodName)
% % % % open this table
% % % try
% % %     winopen(aGoodName)
% % % end;
% % 
% % % Insert table
% % table_data = [Unit_Sorted, Chs, Ch_Units, Unit_Quality_Num SignificantMod];
% % htable = uitable(hf, 'unit', 'centimeters','Data', table_data,...
% %     'ColumnName', {'Unit|Sorted', 'Channel', 'Unit|Channel',  'Quality', 'StatSignificant'},...,
% %     'Position', [1 ylevel_now 8 3], 'ColumnWidth',{50, 50, 50, 50, 100});
% % 

%% Add information
uicontrol('Style','text','Units','centimeters','Position', [1 ylevel_now 6 2],...
    'string', Pop.Name+' | '+strrep(Pop.Session, '_', '-'), ...
    'FontName','Dejavu Sans', 'fontweight', 'bold','fontsize', 12,'BackgroundColor','w')
% Re-adjust figure height

fig_height = ylevel_now + 3+0.5;

figsize = get(hf, 'Position');
figsize(4) = fig_height;
set(hf, 'Position', figsize)

%% Save this figure
thisFolder = fullfile(pwd, 'Fig');
if ~exist(thisFolder, 'dir')
    mkdir(thisFolder)
end
tosavename= fullfile(thisFolder, 'PopulationActivity_'+Pop.Name+'_'+strrep(Pop.Session, '_', ''));
disp('########## making figure ########## ')
tic
print (gcf,'-dpdf', tosavename)
print (gcf,'-dpng', tosavename, '-r300')
print (gcf,'-depsc2', tosavename);
toc
% % try
% %     thisFolder = fullfile(findonedrive, '00_Work' , '03_Projects', '05_Physiology', 'Data', 'PopulationPSTH', Pop.Name);
% %     if ~exist(thisFolder, 'dir')
% %         mkdir(thisFolder)
% %     end
% %     disp('##########  copying figure ########## ')
% %     tic
% %     copyfile([tosavename '.png'], thisFolder)
% %     copyfile([tosavename '.eps'], thisFolder)
% %     toc
% % end