function [IndOut, FigAlign] = findseqmatchrev(seqmom, seqson, man, toprint, toprintname, threshold)
% Jianing Yu 4/19/2021
% find matching index between seqmom and seqson
% output indout is the index of each seqson's element in seqmom;
% NaN means no reliable match is found.
% time unit should be ms
% length of seqmom should be greater than that of seqson

% revised to facilitate fast alignment

if nargin<6
    threshold = 5;
    if nargin<5
        toprintname = 'Alignment';
        if nargin<4
            toprint = 0;
            if nargin<3
                man = 1;
            end
        end
    end
end

tic
IndOut = zeros(1, length(seqson));
maxcorr = zeros(1, length(seqson));
basecorr = zeros(1, length(seqson));

disp('#############################')
disp('########## Please wait #########')
disp('#############################')

lastonechecked = 0;

for i = 1:length(seqson)
    if i > 1
        if IndOut(i-1)~=0
            lastonechecked = 1;
        else
            lastonechecked = 0;
        end
    end

    if ~rem(i, 100)
        sprintf('check seqson %2.0f of %2.0f', i, length(seqson))
    end

    corr_momson = zeros(1, length(seqmom));

    if ~lastonechecked
        for k = 1:length(seqmom)
            corr_momson(k) = toalign(seqmom-seqmom(k), seqson-seqson(i));
        end
    else
        for k = IndOut(i-1):length(seqmom)
            corr_momson(k) = toalign(seqmom-seqmom(k), seqson-seqson(i));
        end
    end

    [maxcorr(i), indmax] =  max(corr_momson);
    basecorr(i) = threshold*std(corr_momson(setdiff(1:length(corr_momson), indmax)));
    if maxcorr(i) > basecorr(i) % max correction has to be much larger than the baseline
        IndOut(i) = indmax;
    end
end

% in addition, indout must increase monotomically

FigAlign = figure(38); clf;
set(gcf, 'name', 'Alignment', 'units', 'centimeters', 'position', [5 1 15 25]);
ha1 = subplot(3, 1, 1);
indorg = find(IndOut>0, 1, 'first');
plot(seqmom - seqmom(IndOut(indorg)), 4, 'ro');
hold on
plot(seqson - seqson(indorg), 3.5, 'k*')
text(0, 4.2, 'seqmom', 'color', 'r')
text(0, 3.7, 'seqson', 'color', 'k')
set(ha1, 'ylim', [2.5 5])

toc
ha2 = subplot(3, 1, 2);
plot(IndOut, maxcorr, 'linewidth', 1);
hold on
plot(IndOut, basecorr, 'r-', 'linewidth', .5)
xlabel('Indout')
ylabel('Corr')

ha3 = subplot(3, 1, 3);
set(gca, 'nextplot', 'add')
xlabel('Corr')
ylabel('Count')
hhist1 = histogram(maxcorr, 100); set(hhist1, 'facecolor', 'b');
hhist2 = histogram(basecorr, 100); set(hhist2, 'facecolor', 'r');
text(min(get(ha3, 'xlim')), max(get(ha3, 'ylim')), 'Select a threshold!', 'fontsize', 15)

if man
    % set a threshold to get rid of outliers.
    disp('Select one point to define min corr, end seleciton by right click');
    [x_thrh, ~] = getpts(gcf);
    corr_min = min(x_thrh); % this is the threshold to extract LED_on times
else
    corr_min = mean(basecorr) + 5*std(basecorr);
end

line([corr_min corr_min], get(gca, 'ylim'), 'color', 'm', 'linestyle', ':', 'linewidth', 5)

ind_outlier = find(maxcorr<=corr_min);
if ~isempty(ind_outlier)
    plot(ha2, IndOut(ind_outlier), maxcorr(ind_outlier), 'mo', 'markersize', 10)
    IndOut(ind_outlier) = NaN;
    plot(ha2, IndOut, maxcorr, 'linewidth', 1, 'color', [.1 .6 .1]);
end

if toprint
    print (gcf, '-dpng', toprintname);
end


function seqcorr = toalign(seq1, seq2)
% seq 1 is a subset of seq2
% based on correlation analysis

if size(seq1, 1)>size(seq1, 2)
    seq1 = seq1';
end
if size(seq2, 1)>size(seq2, 2)
    seq2 = seq2';
end

tmax = max([seq1, seq2]);
tmin = min([seq1, seq2]);

% tdiv = (prctile([diff(seq1) diff(seq2)], 10));
tdiv = (tmax - tmin) / (25*max([length(seq1), length(seq2)]));
edges = tmin:tdiv:tmax;

% won't assume
nseq1 = histcounts(seq1, edges);
nseq2 = histcounts(seq2, edges);
seqcorr = nseq1*nseq2';
