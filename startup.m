warning('off', 'MATLAB:print:ContentTypeImageSuggested');
warning('off', 'MATLAB:stats:pca:ColRankDefX');
try
    set_matlab_default
catch
    disp('You do not have "set_matlab_default"');
end

% Call Psychtoolbox-3 specific startup function:
if exist('PsychStartup'), PsychStartup; end;

