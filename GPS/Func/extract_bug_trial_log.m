function bug_table = extract_bug_trial_log(file_path)
%EXTRACT_BUG_TRIAL_LOG Read bug types and trial IDs from a text log.
%   BUG_TABLE = EXTRACT_BUG_TRIAL_LOG(FILE_PATH) reads nonempty lines in
%   the following format:
%
%       bug_type: trial_id_1, trial_id_2, ...
%
%   The returned table contains one row per trial, with variables
%   "bug_type" (string) and "trial_id" (double). Different lines may use
%   different bug-type names. Commas, semicolons, and whitespace are
%   accepted as separators between trial IDs.

    if isstring(file_path)
        if ~isscalar(file_path)
            error('extract_bug_trial_log:InvalidPath', ...
                'file_path must be a character vector or string scalar.');
        end
        file_path = char(file_path);
    elseif ~ischar(file_path)
        error('extract_bug_trial_log:InvalidPath', ...
            'file_path must be a character vector or string scalar.');
    end

    fid = fopen(file_path, 'rt');
    if fid == -1
        error('extract_bug_trial_log:FileOpenFailed', ...
            'Could not open file: %s', file_path);
    end
    file_cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    file_data = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = file_data{1};

    bug_types = strings(0, 1);
    trial_ids = zeros(0, 1);

    for line_number = 1:numel(lines)
        line_text = strtrim(strrep(lines{line_number}, char(13), ''));

        % Remove a UTF-8 byte-order mark if it appears at the beginning.
        if line_number == 1 && ~isempty(line_text) && line_text(1) == char(65279)
            line_text = line_text(2:end);
        end

        if isempty(line_text)
            continue;
        end

        colon_position = find(line_text == ':', 1, 'first');
        if isempty(colon_position)
            error('extract_bug_trial_log:MissingColon', ...
                'Line %d does not contain a colon: %s', ...
                line_number, line_text);
        end

        bug_type = strtrim(line_text(1:colon_position - 1));
        id_text = strtrim(line_text(colon_position + 1:end));
        if isempty(bug_type)
            error('extract_bug_trial_log:MissingBugType', ...
                'Line %d has an empty bug type.', line_number);
        end
        if isempty(id_text)
            error('extract_bug_trial_log:MissingTrialID', ...
                'Line %d has no trial IDs.', line_number);
        end

        % Accept English/Chinese commas and semicolons, plus whitespace.
        id_parts = regexp(id_text, '[,;，；\s]+', 'split');
        id_parts = id_parts(~cellfun('isempty', id_parts));
        line_trial_ids = str2double(id_parts);

        if any(isnan(line_trial_ids))
            bad_parts = id_parts(isnan(line_trial_ids));
            error('extract_bug_trial_log:InvalidTrialID', ...
                'Line %d contains a nonnumeric trial ID: %s', ...
                line_number, strjoin(bad_parts, ', '));
        end

        line_trial_ids = line_trial_ids(:);
        bug_types = [bug_types; repmat(string(bug_type), numel(line_trial_ids), 1)]; %#ok<AGROW>
        trial_ids = [trial_ids; line_trial_ids]; %#ok<AGROW>
    end

    bug_table = table(bug_types, trial_ids, ...
        'VariableNames', {'bug_type', 'trial_id'});
end
