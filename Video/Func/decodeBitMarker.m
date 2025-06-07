function [trial_id, ind_header] = decodeBitMarker(tOn, tOff, varargin)

% Default params
% led bit code information
bit_len = 100; % ms
bit_num = 9; % number of bits used
% header information
header_on_dur  = 200; % ms
header_off_dur = 100; % ms
% sample rate
fs = 50; % ms
if nargin>2
    for i=1:2:size(varargin,2)
        switch varargin{i}
            case 'bit_len'
                bit_len = varargin{i+1};
            case 'bit_num'
                bit_num = varargin{i+1};
            case 'header_on_dur'
                header_on_dur = varargin{i+1};
            case 'header_off_dur'
                header_off_dur = varargin{i+1};
            case 'fs'
                fs = varargin{i+1};
            otherwise
                errordlg('unknown argument')
        end
    end
end

% led bit code information
bit_dur = bit_len*bit_num;

% header information
header_dur  = header_on_dur + header_off_dur;
header_code = [repmat('1', 1, header_on_dur/bit_len) repmat('0', 1, header_off_dur/bit_len)];
header_num  = length(header_code);

% entire signal information
signal_num = header_num + bit_num;
signal_dur = header_dur + bit_dur;

% exposure time
exposure_time = round(1000 / fs);

% check signal from the first LED-on
trial_id   = nan(length(tOn), 1); % set a larger initial container
ind_header = (1:length(tOn))';
i = 1;
while i < length(tOn)
    % find time of LED-on within signal duration
    t_on_follow = tOn(i:end) - tOn(i);
    t_on_follow = t_on_follow(t_on_follow<=signal_dur);

    % check if there is a bit sequence following, otherwise go to next LED-on
    if length(t_on_follow)<2
        i = i+1;
        continue;
    end

    % find time of LED-off within signal duration
    t_off_follow = tOff(i:end) - tOn(i);
    t_off_follow = t_off_follow(t_off_follow<=header_dur+bit_dur+1000); % include extra 1000 ms, in case of timestamp shift

    % transfer time points to signal sequence
    entire_signal = zeros(1, signal_dur);
    t_signal   = 1:signal_dur;
    for j = 1:length(t_on_follow)
        entire_signal(t_signal>=t_on_follow(j) & t_signal<t_off_follow(j)+exposure_time) = 1;
    end

    % transfer signal sequence to bit code
    entire_code = zeros(1, signal_num);
    for j = 1:signal_num
        j_loc = (j-1)*bit_len + (1:bit_len);
        entire_code(j) = round(mean(entire_signal(j_loc)));
    end

    % check if there is a header code in front
    if any(num2str(entire_code(1:header_num), '%d') ~= header_code)
        i = i+1;
        continue;
    end

    % turn bit code to trial number
    bit_code = entire_code(header_num+1:end);
    i_trial  = bin2dec(num2str(bit_code, '%d'));

    trial_id(i) = i_trial;
    % go to next header
    i = i + length(t_on_follow);
end
ind_header(isnan(trial_id)) = [];
trial_id(isnan(trial_id))   = [];

end