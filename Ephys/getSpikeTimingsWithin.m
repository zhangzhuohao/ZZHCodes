function t_spk_out = getSpikeTimingsWithin(t_spk, t_event, t_window)
%getSpikeTimingsWithin Extract timing of spikes within a time-window around a certain event.
%   Input: 
%       t_spk: 1*n array, timing of spikes during recording period
%       t_event: n*1 array, timing of events
%       t_window: 1*2 or n*2 array (n is equal to the length of t_event), time-window for extracting spike timings

if isempty(t_event)
    t_spk_out = [];
    return
end

n_event = length(t_event);
sz_w  = size(t_window);
if sz_w(1)~=2 && sz_w(2)~=2
    error('Time-window must be a 1*2 or n*2 array (n is equal to the length of t_event).')
elseif sz_w(1)==2 && sz_w(2)~=2
    t_window = t_window';
end
sz_w = size(t_window);

assert(ismember(sz_w(1), [1 n_event]), 'The number of t_window rows must be 1 or equal to the length of t_event.')
assert(all(t_window(:,2)>t_window(:,1)), 'Each row of t_window must be a 2-element vector of increasing numeric values.')

t_spk_out = cell(n_event, 1);
if sz_w(1)==1
    t_window = repmat(t_window, n_event, 1);
end

for i = 1:n_event
    t_spk_out{i} = t_spk(t_spk>=t_event(i)+t_window(i,1) & t_spk<t_event(i)+t_window(i,2)) - t_event(i);
end

end
