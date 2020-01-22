function eventlist = patternToEventList(pattern)
% function eventlist = patternToEventList(pattern)
% convert a pattern to an event list

% determine the number of dimensions of the pattern
% assume that the pattern is a cell array of spikes

s = size(pattern);
ndim = length(s);
npat = prod(s);
try
    nevent = length(cat(2, pattern{:}));
catch
    nevent = length(cat(1, pattern{:})); % ad hoc change on 8/23. Can be removed after spot analysis for kiki is redone.
end

indholder = cell(ndim, 1);

eventlist = zeros(nevent, ndim + 1); % last one is the time


cloc = 0;
for i = 1:npat
    [indholder{:}] = ind2sub(s, i);
    nevent_each = length(pattern{i});
    eventlist(cloc+1:cloc+nevent_each, 1:end-1) = ...
        repmat([indholder{:}], nevent_each, 1);
    eventlist(cloc+1:cloc+nevent_each, end) = pattern{i} / 1000; % converting ms to s
    cloc = cloc + nevent_each;
end

