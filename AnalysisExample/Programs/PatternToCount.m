function FR = PatternToCount(pat, offset, maxtime)
% function FR = PatternToCount(pat, offset, maxtime)
% fixed critical bug on 12/22/15. Any analysis before this date is wrong.

if nargin < 2
    offset = 0;
end
if nargin < 3
    maxtime = inf;
end

FR = zeros(size(pat));
if isinf(maxtime)
    for j = 1:length(pat(:))
        FR(j) = nnz(pat{j} >= offset);
    end
else
    for j = 1:length(pat(:))
        FR(j) = nnz(pat{j} >= offset & pat{j} < maxtime);
    end
end