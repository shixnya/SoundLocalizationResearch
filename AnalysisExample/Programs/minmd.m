function [minval, inds] = minmd(data)
% function [minval, inds] = minmd(data)
% multi-dim minimum finding.
[maxval, inds] = maxmd(-data);
minval = - maxval;