function [maxval, inds] = maxmd(data)
% function [maxval, inds] = maxmd(data)
% multi-dimensional max

ds = size(data);
nd = length(ds);
[maxval, lind] = max(data(:));
inds = zeros(1,nd);

outarg = '[';
for i = 1:nd
    outarg = [outarg sprintf('inds(%d),', i)];
end
outarg = [outarg ']'];
eval([outarg ' = ind2sub(ds, lind);']);