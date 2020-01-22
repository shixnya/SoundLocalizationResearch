function ndata = normalizeRows(data, mode)
% function normalizeRows(data, mode)
% mode can be 'max', 'norm', 'mean'


if nargin < 2
    mode = 'norm'; % changing the default mode
end

dsize = size(data);

if strcmp(mode, 'max')
    ndata = data ./ repmat(max(data, [], 2), 1, dsize(2));
elseif strcmp(mode, 'norm')
    ndata = zeros(dsize);
    for i = 1:dsize(1)
        ndata(i,:) = data(i,:) / norm(data(i,:));
    end
elseif strcmp(mode, 'mean')
    ndata = data ./ repmat(mean(data, 2), 1, dsize(2));
elseif strcmp(mode, 'std')
    ndata = data ./ repmat(std(data, 0, 2), 1, dsize(2));
else
    mode
    error('Mode unrecognized')
end
