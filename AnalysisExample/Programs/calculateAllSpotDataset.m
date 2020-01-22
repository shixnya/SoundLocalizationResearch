function res = calculateAllSpotDataset(patterns, sf, mode)
%

if nargin < 3
    mode = 1;
end

nNeu = size(patterns.fullpat, 5);
nCol = size(patterns.fullpat, 1);


if nNeu == 0
    res = [];
    return
end
res(nNeu, nCol).params = nan;
res(nNeu, nCol).nllf = nan;
res(nNeu, nCol).cov = nan;
res(nNeu, nCol).nSpike = nan;


for i = 1:nNeu
    res(i,:) = SpotFit_onecell(patterns.fullpat(:,:,:,:,i), sf, mode);
end