function psths = fastttldivide(sptrain_or_asdf, ttl, endtime, offset)
% function psths = fastttldivide(sptrain_or_asdf, ttl, endtime, offset)
% this function does something like
% psth{elem} = sp(sp >= ttl(i) & sp < ttl(i) + endtime);
% but using histcounts (so that it does not double count spikes)
% and faster.
% if the first element is asdf, it returns {nTTL, nNeu} cell array.
% otherwise, it returns {nTTL} cell array
% endtime is in ms

if nargin < 4
    offset = 0;
end

if iscell(sptrain_or_asdf)
    % asdf mode
    nNeu = sptrain_or_asdf{end}(1);
    psths = cell(length(ttl), nNeu);
    for n = 1:nNeu
        psths(:, n) = fastttldivide_sub(sptrain_or_asdf{n}, ttl, endtime, offset);
    end
else
    % spike train mode
    psths = fastttldivide_sub(sptrain_or_asdf, ttl, endtime, offset);
end



function psths = fastttldivide_sub(sptrain, ttl, endtime, offset)
hc = histcounts(sptrain, [0; ttl; max(sptrain(end), ttl(end))]);
cind = hc(1);
nttl = length(ttl);
psths = cell(nttl, 1);

for i = 1:nttl
    if hc(i+1) == 0
        continue
    end
    psths{i} = sptrain(cind + 1: cind + hc(i+1)) - ttl(i);
    if psths{i}(end) >= endtime
        psths{i} = psths{i}(psths{i} < endtime);
    end
    
    if ~isempty(psths{i}) & psths{i}(1) < offset
        psths{i} = psths{i}(psths{i} >= offset);
    end
    cind = cind + hc(i+1);
end
