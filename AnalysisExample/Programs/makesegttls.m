function makesegttls(datafolder)
% operate on a datafolder, and make a segment ttl file.

datafolder = slashappend(datafolder);
ttlTimes = getTTLTimes(datafolder);
sl = load([datafolder 'segmentlengths']);
ssep = [0 sl.segmentseparations];

segttls = {};
nstim = length(sl.segmentseparations);
for i = 1:nstim
    segttl = ttlTimes(ttlTimes > ssep(i) & ttlTimes < ssep(i+1));
    %segttl = TTLGlitchRemoval(segttl, sf{i});
    segttls{i} = segttl;
end
save([datafolder, 'segttls.mat'], 'segttls');
