function [segfr, segdisper] = SegmentFR(asdf, ttls, starttime, endtime)
% calculate firing rate and dispersion of a segment between start time and end time.

nNeu = asdf{end}(1);

% characterize the spontaneous firing rate
%endtime = 2000;
%offset = 1000;
offset = starttime;
duration = (endtime - offset) / 1000; % making it second

%tic
frcells = fastttldivide(asdf, ttls, endtime, offset);
%toc

%segfr = zeros(nNeu, 1);
%nseg = size(frcells, 1);

counts = cellfun(@length, frcells);
meanfr = mean(counts)';
segdisper = var(counts)' ./ meanfr;
segfr = meanfr / duration;
%for i = 1:nNeu
%    segfr(i) = length(cat(2, frcells{:, i})) / (nseg * duration);
%end
