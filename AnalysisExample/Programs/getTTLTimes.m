function ttls = getTTLTimes(mainpath)
%function ttls = getTTLTimes(mainpath)
% given mainpath, get TTL timings


LoadVision2;
spath = slashappend(mainpath);
sfile = LoadVisionFiles([spath '*.spikes']);

rawttl = sfile.getTTLTimes;
ttls = double(rawttl) / 20; % double in ms
sfile.close;
