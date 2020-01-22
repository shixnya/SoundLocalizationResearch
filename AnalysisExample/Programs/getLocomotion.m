function [asdf_loco, locoinfo_all] = getLocomotion(datafolder)
if nargin < 1
    datafolder = pwd;
end
datafolder = slashappend(datafolder);
locofile = [datafolder, 'EncoderLocomotion.mat'];

if exist(locofile, 'file')
    load(locofile);
    return
end


rpath = readpath(datafolder);


removeind = logical(zeros(1, length(rpath)));
dataind = 0;
for i = 1:length(rpath)
    ri = rpath{i};
    removeind(i) = isempty(ri) | contains(ri, 'merged') | contains(ri, 'processed');
end
rpath(removeind) = [];
datasep = strsplit(strrep(rpath{end}, '/', ''), ',');
folders = multicat(strcat(rpath{1:end-1}), datasep);
locoinfos = cellfun(@(x) load([slashappend(x) 'locoinfo.mat']), folders);
locoinfo_all = cat(2, locoinfos.locoinfo);
% make ASDF based on this locoinfo.

asdf_loco{1} = find(diff(sum(locoinfo_all))) / 20;
asdf_loco{2} = 1;
asdf_loco{3} = [1, ceil(size(locoinfo_all, 2) / 20)];
%find(diff(locoinfo_all(1,:)));
%find(diff(locoinfo_all(2,:)));
save('-v7.3', locofile, 'asdf_loco', 'locoinfo_all');