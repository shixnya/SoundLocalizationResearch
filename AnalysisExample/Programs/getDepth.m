function [depth, cdepth, cx, cy, strct] = getDepth(datafolder, electrodes)
% function [depth, cdepth, cx, cy] = getDepth(datafolder)
% get neurons' heights using SurfaceContour
%
% optionally, it gives the height from the top neurons if SurfaceContour
% does not exist in the folder, giving a warning.
% 
% if the second argument is 1, it gives the location of the electrodes.

if nargin < 2
    electrodes = 0;
end

datafolder = slashappend(datafolder);


if electrodes
    % define x and y for the electrodes
    load([datafolder 'basicinfo'])
    emap = loademap(arrayID);
    x = emap(:,1);
    y = emap(:,2);
else
    load([datafolder 'xy']);
end
nNeu = length(x);
depth = zeros(nNeu, 1);

cname = [datafolder 'SurfaceContour.mat'];

if exist(cname, 'file')
    load(cname)
    for i = 1:nNeu
        [~, minind] = min(abs(cont(:,1) + y(i)));
        depth(i) = cont(minind, 2) - x(i);
    end
else % in case it is not calculated
    %warning('Surface definition file %s does not exist.', cname);
    %depth = max(x) - x;
    depth(:) = nan;
end


if nargout > 1 
    % calculate corredted depth
    insang = getInsertionAngle(datafolder);
    cdepth = depth * cos(insang);
    cx = x * cos(insang) + y * sin(insang);
    cy = y * cos(insang) - x * sin(insang);
    strct.depth = depth;
    strct.cdepth = cdepth;
    strct.x = x;
    strct.y = y;
    strct.cx = x;
    strct.cy = y;
end