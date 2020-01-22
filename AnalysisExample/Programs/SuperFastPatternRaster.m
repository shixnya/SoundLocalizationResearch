function SuperFastPatternRaster(pattern2d, timerange, running)
%function SuperFastPatternRaster(pattern2d, timerange)
%
%    pattern2d - {ny, nx, nrep} 2d PSTH pattern with nrep.
%    timerange (1, 1) or (1, 2) - range of time in ms. if length is 1, it means 
%    max and min is 0.
%
% Returns:
%
%
% Description:
%    Put all of the patterns into one big image to make it faster.
%
% Example:
%    
%
% Requires:
%
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/10/2018
% Modified: 

%pattern2d = permute(patterns_2s.fullpat(1, 1, :, :, 1), [2 3 4 1]);
%timerange = 1000;

if length(timerange) == 1
    mintime = 0;
    maxtime = timerange;
else
    mintime = timerange(1);
    maxtime = timerange(2);
end
dur = maxtime - mintime;



% first, plot the running state
if ~isempty(running)
    rs = size(running);
    rp = reshape(permute(running, [4 2 3 1]), rs(2) * rs(4), rs(3));
    rp(end+1, end+1) = 0; % push the edge to use with pcolor
    x = mintime:dur:(mintime + dur * rs(3));
    y = 0:rs(2)*rs(4);
    pc = pcolor(x, y, double(rp));
    set(pc, 'facealpha', 0.2);
    shading flat
    colormap viridis
    hold on
end

% get information from one panel
ns = size(pattern2d); % [ny, nx, nrep];

plx_b = cell(ns(1), ns(2));
ply_b = cell(ns(1), ns(2));
for i = 1:ns(1)
    for j = 1:ns(2)
        currentpat = squeeze(pattern2d(i, j, :));
        [plx_b{i, j}, ply_b{i, j}] = getOnePattern(currentpat, mintime, maxtime);
        plx_b{i, j} = plx_b{i, j} + (j-1) * dur;
        ply_b{i, j} = ply_b{i, j} + (i-1) * ns(3);
    end
end

grids_x = zeros(3, ns(1) + ns(2));
grids_x(3, :) = nan;
grids_y = zeros(3, ns(1) + ns(2));
grids_y(3, :) = nan;
for i = 1:ns(1)
    grids_x(1, i) = 0;
    grids_x(2, i) = dur * ns(2);
    grids_y(1:2, i) = (i-1) * ns(3);
end
for j = 1:ns(2)
    grids_x(1:2, ns(1) + j) = dur * (j - 1);
    grids_y(1, ns(1) + j) = 0;
    grids_y(2, ns(1) + j) = ns(3) * ns(1);
end

plx = cat(2, plx_b{:});
ply = cat(2, ply_b{:});

%figure(510);clf
plot(plx(:), ply(:), 'linewidth', 2.5);
hold on
plot(grids_x(:), grids_y(:), 'k', 'linewidth', 1);
xlim(mintime + [0, dur * ns(2)]);
ylim([0, ns(3) * ns(1)]);
set(gca, 'XTick', [mintime, maxtime]);
set(gca, 'YTick', [0, ns(3)]);





function [plx, ply] = getOnePattern(onepattern, mintime, maxtime, running)
redpattern = cellfun(@(x) limtime(x, mintime, maxtime), onepattern, 'UniformOutput', false);

lens = cellfun('length', redpattern);
inds = [0; cumsum(lens)];
plx = zeros(3, inds(end));
ply = zeros(3, inds(end));
for i = 1:size(onepattern, 1)
    plx(1:2, inds(i)+1:inds(i+1)) = repmat(redpattern{i}, 2, 1);
    plx(3, inds(i)+1:inds(i+1)) = nan;
    ply(1, inds(i)+1:inds(i+1)) = i-1;
    ply(2, inds(i)+1:inds(i+1)) = i;
    ply(3, inds(i)+1:inds(i+1)) = nan;
end

function lim = limtime(list, mintime, maxtime)
lim = list(list > mintime & list < maxtime);

