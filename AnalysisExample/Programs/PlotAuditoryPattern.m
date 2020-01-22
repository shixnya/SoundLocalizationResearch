function PlotAuditoryPattern(patterns, windows, plothemi, show)

%windows = [5 25; 25 60; 60 200; 5 200; 0 0; 0 200];
%pat = pattern2d;

if nargin < 3
    plothemi = 1;
end
if nargin < 4
    show = [];
end

if plothemi
    nrow = 5;
else
    nrow = 4;
end
rasterrange = windows(end, :);

clf;
subplot(nrow, 1, 1:2);
SuperFastPatternRaster(patterns, rasterrange, show);

wsize = diff(windows(1:end-1, :), 1, 2);
windows(wsize == 0, :) = [];

nwin = size(windows, 1) - 1;
kp = [];
for i = 1:nwin
    spatpat = PatternToCount(patterns, windows(i, 1), windows(i, 2));
    spatpat(squeeze(~show)) = nan; % delete non show info
    meanpat = mean(spatpat, 3, 'omitnan');
    subplot(nrow, nwin, nwin * 2 + i);
    plot(meanpat');
    title(sprintf('IW%d', i));

    subplot(nrow, nwin, nwin * 3 + i);
    imagesc(meanpat);
    set(gca, 'ydir', 'normal');
    set(gca, 'xtick', [0.5, 9.5, 17.5]);
    set(gca, 'xticklabel', [-144, 0, 144]);
    set(gca, 'ytick', [1, 3, 5]);
    set(gca, 'yticklabel', [0, 40, 80]);
    axis image
    colorbar
    
    if plothemi
        kp(i) = subplot(nrow, nwin, nwin * 4 + i);
        PlotFullHemisphere(meanpat);
    end
end

colormap viridis

Link = linkprop(kp, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);
