function [initp, bounds] = AuditorySetInits(inputdata)
% setting the initial values and boundary for the auditory spots.

medy = median(inputdata.data(:, 1));
medx = median(inputdata.data(:, 2));
medt = median(inputdata.data(:, 3));

stdy = std(inputdata.data(:, 1));
stdx = std(inputdata.data(:, 2));
stdt = std(inputdata.data(:, 3));

initp = [medy, medx, stdy, stdx, 0, medt, stdt, stdt, 0.1, 0.1];
%initp = [2, 3, stdy, stdx, 0, 0.01, 0.01, 0.01, 0.1, 0.1];


% set step sizes
stepsize = initp / 100;
stepsize(5) = pi / 180;
stepsize(10) = 0.1 / 100;

% if there is only one dimension for Y, move to 1D space mode.
if stdy == 0
    initp(3) = 1; % any value is fine.
    stepsize([1, 3, 5]) = 0; % turn off Y-value changes
end

% set bounds
% [y, x, sy, sx, theta, t, stl, str, 
bound1 = [             1,              1, 0.5, 0.5, -2*pi, 0.003, 0.0001, 0.0001,  0, 0];
bound2 = [inputdata.maxy, inputdata.maxx,  10,   9,  2*pi,  0.15,   0.1,   0.1, 10, 1];
bounds = [1:length(initp);stepsize;bound1;bound2]';
