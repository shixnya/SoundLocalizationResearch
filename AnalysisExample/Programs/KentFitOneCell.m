function res = KentFitOneCell(pattern, mode, ignorelist)
% do one cell fit with Kent distribution
% Binned maximum likelihood fit with quasi-Poisson distribution.
% (means that the bursty neurons are corrected by fano factor.
% mode 1, 2 use front-z coordinates.
% mode 3, 4 are for top-z coordinates.
% mode 9, 10 are new auditory model
% modified on 6/11, fixing the issue with inflated error at zero-bins and
% some others

if nargin < 3
    ignorelist = [];
end

if ismember(mode, [1, 3, 9, 10])
    starttime = 5;
    endtime = 20;
elseif ismember(mode, [2, 4])
    starttime = 20;
    endtime = 100;
end

c = squeeze(PatternToCount(pattern(:, :, :, :), starttime, endtime));
c(squeeze(ignorelist)) = nan;
cm = mean(c, 3, 'omitnan');
valtrinum = sum(~isnan(c), 3);
ce = std(c, 0, 3, 'omitnan') ./ sqrt(valtrinum); % get error
valind = valtrinum > 0;
oneind = valtrinum == 1;
ce(oneind) = sqrt(ce(oneind) .* valtrinum(oneind)) ./ valtrinum(oneind); % replace with poisson error
ce(ce == 0) = 1 ./ valtrinum(ce == 0); % replace zero-bins with finite error.
% this procedure should not affect the poisson statistics.
nzerobins = nnz(ce == 0);

[thetas, phis] = meshgrid((0:18:144) * pi / 180, (0:20:180) * pi / 180);
thetas = [thetas(end:-1:6, end:-1:2), thetas(1:5, :)];
phis = [phis(end:-1:6, end:-1:2), phis(1:5, :)];
maxfr = max(cm(:));

if ismember(mode, [3, 4, 5, 6])
    [nthetas, nphis] = flipspherical(thetas, phis); % flip the coordinate.
    angles = {nthetas(valind), nphis(valind)};
    % kappa, beta, elev, azim, alpha, amp, baseline
    initp = [1, 0.2, pi/3, pi/6, pi/2, 1, 0.1];
    b1 = [0, -0.5, 0, -144 / 180 * pi, -2 * pi, 0, 0];
    b2 = [100, 0.5, pi/2, 144 / 180 * pi, 2 * pi, maxfr * 1.5, maxfr * 1.5];
elseif ismember(mode, [9, 10])
    if mode == 10
        valind = valind & (phis == pi | phis == 0);
    end
    ifile = load('ILD_mini_082018.mat');
    angles = {thetas(valind), phis(valind), ifile.ILD_mini(valind)};
    initp = [1, 1, 1, 1];
    % softplus
    %b1 = [-maxfr, -maxfr, -maxfr, 0.1, 0] * 3;
    %b2 = [maxfr * 3, maxfr * 3, maxfr * 3, 100, 2];
    b1 = [-maxfr, -maxfr, -maxfr, 0.1] * 3;
    b2 = [maxfr * 3, maxfr * 3, maxfr * 3, 100];
else
    angles = {thetas(valind), phis(valind)};
    initp = [1, 0.5, pi / 3, pi / 6, pi / 2, 1, 0.1];
    b1 = [0, -0.5, 0, 0, -2 * pi, 0, 0];
    b2 = [100, 0.5, pi, pi, 2 * pi, maxfr * 1.5, maxfr * 1.5];
end


if ismember(mode, [9, 10])
    fitfunc = @AuditoryModelFunc;
else
    fitfunc = @kentdist;
end

res = MinuitFunctionFit(fitfunc, angles, cm(valind), initp, ce(valind), 0, b1, b2, 0);

%res
%res.params_err
%diag(sqrt(res.cov_err))'

flatdist = @(p, x) p(1);
initp = 1;
b1 = 0;
b2 = 1000;
res2 = MinuitFunctionFit(flatdist, angles, cm(valind), initp, ce(valind), 0, b1, b2, 0);


if res2.goodfit_err
    res.flatnll2 = res2.nll2;
    res.flatnll2red = res2.nll2red;
    res.flatnll2sig = res2.nll2sig;
else
    res.flatnll2 = inf;
    res.flatnll2red = inf;
    res.flatnll2sig = inf;
    %keyboard; % do investigation
    %error('fit failed');
end
res.passdata.valind = valind;
res.nzerobins = nzerobins;
