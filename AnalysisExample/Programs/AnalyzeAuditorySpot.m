function ret = AnalyzeAuditorySpot(asdf, ttls, varargin) % stimtype, fittype, siglevel, minspike, limitelev)
% function ret = AnalyzeAuditorySpot(asdf, ttls, stimtype, fittype, siglevel, minspike, limitelev)
%
%    a - (c,1)
%    b - {d,2}
%
% Returns:
%
%
% Description:
%
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
% Created: some time in 2018
% Modified: 8/14/2018
p = inputParser;
%p.addRequired('asdf', @iscell);
%p.addRequired('ttls', @isnumeric);
p.addOptional('stimtype', [], @ischar);
p.addOptional('fittype', 'None', @ischar);
p.addOptional('siglevel', 0.001, @isnumeric);
p.addOptional('minspike', 10, @isnumeric);
p.addOptional('limitelev', 0, @isnumeric);
p.addOptional('windows', [5, 20; 20, 100; 105, 120], @isnumeric);
p.addOptional('spontwin', [950, 1950], @isnumeric);
p.addParameter('asdf_loco', [], @iscell)

p.parse(varargin{:});
pr = p.Results;


%si = DetermineAuditoryStim(ttls, stimtype); % getting stimulation info
si = AuditoryStimInfo(pr.stimtype);
nNeu = asdf{end}(1);
ntypes = length(ttls) / si.npattern / si.nreps;

% calculate firing rate in each integration windows.
%windows = [5, 25; 25, 60; 60, 200];
nwin = size(pr.windows, 1);
win_fr = cell(nwin, 1);
nseg = length(ttls);

[spont_fr, spont_disper] = SegmentFR(asdf, ttls, pr.spontwin(1), pr.spontwin(2));
for i = 1:nwin
    win_fr{i} = SegmentFR(asdf, ttls, pr.windows(i, 1), pr.windows(i, 2));
end


% get p-values for each window
winsize = diff(pr.windows, 1, 2);
pvals = zeros(nNeu, nwin);
ac_sp = zeros(nNeu, nwin);
%minspike = 10; % a parameter. minimum number of spike to consider

for i = 1:nwin
    lambda = spont_fr * winsize(i) / 1000 * nseg; % expected number of spikes
    ac_sp(:, i) = win_fr{i} * winsize(i) / 1000 * nseg; % actual number of spikes
    %pvals(:, i) = poisscdf(ac_sp(:, i), lambda);
    pvals(:, i) = poisscdf(ac_sp(:, i) ./ spont_disper, lambda ./ spont_disper);
    % divide by dispersion to get overdispersed poisson statistics.
end

% get cells with a positive response
%figure;histogram(pvals, 200)
pospval = 1-pvals;
possig = pospval < pr.siglevel / nwin & ac_sp >= pr.minspike;
negsig = pvals < pr.siglevel / nwin & ac_sp >= pr.minspike;
posneu = any(possig, 2);
negneu = any(negsig, 2);
winnumsig = zeros(nwin, 1);
winnumsig_n = zeros(nwin, 1);

posneu1 = find(pospval(:, 1) < pr.siglevel);
%posneu2 = find(any(pospval(:, 2:3) < pr.siglevel / 2, 2));
posneu2 = find(pospval(:, 2) < pr.siglevel);
posneu3 = find(pospval(:, 3) < pr.siglevel);

for i = 1:nwin
    winnumsig(i) = nnz(possig(:, i));
    winnumsig_n(i) = nnz(negsig(:, i));
end
%winnumsig;
%nnz(posneu);
%nnz(negneu);

%% do fit analysis
asdf_pos = ASDFSubsample(asdf, find(posneu));
patterns = FormatAuditoryPSTH(asdf_pos, ttls, si, 200);
patterns_2s = FormatAuditoryPSTH(asdf_pos, ttls, si, 2000);
patterns_all = FormatAuditoryPSTH(asdf, ttls, si, 2000);

[projs, projs_pn] = AuditoryPSTHAnalysis(patterns_2s, spont_fr(posneu), 13, 0.05, 500);

%% after getting the patterns, do running analysis.

if ~isempty(pr.asdf_loco)
    patterns_loco = FormatAuditoryPSTH(pr.asdf_loco, ttls, si, 500);
    running = PatternToCount(patterns_loco.fullpat, 0, 500) > 3.33;
    
    % define spontaneous firing rate for each condition
    nPos = nnz(posneu);
    spont_fr_run = zeros(nPos, 1);
    spont_fr_stat = zeros(nPos, 1);
    for i = 1:nPos
        wsize = diff(pr.spontwin);
        pat = patterns_2s.fullpat(:, :, :, :, i);
        spont_fr_run(i) = mean(cell2mat(PatternToFRhist(pat(running), pr.spontwin))) * 1000 / wsize;
        spont_fr_stat(i) = mean(cell2mat(PatternToFRhist(pat(~running), pr.spontwin))) * 1000 / wsize;
    end
    
    [projs_stat, projs_pn_stat] = AuditoryPSTHAnalysis(patterns_2s, spont_fr_stat, 13, 0.05, 500, ~running);
    [projs_run, projs_pn_run] = AuditoryPSTHAnalysis(patterns_2s, spont_fr_run, 13, 0.05, 500, running);
    % a little special treatment here for renormalizing the stat/run
    normfactor = repmat(1 ./ mean(max(projs, 0), 2), 1, size(projs, 2));
    projs_renorm_stat = max(projs_stat, 0) .* normfactor;
    projs_renorm_run = max(projs_run, 0) .* normfactor;
    
else
    patterns_loco = [];
    running = [];
    spont_fr_stat = [];
    projs_stat = [];
    projs_pn_stat = [];
    projs_renorm_stat = [];
    spont_fr_run = [];
    projs_run = [];
    projs_pn_run = [];
    projs_renorm_run = [];
end
%%

%spatpat = PatternToFRhist(patterns.fullpat(1, :, :, :, :), [5, 25, 60, 200]);
%newspatmat = permute(cell2mat(permute(spatpat, [2,1,3,4,5])), [2,1,3,4,5])

sf.nSteps1 = length(si.meta2); % nazims
sf.nSteps2 = length(si.meta1); % nelevs
%patterns.fullpat = psthmat;
if strcmp(pr.fittype, 'UMLE') % normal unbinned MLE fit
    patterns_fit = patterns;
    if pr.limitelev
        patterns_fit.fullpat = patterns.fullpat(:, pr.limitelev, :, :, :);
    end
    fitresult = calculateAllSpotDataset(patterns_fit, sf, 5);
elseif strcmp(pr.fittype, 'Kent') % kent distribution fit
    asdf_pos1 = ASDFSubsample(asdf, posneu1);
    patterns1 = FormatAuditoryPSTH(asdf_pos1, ttls, si, pr.windows(1,2));
    nNeu1 = asdf_pos1{end}(1);
    nstim = size(patterns1.fullpat, 1);
    fitresult_stims = cell(nstim, 3);
    for i = 2:nstim
        fitresult_stims{i} = cell(nNeu1, 1);
    end
    fitresult = cell(4, 1);
    fitresult{1} = cell(nNeu1, 1);
    fitresult{3} = cell(nNeu1, 1);
    fitresult{9} = cell(nNeu1, 1);
    fitresult{10} = cell(nNeu1, 1);
    if ~isempty(running)
        % 5: stationary, 6: running;
        fitresult{5} = cell(nNeu1, 1);
        fitresult{6} = cell(nNeu1, 1);
        fitresult{7} = cell(nNeu1, 1);
        fitresult{8} = cell(nNeu1, 1);
    end
    
    
    asdf_pos2 = ASDFSubsample(asdf, posneu2);
    patterns2 = FormatAuditoryPSTH(asdf_pos2, ttls, si, pr.windows(2,2));
    nNeu2 = asdf_pos2{end}(1);
    fitresult{2} = cell(nNeu2, 1);
    fitresult{4} = cell(nNeu2, 1);
    firsthalf = running;
    firsthalf(1, :, :, 1:(end/2)) = 1;
    firsthalf(1, :, :, (end/2+1):end) = 0;
            
    for i = 1:nNeu1
        fitresult{1}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 1);
        fitresult{3}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 3);
        fitresult{9}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 9);
        fitresult{10}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 10);
        if ~isempty(running)
            fitresult{5}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 3, running(1,:,:,:));
            fitresult{6}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 3, ~running(1,:,:,:));
            fitresult{7}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 3, ~firsthalf(1,:,:,:));
            fitresult{8}{i} = KentFitOneCell(patterns1.fullpat(1,:,:,:,i), 3, firsthalf(1,:,:,:));
        end
        if nstim > 1 % do additional fits
            for j = 2:nstim
                fitresult_stims{j, 1}{i} = KentFitOneCell(patterns1.fullpat(j,:,:,:,i), 3);
                fitresult_stims{j, 2}{i} = KentFitOneCell(patterns1.fullpat(j,:,:,:,i), 9);
                fitresult_stims{j, 3}{i} = KentFitOneCell(patterns1.fullpat(j,:,:,:,i), 10);
            end
        end
    end
    
    fitresult_stims{1,1} = fitresult{3};
    fitresult_stims{1,2} = fitresult{9};
    fitresult_stims{1,3} = fitresult{10};
    for i = 1:nNeu2
        fitresult{2}{i} = KentFitOneCell(patterns2.fullpat(1,:,:,:,i), 2);
        fitresult{4}{i} = KentFitOneCell(patterns2.fullpat(1,:,:,:,i), 4);
    end
elseif strcmp(pr.fittype, 'vonMises') % von Mises fit.
    asdf_pos1 = ASDFSubsample(asdf, posneu1);
    patterns1 = FormatAuditoryPSTH(asdf_pos1, ttls, si, 25);
    nNeu1 = asdf_pos1{end}(1);
    ntype = size(patterns1.fullpat, 1);
    fitresult = cell(2, ntype);
    fitresult_stims = {};
    for i = 1:ntype
        fitresult{1, i} = cell(nNeu1, 1);
        fitresult{2, i} = cell(nNeu1, 1);
        for n = 1:nNeu1
            theta = ((-144:18:144) / 180 * pi)';
            fitresult{n, i} = PatternvonMisesFit(theta, ...
                patterns1.fullpat(i, :, :, :, n), [5, 25]);
        end
    end
elseif strcmp(pr.fittype, 'None')
    fitresult = [];
    fitresult_stims = [];
else
    error('fittype is wrong. fittype: %s', pr.fittype);
end


%% display some results
fprintf('nNeu: %3d, IW1 %3d neurons,     IW2: %3d neurons,      IW3: %3d neurons,',...
    nNeu, winnumsig(1), winnumsig(2), winnumsig(3));
fprintf('     %3d total neurons with a positive response\n', nnz(posneu));

fprintf('           IW1 %3d neurons,     IW2: %3d neurons,      IW3: %3d neurons,',...
    winnumsig_n(1), winnumsig_n(2), winnumsig_n(3));
fprintf('     %3d total neurons with a negative response\n', nnz(negneu));


%else%patterns;

%% return packed information
ret = v2struct(pr, spont_fr, spont_disper, ac_sp, possig, negsig,...
    posneu, negneu, pvals, pospval, si,  fitresult, projs, projs_pn,...
    posneu1, posneu2, posneu3, patterns, patterns_2s, patterns_all,...
    patterns_loco, running, projs_stat, projs_run, projs_pn_stat,...
    projs_pn_run, projs_renorm_stat, projs_renorm_run,...
    spont_fr_stat, spont_fr_run, fitresult_stims);
%saveStruct('AuditorySpotSummary', ret);


