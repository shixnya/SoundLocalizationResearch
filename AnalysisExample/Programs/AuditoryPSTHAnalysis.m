function [projs, pos_norm_projs] = AuditoryPSTHAnalysis(patterns, spont_fr, ord, binsize, maxtime, running)
% function projs = AuditoryPSTHAnalysis(patterns, ord, binsize, maxtime)
% project the auditory PSTH to the cosyne bump basis
% spont_fr needs to be prepared for neurons selected for patterns.


if nargin < 3
    ord = 12;
end
if nargin < 4
    binsize = 0.05; % default 20 kHz data.
end
if nargin < 5
    maxtime = 1500; % maxtime is set to 500 ms.
end
if nargin < 6
    running = [];
end

% analysis of PSTH of auditory neurons.

%datafolder = slashappend(pwd);
%load([datafolder, 'segttls']);
%stiminfos = cellfun(@DetermineAuditoryStim, segttls, 'UniformOutput', 0);

%%
%load([datafolder 'asdf']);
%ret = AnalyzeAuditorySpot(asdf_raw, segttls{1}, 0, 0.001, 10);


%% make a PSTH of all the neurons. (all positions together)

%nNeuPos = size(ret.patterns.fullpat, 5);
nNeuPos = size(patterns.fullpat, 5);
%ord = 12;
%binsize = 0.05; % ms
%maxtime = 500;



timerange = 0:binsize:(maxtime+binsize); % characterize between 0 and 300 ms.
timerange(end) = timerange(end) - 0.0001; % to avoid right edge issue.
%projs = zeros(nNeuPos, ord + 1); % 12 dim vector for each neuron.
projs = zeros(nNeuPos, ord); % 12 dim vector for each neuron.
%projs_i = zeros(ord + 1, nNeuPos);
%[ibasis, obasis, basis] = cosyneBumpBasis(pi, -5, 0:binsize:maxtime, ord, 1); % orthonormal cosyne bump basis with default parameters.
basis = cosyneBumpBasis(pi, -5, 0:binsize:maxtime, ord, 0); % orthonormal cosyne bump basis with default parameters.
fpsize = size(patterns.fullpat);
nstim = prod(fpsize(2:4));
%posneuind = find(ret.posneu);
for i = 1:nNeuPos
    p = patterns.fullpat(:, :, :, :, i);
    if isempty(running)
        allspikes = cat(2, p{:});
    else
        allspikes = cat(2, p{running});
    end
    %spont_spikes = ret.spont_fr(posneuind(i)) * 0.05 * 0.001 * nstim; % expected number of spikes per bin.
    spont_spikes = spont_fr(i) * binsize * 0.001 * nstim;
    hc = histcounts(allspikes, timerange);
    hc = hc - spont_spikes;
    %hc = max(hc, 0);
    if sum(hc) < 0
        %warning(sprintf('Neuron:%3d, negative response', i));
    end
    %hc = hc / sum(hc); % normalization by spike count.
    %hc = hc / sum(max(hc, 0)); % normalization only by positive spike count.
    hc = hc / nnz(allspikes);
    %projs(i, :) = hc * obasis;
    projs(i, :) = hc * basis;
    %projs_i(:, i) = ibasis * hc';
end

% make positive-normalized projections.
posprojs = max(projs, 0);
pos_norm_projs = normalizeRows(posprojs, 'mean');



%% Do clustering analysis

%ClusterBrusher(projs_b)