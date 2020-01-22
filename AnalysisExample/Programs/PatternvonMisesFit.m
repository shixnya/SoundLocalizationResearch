function fitres = PatternvonMisesFit(theta, singlepattern, trange)
% function PatternvonMisesFit(singlepattern)
%
%    singlepattern - {cell} A cell array that contains trial pattern in 4th dimension.
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
% Created: 1/1/2018
% Modified: 

pat_hist = PatternToFRhist(singlepattern, [5, 25]);

mat_hist = cell2mat(pat_hist);
nTrial = size(mat_hist, 4);

minallowed = 1 / nTrial; % minimum allowed value of SEM

errvals = max(std(mat_hist, 1, 4) / sqrt(nTrial), minallowed); % SEM capped lower side
errpat = ErrorNum(mean(mat_hist, 4), errvals);

y = squeeze(errpat);
fitres = vonMisesFit1DEN(theta, y);

if isfield(fitres, 'params_err')
    fitres.resEN = ErrorNum(fitres.params_err, sqrt(diag(fitres.cov_err))');
end
