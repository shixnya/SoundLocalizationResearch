function FRhist = PatternToFRhist(pat, duration, binsize)
% function FRhist = PatternToFRhist(pat, duration, binsize)
% OR
% function FRhist = PatternToFRhist(pat, range)
%
%
%    pat - {nD cell array} PSTH pattern
%    duration - (scalar) endtime of the PSTH in ms
%    binsize - (scalar) bin size in ms
%    OR
%    range - (array) edges of the histogram used for histcounts in ms
%
% Returns:
%    FRhist - Histcounts of firing rates
%
%
% Description:
%    in 3 argument mode where you define one bin size, it returns firing
%    rate in Hz. Otherwise, it does not normalize
%
% Example:
%    
%
% Requires:
%    none
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/21/2018
% Modified: 9/21/2018 added explanation

if nargin == 2
    range = duration;
    normfac = 1;
else
    range = (0:binsize:duration);
    normfac = binsize / 1000;
end
histhandle = @(x) histcounts(x, range) / normfac;
FRhist = cellfun(histhandle, pat, 'UniformOutput', 0);
