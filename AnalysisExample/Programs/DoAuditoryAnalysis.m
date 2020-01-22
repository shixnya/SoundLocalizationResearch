function DoAuditoryAnalysis(datafolder, doind)
% function DoAuditoryAnalysis(datafolder, doind)
% Do Auditory Analysis by automatic stimulus detection
if nargin < 1
    datafolder = pwd;
end
if nargin < 2
    doind = 0; % do all
end
datafolder = slashappend(datafolder);


% this analysis requires the user to put their analysis files in
%auditorystimpath = '~/AuditoryAnalysisParams/';
auditorystimpath = pwd;

% get the time stamp from the current directory.
load([datafolder, 'segmentlengths.mat']);
% this load 'timestamps', 'segmentlengths', 'segmentseparations'
sttlname = [datafolder, 'segttls.mat'];
if ~exist(sttlname) % if not there, make it.
    makesegttls(datafolder);
end
load([datafolder, 'segttls']) % load it.

stimfiles = arrayfun(@(x) FindStimFile(x, auditorystimpath, 'A'),...
    timestamps, 'UniformOutput', false);
stimlab = cellfun(@ReadAuditoryLabViewFile, stimfiles);


if doind == 0
    doind = 1:length(stimlab);
end

for i = doind
    makeAuditorySpotSummary(datafolder, stimlab(i), i);
end

