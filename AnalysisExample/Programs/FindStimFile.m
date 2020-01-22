function [stimfilename trackballfilename] = FindStimFile(datatime, stimfolder, prefix, maxtime)
%function [stimfilename stimindex] = FindStimFile(datatime, stimfolder);
%
%    datatime - (1,1) Matlab serial date number of the data recording
%    stimfolder - (String) a folder name that contains a bunch of stim
%    files (P_*)
%    prefix - (String) One letter that defines stimulus type. For visual
%    stimulus, it is P, trackball is T, auditory stim is A.
%
% Returns:
%    stimfilename - (String) exact file name that is closest in time with
%    data.
%    trackballfilename - (String) trackball file name for spherical
%    treadmill. Starts with T_.
%    stimindex - (1,1) index of the stimfile in that folder
%
% Description:
%    
%
% Example:
%    
%
% Requires:
%    slashappend.m
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/11/2014
% Modified: 10/12/2014 for adding trackball file
% Modified: 9/10/2018 for adding Auditory analysis and max time


if nargin < 3
    prefix = 'P'; % the default is the visual experiment from PsychStimController.
end

%stimfolder = '/Volumes/SSD/HouseInVivo/2014-09-08-0/08-Sep-2014';
%datatime = serialtime;

debug = 1;


stimfolder = slashappend(stimfolder);


% find stim files

stimfiles = dir([stimfolder prefix '_*.*']); % all the parameters files

if isempty(stimfiles)
    error('Stim Files not found!!!!')
end
for i = 1:length(stimfiles)
    stimtimes(i) = FilenameToTime(stimfiles(i).name);
    difftime(i) = stimtimes(i) - datatime;
end

[minval, stimindex] = min(abs(difftime));

if minval > 0.003
    warning('The minimum difference of time between stim and data is more than minutes. Check if it is OK.');
elseif debug
    fprintf('Stim file chosen: %s, time difference: %ss\n', stimfiles(stimindex).name, num2str(minval * 24 * 60 * 60));
end


trackname = stimfiles(stimindex).name;
trackname(1) = 'T';
stimfilename = [stimfolder stimfiles(stimindex).name];
trackballfilename = [stimfolder trackname];
if ~exist(trackballfilename, 'file')
    disp('Trackball file does not exist.')
    trackballfilename = 0;
end
