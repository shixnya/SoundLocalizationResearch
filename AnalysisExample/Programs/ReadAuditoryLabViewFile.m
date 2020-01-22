function stimlab = ReadAuditoryLabViewFile(filename)
% function auditorystim = ReadAuditoryLabViewFile(filename)
%
%    filename - (String) .txt filename from labview.
%
% Returns:
%    stiminfo.
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
% Created: 9/10/2018
% Modified: 

%filename = stimfiles{1};

% this is a 4 line file with 
% stimulus folder name (also corresponds with stim file name)
% stimulus 
% minimum interval of the stimulus (usually 2s)
% type of LabView Filter (0: no filter, 1: 4k high-pass, 2: 4-50k 3: 4-20k)

fid = fopen(filename);
l1 = fgetl(fid);
stimlab.stimname = strrep(l1(18:end), '\', filesep);
l2 = fgetl(fid);
stimlab.multipliers = str2num(l2);
l3 = fgetl(fid);
stimlab.min_interval = str2num(l3);
l4 = fgetl(fid);
stimlab.filtertype = str2num(l4);
fclose(fid);


