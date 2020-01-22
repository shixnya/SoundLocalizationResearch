function serialtime = FilenameToTime(filename)
%function serialtime = FilenameToTime(filename);
%
%    filename - (String) file that has Intan date format
%
% Returns:
%    serialtime - (1,1) matlab serial date number
%
% Description:
%    a simple function to convert date fileformat that we use for
%    parameters files and recording files to serial date format in matlab.
%    This will be used to compare time stamps and choose appropriate stim
%    information files.
%
% Example:
%    filename = 'R_UCSC355A_140908_145017';
%    timestamp = FilenameToTime(filename);
%
% Requires:
%    nothing
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/11/2014
% Modified: 

%filename = 'R_UCSC355A_140908_145017.rhd';
stampname = filename(end-16:end-4);
serialtime = datenum(stampname, 'yymmdd_HHMMSS');