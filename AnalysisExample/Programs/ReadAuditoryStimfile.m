function stiminfo = ReadAuditoryStimfile(filename)

file = fopen(filename);
% first two lines should be stim meta info
firstline = fgetl(file);

if strcmp(firstline, 'Auditory stimulus definition file, version 2')
    % version 2 analysis
    stiminfo.ntypes = str2num(fgetl(file));
    stiminfo.meta1 = str2num(fgetl(file));
    stiminfo.meta2 = str2num(fgetl(file));
else % version 1 without definition head line
    stiminfo.ntypes = 1;
    stiminfo.meta1 = str2num(firstline);
    stiminfo.meta2 = str2num(fgetl(file));
end

stimnums = [];
while 1
    tline = fgetl(file);
    if ~ischar(tline), break, end
    stimnums = [stimnums str2num(tline)];
end
fclose(file);

stiminfo.stimnums = stimnums;
l1 = length(stiminfo.meta1);
l2 = length(stiminfo.meta2);
lt = stiminfo.ntypes;


if l1 == l2 % this is an ad-hoc code. likely be flattened (serial) stim.
    stiminfo.nreps = length(stimnums) / l1 / lt;
    stiminfo.serial = 1;
    stiminfo.npattern = l1 * lt;
else
    stiminfo.nreps = length(stimnums) / l1 / l2 / lt;
    stiminfo.serial = 0;
    stiminfo.npattern = l1 * l2 * lt;
end