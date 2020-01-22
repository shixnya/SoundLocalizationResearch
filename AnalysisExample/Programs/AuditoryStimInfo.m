function stiminfo = AuditoryStimInfo(stimtype)
%function AuditoryStimInfo(ttls, stimtype)
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
% Created: 8/14/2018
% Modified: 

%%
%stimtype = 'gen4/gen4_horizontal_pupcall_1'
%stimtype = 'gen4/loomingspot_gen4_noise_triangle'
%auditoryfolder = '~/python/Auditory/stimgen/';
auditoryfolder = slashappend(pwd);
sname = [auditoryfolder stimtype '.txt'];
stiminfo = ReadAuditoryStimfile(sname);

