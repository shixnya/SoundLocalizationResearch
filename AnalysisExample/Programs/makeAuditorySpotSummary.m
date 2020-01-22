function makeAuditorySpotSummary(datafolder, stimlab, segnum)
% make a summary file.
datafolder = slashappend(datafolder);
usefilename = [datafolder, 'SpotInd.mat'];

if exist(usefilename, 'file')
    load(usefilename);
else
    useind = 1;
end
load([datafolder, 'asdf']);
load([datafolder, 'segttls']);
%ttls = segttls{useind};
ttls = segttls{segnum};


% here it defines overall default parameters.
stimtype = stimlab.stimname;
fittype = 'None';
siglevel = 0.001;
minspike = 10;
limitelev = 0;
windows = [5, 20; 20, 100; 105, 120];
spontwin = [950, 1950];
asdf_loco = getLocomotion(datafolder);

% here it defines special override parameters.
if contains(stimtype, 'fullfield')
    % full field stimulus.
    fittype = 'Kent';
elseif contains(stimtype, 'horizontal_pupcalls_5')
    % horizontal pupcall 5 stimulus use different length.
    windows = [5, 25; 25, 60; 60, 200; 5, 1200]; % add a longer one
    spontwin = [1450, 1950];
    fittype = 'vonMises';
elseif contains(stimtype, 'horizontal')
    % other horizontal stimulus. use the same default structures
    % actually, there is nothing to override. if 2d fit is implemented.
    % it will be listed here.
    if mod(length(ttls), 17) == 0
        fittype = 'vonMises';
    end
end


%ret = AnalyzeAuditorySpot(asdf_raw, ttls, stimtype, dofit, siglevel, minspike);
ret = AnalyzeAuditorySpot(asdf_raw, ttls, stimtype, fittype, siglevel,...
    minspike, limitelev, windows, spontwin, 'asdf_loco', asdf_loco');
saveStruct([datafolder, 'AuditorySpotSummary_' num2str(segnum)], ret);
