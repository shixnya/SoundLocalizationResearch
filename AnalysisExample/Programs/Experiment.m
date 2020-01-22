classdef Experiment < handle
    % A class for managing expriment
    % This class has number of methods for loading/calculating properties
    % of the cells. The class also implements a number of filters and limit
    % the use of neurons by those filters.
    % 
    
    
    
    properties
        datafolder
        movementfiles
        arrayID
        segttls
        neurons; % properties of the neurons
        expprop; % properties of the experiment (stim num, etc)
        useIDs; % indicies used for the analysis (not logical)
        nNeuTotal;
        filenamemap; % mapping of the features to filename associated with it.
        % ok. undefined properties are initialized as [] (0x0 double).
        wpobj; % wavelet peak object. Set by loading it.
        spot_patterns;
        exploratory; % exploratory data before unblinding
        blinded; % blinded data. Unblind
    end
    methods
        function obj = Experiment(datafolder)
            % Make an object out of data folder.
            % It tries to read the depth from thefolder and return only
            % depth >= 0 neurons by default
            % If you want all the neurons you can use resetUseIDs().
            obj.datafolder = slashappend(datafolder);
            % by default, do not use neurons above the surface

            depth = getDepth(datafolder);
            obj.nNeuTotal = length(depth);
            obj.useIDs = find(depth >= 0);
            obj.neurons = struct('depth', depth);
            obj.wpobj = [];
            obj.spot_patterns = [];
        end
        
        function n = nNeu(obj)
            % Returns the number of neurons currently in use
            n = length(obj.useIDs);
        end
        
        function resetUseIDs(obj)
            % reset use IDs and use all the neurons
            obj.useIDs = (1:obj.nNeuTotal)';
        end
        
        function blind(obj)
            % hide half of the neurons from this dataset.
            nNeu = obj.nNeu;
            nNeuh = floor(nNeu/2);
            perm = randperm(nNeu);
            
            obj.exploratory = sort(obj.useIDs(perm(1:nNeuh)));
            obj.blinded = sort(obj.useIDs(perm(nNeuh+1:end)));
            obj.useIDs = obj.exploratory;
        end
        
        function blindswitch(obj, state)
            if state == 1
                obj.useIDs = obj.exploratory;
            elseif state == 2
                obj.useIDs = obj.blinded;
            elseif state == 3
                obj.useIDs = sort([obj.exploratory; obj.blinded]);
            else
                error('Unknown state specified. Choose from 1 to 3.');
            end

        end
        
        
        function class = classifyLocation(obj, border, th_depth)
            % classify neurons into two categories
            % 1: anterior map, double zone (superficial)
            % 2: posterior map, double zone (superficial)
            % 3: anterior map, single zone (deep)
            % 4: posterior map, signel zone (deep)
            % 5: border area
            dep = obj.load('cdepth');
            shanks = obj.getShanks;
            class = zeros(obj.nNeu, 1);
            
            if length(border) == 1 % border on the shank
                class(shanks < border) = 1;
                class(shanks > border) = 2;
                class(dep > th_depth) = class(dep > th_depth) + 2;
                class(shanks == border) = 5;
            else % border between the shanks
                class(shanks <= border(1)) = 1;
                class(shanks >= border(2)) = 2;
                class(dep > th_depth) = class(dep > th_depth) + 2;
            end
        end
        
        
        
        % get movement, given the data number
        function mov = getMovement(obj, datanumber)
            obj.readBasicInfo;
            mov = obj.movementfiles{datanumber};
        end
        
        % run this before using segttls.
        function readBasicInfo(obj)
            if isempty(obj.movementfiles)
                load([obj.datafolder 'basicinfo'])
                obj.movementfiles = movement;
                obj.arrayID = arrayID;
                obj.segttls = segttls;
            end
        end
        
        % get vision ID of a neuron
        function vIDs = getVisionID(obj, nnum)
            load([obj.datafolder, 'asdf'], 'IDs');
            vIDs = IDs(obj.useIDs(nnum));
        end
        
        % get stimulus based time series?
        function stimTS(obj, stimulus)
            fullpatfilename = [stimulus '_fullpat'];
            fp = load([obj.datafolder fullpatfilename], 'patterns', 'sf');
            % incomplete. TODO?
        end
        
        % simple function for getting file saved in a folder
        function save(obj, variable, variablename, filename)
            eval([variablename ' = variable;']);
            save([obj.datafolder, filename], variablename)
        end
        
        % returns if the neurons are properly measured.
        function bad = loadp(obj, variablename, filename)
            if nargin < 3
                if strcmp(variablename, 'nlind')
                    filename = 'Cphasegratings_fullpat';
                else
                    filename = obj.getFilename(variablename);
                end
            end
            
            
            try
                load([obj.datafolder filename], 'stimnum');
            catch
                warning('file missing: %s', filename);
                bad = logical(ones(length(obj.useIDs), 1));
                return
            end
            if exist('stimnum', 'var')
                load([obj.datafolder 'FRACresults.mat'], 'badsegments_startend');
                bad = logical(badsegments_startend(obj.useIDs, stimnum));
            else
                warning('stimnum is strange: %s', filename);
                bad = logical(ones(length(obj.useIDs), 1));

                %error('stimnum was not found');
            end
        end
        
        function bad = stimulusFilter(obj, stimnum)
            load([obj.datafolder 'FRACresults.mat'], 'badsegments_startend');
            bad = logical(badsegments_startend(obj.useIDs, stimnum));
            obj.useIDs(bad) = [];
        end
        

        % returns the stimhulation number, given the name of the data such
        % as 'DriftGratings'
        function stimnum = getStimnum(obj, dataname)
            load([obj.datafolder dataname], 'stimnum');
        end
        
        % get spontaneous firing rate of the cells, given dataname
        % not working
        function [sfr, sfr_e] = getSFR(obj, dataname, offset)
            pats = load([obj.datafolder dataname], 'patterns', 'sf');
            [sfr, sfr_e] = CalculateSFR(pats.patterns, pats.sf.WaitInterval, offset);
        end
            
        function indv_patterns = getPatterns(obj, filename, varname)
            if nargin < 3
                varname = 'patterns';
            end
            indv_patterns = cell(obj.nNeu, 1);
            orig = obj.load(varname, filename);
            for i = 1:obj.nNeu

                indv_patterns{i} = orig.fullpat(:, :, :, :, obj.useIDs(i));
            end
        end
        
        function fullpat = getAuditoryAllPatterns(obj, patnum)
            if nargin < 2
                patnum = 1;
            end
            orig = obj.load('patterns_all', sprintf('AuditorySpotSummary_%d', patnum));
            fullpat = permute(orig.fullpat(:,:,:,:,obj.useIDs), [5 1 2 3 4]);
            % bring the neurons first for concatenation
            %fullpat = permute(fullpat, [5 1 2 3 4]);
        end
        
        function sfs = getSFs(obj, filename)
            sfs = cell(obj.nNeu, 1);
            sf = obj.load('sf', filename);
            [sfs{:}] = deal(sf);
        end

        
        % getting partial ASDF
        function [asdf, ttls] = getPartialASDF(obj, stimnum)
            load([obj.datafolder 'segmentlengths.mat'], 'segmentseparations');
            load([obj.datafolder 'asdf.mat'], 'asdf_raw');
            load([obj.datafolder 'ttlTimes'])

            offsets = [0 segmentseparations];
            offset = offsets(stimnum);
            
            st = offsets(stimnum); % starting time
            et = offsets(stimnum + 1); % ending time
            
            asdf = ASDFChooseTime(asdf_raw, st, et);
            asdf = ASDFSubsample(asdf, obj.useIDs);
            
            ttls = ttlTimes(ttlTimes > st & ttlTimes < et) - st;
        end
        
        
        % load wavelet peak object.
        % if an argument is not given, load it from base.
        % the argument can be a number (stimnum) or
        % String (filename that stores stimnum)
        function loadWaveletPeaks(obj, loadlocation)
            if nargin < 2 % read from base
                obj.wpobj = WaveletPeaks([obj.datafolder 'wpeaks8.mat']);
            else % determine what wpeaks it wants to load
                % get a list of names of the breakdown folders
                indv = obj.getIndividualFolders;
                basefol = obj.datafolder;

                if isnumeric(loadlocation) % stimnum
                    obj.wpobj = WaveletPeaks([basefol indv{loadlocation} '/']);
                elseif ischar(loadlocation) % either filename or direct variable name
                    load([obj.datafolder loadlocation], 'stimnum');
                    obj.wpobj = WaveletPeaks([basefol indv{stimnum} '/']);
                end
            end
            % apply filter with current useID.
            obj.wpobj = obj.wpobj.neuronFilter(obj.useIDs);
        end
        
        function indvfolders = getIndividualFolders(obj)
            bpath = readpath(obj.datafolder);
            lastone = strrep(bpath{end-1}, '/', ''); % chop the last slash
            indvfolders = strsplit(lastone, ',');
        end

        function WFstruct = findWNcells(obj, namestring)
            % only normal mode where txt file is already provided
            wn = FindSpecificWNCells(obj.datafolder, namestring);
            names = fieldnames(wn);
            
            % make correspondence of the cellIDs.
            [commonIDs, commoninds, commoninds2] = intersect(wn.cellIDs, obj.useIDs);
            for i = 1:length(names);
                WFstruct.(names{i}) = wn.(names{i})(commoninds, :);
            end
            WFstruct.cellIDs = commoninds2;
        end
        
        function [goodcells, linearinds, visIDs] = chooseCellandRF(obj)
            [goodcells, linearinds, visIDs] = ChooseCellandRF(obj);
        end
        
        function possig = getAuditoryPossig(obj, num, dnum)
            if nargin < 3
                dnum = 1;
            end
            if nargin < 2
                num = 1;
            end
            posnum = sprintf('posneu%d', num);

            lob = load([obj.datafolder sprintf('/AuditorySpotSummary_%d.mat', dnum)], posnum);
            possig = false(obj.nNeuTotal, 1);
            possig(lob.(posnum)) = 1;
        end
        
        % this method ensures that the properties are read just once.
        % I don't think it is implemented. Should it check fields of
        % neurons and expprop?
        function property = load(obj, variablename, filename, ignoremissing)
            if nargin < 4
                ignoremissing = 1; % give only warning by default.
            end
            
            % isn't there a better way?
            if nargin > 2 & isempty(filename)
                clear filename
            end

            try
                % deal with some special cases
                if strcmp(variablename, 'depth')
                    property = getDepth(obj.datafolder);
                elseif strcmp(variablename, 'cdepth')
                    [~, property] = getDepth(obj.datafolder);
                elseif strcmp(variablename, 'nlind')
                    property = getNonLinearityIndex(obj.datafolder);
                elseif strcmp(variablename, 'cx');
                    [~,~,property] = getDepth(obj.datafolder);
                elseif strcmp(variablename, 'cy');
                    [~,~,~,property] = getDepth(obj.datafolder);
                else
                    % all the other properties
                    % if there is a special instruciton for defining the
                    % filename for the given variable, use it.
                    if nargin < 3 & exist([obj.datafolder 'FilenameInstruction.mat'], 'file')
                        load([obj.datafolder 'FilenameInstruction.mat'])
                        map = containers.Map(hashKeys, hashValues);
                        if map.isKey(variablename)
                            filename = map(variablename);
                        end
                    end
                    
                    % otherwise, find from the pre-defined list
                    if ~exist('filename', 'var')
                        filename = obj.getFilename(variablename);
                    elseif strcmp(filename, 'looming')
                        filename = getWildName([obj.datafolder 'LoomingSpots_fullpat*.mat']);
                    end
                    load([obj.datafolder filename], variablename)
                    obj.neurons.(variablename) = eval(variablename);
                    
                    
                    %if nargin > 2 % if filename is specified, read a specific file
                    %    load([obj.datafolder filename], variablename)
                    %    obj.neurons.(variablename) = eval(variablename);
                    %else % try to find a file from known file list.
                    %    filename = obj.getFilename(variablename);
                    %    load([obj.datafolder filename], variablename)
                    %    obj.neurons.(variablename) = eval(variablename);
                    %end
                    property = obj.neurons.(variablename);
                end
            catch ME
                if ignoremissing
                    % if it is SamsungStryker, there might be an
                    % alternative file.
                    %if strcmp(filename, 'SamsungStryker_TBR.mat')
                        
                    %else % normal case
                    property = nan(obj.nNeuTotal, 1); % returning a column vector.
                    warning('variable %s not found in %s', variablename, obj.datafolder);
                    %end
                else
                    rethrow(ME)
                end
            end
            
            
            
            % trim down unused neurons
            if isvector(property)
                if length(property) == obj.nNeuTotal
                    property = property(obj.useIDs);
                    property = property(:); % ensure column vec
                end
            elseif size(property, 1) == obj.nNeuTotal
                property = property(obj.useIDs, :);
            end
        end
        
        function npanel = getSpotPanels(obj)
            sf_spot = obj.load('sf', 'Spot_fullpat');
            if isempty(sf_spot)
                npanel = nan(obj.nNeu, 1);
                return
            end
            np_single = length(sf_spot.orient);
            npanel = np_single * ones(obj.nNeu, 1);
        end
        
        % plot the position of the cells that are specified.
        % This method draw on gca.
        function plotPosition(obj, nids, numtext)
            if nargin < 3
                numtext = 0; % no text by default
            end 
            x = obj.load('x');
            y = obj.load('y');
            aid = obj.load('arrayID');
            try
                cont = obj.load('cont', 'SurfaceContour.mat');
            catch
                cont = [];
            end
            
            LoadVision2;
            cla;
            if isnan(aid)
                aid = 10901;
            end
            TNB_PlotArray(aid, gca, 1); hold on;
            
            if aid == 11001 | aid == 10901
                sign = 1;
            else
                sign = -1;
            end
            if ~(isempty(cont) | all(isnan(cont)))
                plot(-sign * cont(:,1), cont(:,2), '--k');
            end
            plot(sign * y, x, 'b.', 'markersize', 20);
            plot(sign * y(nids), x(nids), '.r', 'markersize', 20);
            if numtext
                x = obj.load('x');
                y = obj.load('y');
                for i = 1:obj.nNeu
                    text(sign * y(i) + 10, x(i), num2str(i));
                end
            end
            

        end
            
        
        % get the firing rate in response to spots.
        % first output is the mean firing rate for each color
        % second output is the total firing rate for each color (mean *
        % npanels)
        function [frate, frate_sum] = getSpotFR_bw(obj)
            frate = zeros(length(obj.useIDs), 2); % initialize before return
            frate_sum = zeros(length(obj.useIDs), 2); % initialize before return
            try
                load([obj.datafolder 'Spot_Result'], 'pxfrate');
            catch
                return
            end
            
            for i = 1:length(obj.useIDs)
                frate(i,1) = mean(mean(pxfrate(:, :, 1, obj.useIDs(i))));
                frate(i,2) = mean(mean(pxfrate(:, :, 2, obj.useIDs(i))));
            end
            
            if nargout > 1
                npx = prod(size(pxfrate(:,:,1,1)));
                frate_sum = frate * npx;
            end
        end
        
        function frate_evoked = getSpotFR_evoked(obj)
            [~, frate_sum] = obj.getSpotFR_bw;
            MLfit = obj.load('MLfit');
            if isempty(MLfit)
                frate_evoked = [];
                return
            end

            inter = SpotFitInterpretation(MLfit);
            sp_b = StructArrayResolve(inter{1});
            sp_w = StructArrayResolve(inter{2});
            
            nNeu = obj.nNeu;
            frate_evoked = 1 / 2 ./ [sp_b.RFarea sp_w.RFarea]...
                .* frate_sum .* [sp_b.response_frac sp_w.response_frac];
        end
            
            
        
        function shanks = getShanks(obj)
            ids = obj.load('IDs');
            arrayID = obj.load('arrayID');
            
            % needs to be defined for each type
            % let's do it just for 128DN and 256A, AN.
            switch arrayID
                case 10701 % 128D
                    shankchan = 32;
                case 10801 % 128DN
                    shankchan = 32;
                case 10901 % 256A
                    shankchan = 64;
                case 11001 % 256AN
                    shankchan = 64;
            end
            
            % depending on how many channels are on the channel, it can
            % determine the shank ID.
            elecnum = ceil(double(ids) / 15);
            shanks = ceil(elecnum / shankchan);
        end

        % a simple method to give a filename given a variablename
        function filename = getFilename(obj, variablename)
            % make a hashmap if it is the first time
            if isempty(obj.filenamemap)
                % load the pre-defined hashmap
                load InVivoDataFileHashMaps
                obj.filenamemap = containers.Map(hashKeys, hashValues);
            end
            filename = obj.filenamemap(variablename);
            
            
            % special treatment for the stryker movie files
            if strcmp(filename, 'SamsungStryker_TBR.mat')
                % check if file exists
                if exist([obj.datafolder filename])
                    return % file exists. OK to return
                % if not, check if other file exists
                elseif exist([obj.datafolder 'StrykerNoise_TBR.mat'])
                    filename = 'StrykerNoise_TBR.mat';
                end
            end
        end
        
        % it returns a list of stim ID that corresponds to the stimname.
        function datanums = findStimulus(obj, stimname)
            load([obj.datafolder 'segmentstims'])
            datanums = find(strcmp(segmentstims, stimname));
        end
        
        function useFRFilterAround(obj, stimnum, thfr)
            % use FR filter around specified area.
            % if stimnum is 3, it apply FR filter for 2 and 4 as well.
            
            obj.load('segfr');
            nSeg = size(obj.neurons.segfr, 2);
            range = stimnum-1:stimnum+1;
            % not complete, but should be practically fine.
            range = setdiff(range, [0 nSeg+1]);
            
            obj.useFRFilter(range, thfr);
        end
        
        function f = AuditoryRunFraction(obj)
            fname = 'AuditorySpotSummary_1';
            load([obj.datafolder, fname], 'running');
            f = mean(running(:));
        end

        
        function h = plotAuditoryRF(obj, neu, patnum, running, windows)
            if nargin < 3
                patnum = 1;
            end
            if nargin < 4
                running = 0; % 0: does not use running
                % 1: show only running
                % 2: show only stationary
            end
            if nargin < 5
                windows = [5 20; 20 100; 105 120; 120 200; 0 200];
            end
            fname = 'AuditorySpotSummary_1';
            nID = obj.useIDs(neu);
            ptn = load([obj.datafolder, fname], 'patterns_all', 'running');
            
            for i = 1:length(patnum)
                figure(1000 + i);clf
                pat1 = ptn.patterns_all.fullpat(patnum(i), :, :, :, nID);
                if running == 1
                    rpat = ptn.running(patnum(i), :, :, :);
                elseif running == 2
                    rpat = ~ptn.running(patnum(i), :, :, :);
                else
                    rpat = [];
                end
                PlotAuditoryPattern(squeeze(pat1), windows, 1, rpat);
            end
        end
        
        function h = plotAuditoryPSTH(obj, neu, dur, incr, fname, runseg, proj)
            % function to plot auditory temporal response.
            if nargin < 5
                fname = 'AuditorySpotSummary_1';
            end
            if nargin < 6
                runseg = 1; % if runseg, segregate running vs stationary
            end
            if nargin < 7
                proj = 0;
            end
            nID = obj.useIDs(neu);
            ptn = load([obj.datafolder, fname], 'patterns_all');
            load([obj.datafolder, fname], 'running');
            if runseg
                states = {~running, running};
            else
                states = {true(size(running))};
            end
            h = {};
            lc = lines(5);
            for i = 1:length(states)
                pat1 = ptn.patterns_all.fullpat(:, :, :, :, nID);
                spike1 = PatternToFRhist(pat1(states{i}), dur, incr);
                cat1 = cat(1, spike1{:});
                en1 = meansem(cat1);
                h{i} = stairsEN(0:incr:dur, en1, 'Color', lc(i,:));
                % spontaneous rate
                counts = PatternToCount(pat1(states{i}), 950, 1950);
                spont = mean(counts);
                h{i}(3) = plot([0, dur], spont * [1 1], '--',...
                    'Color', lc(i,:));
            end
            if runseg
                legend([h{1}(1), h{2}(1)], 'Stationary', 'Running');
                legend boxoff
            end
            if proj
                load([obj.datafolder, fname], 'projs_pn');
                load([obj.datafolder, fname], 'posneu');
                posf = find(posneu);
                ind = find(posf == nID);
                cla
                plot(projs_pn(ind, :));
            end
        end
            
            
        
        function useFRFilter(obj, stimlist, thfr)
            obj.load('segfr');
            % not complete, but should be practically fine.
            badind = find(sum(obj.neurons.segfr(:,stimlist) < thfr, 2));
            obj.useIDs = setdiff(obj.useIDs, badind);
        end

        function variableFilter(obj, varargin)
            bad = obj.loadp(varargin{:});
            obj.useIDs(bad) = [];
        end
        
        function loomingFilter(obj)
            wildloomingname = getWildName([slashappend(obj.datafolder) 'LoomingSpots_fullpat*.mat']);
            %load([slashappend(e.exprs{datanum}.datafolder), wildloomingname], 'patterns_mu', 'sf');
            obj.variableFilter('patterns', wildloomingname);
        end
        
        function getNetworkReport(obj)
            if isempty(obj.wpobj)
                % just give an error.
                error('Wavelet is not loaded. Use loadWaveletPeaks to load it up.')
            end
            % use wpeaks to generate a report.
            xy = [obj.load('x') obj.load('y')];
            load([obj.datafolder 'basicinfo.mat'], 'arrayID');
            obj.wpobj.networkReport(obj.nNeu, xy, arrayID);
        end

        function osstruct = getODSCells(obj, threshold, local)
            if nargin < 2
                threshold = 0.001;
            end
            if nargin < 3
                local = 0;
            end
            
            if local == 0
                filename = 'OSFitResult';
            else
                filename = 'LocalOSFitResult';
            end
            
            
            astat = obj.load('allOSstat', filename);
            ossig = pick(astat, 'ossig');
            dssig = pick(astat, 'dssig');
            responsive = pick(astat, 'maxresponse_sig_corrected');
            
            
            osstruct.pos = ossig < threshold & dssig >= threshold;
            osstruct.ds = dssig < threshold;
            osstruct.responsive = responsive;
            
            osstruct.posres = zeros(obj.nNeu, 1);
            osstruct.negres = zeros(obj.nNeu, 1);
            
            for i = 1:obj.nNeu
                rstat = getFRstat(astat{i});
                if ~isempty(rstat)
                    osstruct.posres(i) = rstat.posres.sig < 0.01 & rstat.posres.value > 0;
                    osstruct.negres(i) = rstat.negres.sig < 0.01 & rstat.negres.value < 0;
                end
            end
            
            osstruct.posres = logical(osstruct.posres);
            osstruct.negres = logical(osstruct.negres);
            
            if isempty(osstruct.pos)
                osstruct.pospos = [];
                osstruct.posneg = [];
                osstruct.dspos = [];
                osstruct.dsneg = [];
                return
            end
            
            osstruct.pospos = osstruct.pos & osstruct.posres & ~osstruct.negres;
            osstruct.posneg = osstruct.pos & osstruct.negres & ~osstruct.posres;
            osstruct.dspos = osstruct.ds & osstruct.posres & ~osstruct.negres;
            osstruct.dsneg = osstruct.ds & osstruct.negres & ~osstruct.posres;
        end
        
        function map(obj, condition)
            msize = 20;
            cx = obj.load('cx');
            cy = obj.load('cy');
            plot(cy, cx, '.', 'Markersize', msize); hold on
            plot(cy(condition), cx(condition), 'r.', 'Markersize', msize);
        end
        
        function nl_drift = getF1F0Linearity(obj)
            astat = obj.load('allOSstat', 'OSFitResult');
            nNeu = obj.nNeu;
            
            nl_drift = zeros(nNeu, 1);
            
            for i = 1:nNeu
                pref_spat = astat{i}.pref_spat;
                preff1 = astat{i}.fourier(2, pref_spat);
                preff0 = astat{i}.fourier(1, pref_spat);
                sfr = astat{i}.interval_poi_fr;
                nl_drift(i) = preff1 * 2 / (preff0 - sfr);
            end
        end
        
        
        
        function [maxpat, a] = spotSimplePSTH(obj, nID, center, dia, color)
            f = load([obj.datafolder 'Spot_fullpat.mat'], 'patterns');
            realnid = obj.useIDs(nID);
            [maxpat, a] = SpotSimplePSTH(f.patterns, realnid, center, dia, color);
        end
        
        function dstring = getDateString(obj)
            ind = regexp(obj.datafolder, '\d{4}-\d{2}-\d{2}-\d');
            dstring = obj.datafolder(ind:ind+11);
        end
        
        function spotparams = getSpotMLEParameters(obj)
            mlfit = obj.load('MLfit');
            inter = SpotFitInterpretation(mlfit);
            nCol = size(mlfit, 2);
            spotparams = {};
            for i = 1:nCol
                spotparams{i} = StructArrayResolve(inter{i});
            end
            %spotparams{i}
            %sp_b = StructArrayResolve(inter{1});
            %sp_w = StructArrayResolve(inter{2});
            %spotparams = {sp_b, sp_w};
        end
        
        
        % clear cache of the spot MLE
        function clearSpotCache(obj)
            obj.spot_patterns = [];
        end
        
        
        function loadSpotCache(obj, negative)
            if isempty(obj.spot_patterns)
                f = load([obj.datafolder 'Spot_fullpat.mat'], 'patterns');
                if negative
                    f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit_negative');
                    obj.spot_patterns.MLfit = f2.MLfit_negative;
                else
                    f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit');
                    obj.spot_patterns.MLfit = f2.MLfit;
                end
                obj.spot_patterns.patterns = f.patterns;
                %obj.spot_patterns.fullpat = f.patterns.fullpat;
                try
                    obj.spot_patterns.pxfrate = f2.pxfrate;
                catch
                    warning('pxfrate cannot be found. calculating on the fly...')
                    obj.spot_patterns.pxfrate = permute(mean(PatternToCount(f.patterns.fullpat), 4), [2, 3, 1, 5, 4]);
                end
                obj.spot_patterns.negative = negative;
            else
                if obj.spot_patterns.negative ~= negative
                    % if the information is inconsistent, reload the
                    % negative part only.
                    if negative
                        f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit_negative');
                        obj.spot_patterns.MLfit = f2.MLfit_negative;
                        obj.spot_patterns.negative = negative;
                    else
                        f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit');
                        obj.spot_patterns.MLfit = f2.MLfit;
                        obj.spot_patterns.negative = negative;
                    end
                end
            end
        end
            
        function sp_pat = getSpotPattern(obj, nID, negative)
            if nargin < 3
                negative = 0;
            end
            obj.loadSpotCache(negative);

            realnid = obj.useIDs(nID);
            sp_pat = obj.spot_patterns.patterns.fullpat(:,:,:,:,realnid);
        end
            
        
        % show MLE parameters for the spots
        function [spotparams, outsidehist, pcdf] = spotMLE(obj, nID, plotnum, negative, inout, cutradius, onesec)
            % let's use data cache for these results to reduce time for
            % individual plot.
            % if plotnum is not given, it plots all of them.
            % plotnum is organized like this
            % 1: spatial, black
            % 2: temporal, black
            % 3: spatial, white
            % 4: temporal, white
            % 5: all
            % if negative is 1, use MLfit for negative response
            outsidehist = [];
            pcdf = 0;
            if nargin < 3
                plotnum = 5;
            end
            
            if nargin < 4
                negative = 0;
            end
            if nargin < 5
                inout = 0; % default.
            end
            if nargin < 6
                cutradius = 2;
            end
            if nargin < 7
                onesec = 0; % default, show only first 0.5s
            end
                
            % cache data. load only for the first time.
            obj.loadSpotCache(negative);
                
%             
%             if isempty(obj.spot_patterns)
%                 f = load([obj.datafolder 'Spot_fullpat.mat'], 'patterns');
%                 if negative
%                     f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit_negative');
%                     obj.spot_patterns.MLfit = f2.MLfit_negative;
%                 else
%                     f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit');
%                     obj.spot_patterns.MLfit = f2.MLfit;
%                 end
%                 obj.spot_patterns.fullpat = f.patterns.fullpat;
%                 obj.spot_patterns.pxfrate = f2.pxfrate;
%                 obj.spot_patterns.negative = negative;
%             else
%                 if obj.spot_patterns.negative ~= negative
%                     % if the information is inconsistent, reload the
%                     % negative part only.
%                     if negative
%                         f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit_negative');
%                         obj.spot_patterns.MLfit = f2.MLfit_negative;
%                         obj.spot_patterns.negative = negative;
%                     else
%                         f2 = load([obj.datafolder 'Spot_Result.mat'], 'pxfrate', 'MLfit');
%                         obj.spot_patterns.MLfit = f2.MLfit;
%                         obj.spot_patterns.negative = negative;
%                     end
%                 end
%             end

            nCol = size(obj.spot_patterns.pxfrate, 3);

            realnid = obj.useIDs(nID);
            if plotnum == 5;
                for col = 1:nCol
                    subplot(2, nCol, col); cla reset
                    obj.plotSpotpx(obj.spot_patterns.pxfrate,...
                        obj.spot_patterns.MLfit, col, realnid, inout);
                    subplot(2, nCol, col + nCol);cla reset
                    obj.plotSpotTime(obj.spot_patterns.patterns,...
                        obj.spot_patterns.MLfit, col, realnid, inout);
                end
            elseif plotnum == 1
                obj.plotSpotpx(obj.spot_patterns.pxfrate,...
                    obj.spot_patterns.MLfit, 1, realnid, inout, cutradius);
            elseif plotnum == 2
                [outsidehist, pcdf] = obj.plotSpotTime(obj.spot_patterns.patterns,...
                    obj.spot_patterns.MLfit, 1, realnid, inout, cutradius, onesec);
            elseif plotnum == 3
                obj.plotSpotpx(obj.spot_patterns.pxfrate,...
                    obj.spot_patterns.MLfit, 2, realnid, inout, cutradius);
            elseif plotnum == 4
                [outsidehist, pcdf] = obj.plotSpotTime(obj.spot_patterns.patterns,...
                    obj.spot_patterns.MLfit, 2, realnid, inout, cutradius, onesec);
            end
            spotparams = obj.spot_patterns.MLfit(realnid,:);
        end
        
        function plotSpotpx(obj, pxfrate, MLfit, col, realnid, inout, cutradius)
            im = pxfrate(:, :, col, realnid);
            imagesc(im);
            set(gca, 'ydir', 'normal');
            colormap gray
            p = MLfit(realnid, col).params.value;
            el = ellipse(p(4), p(3), p(5), p(2), p(1), 'r');
            el.LineWidth = 1.7;
            if inout
                elcut = ellipse(p(4) * cutradius, p(3) * cutradius, p(5), p(2), p(1), 'y');
                elcut.LineWidth = 1.7;
                elcut.LineStyle = '--';
            end
        end
        function [outsidehist, pcdf] = plotSpotTime(obj, patterns, MLfit, col, realnid, inout, cutradius, onesec)
            if size(MLfit, 2) == 3 % auditory
                auditory = 1;
                onesec = 0;
            else
                if nargin < 8
                    onesec = 1;
                end
                auditory = 0;
            end
            
            %binsize = 0.02; % default
            binsize = 0.05; % for some analysis
            fullpat = patterns.fullpat;
            pattern = fullpat(col, :, :, :, realnid);
            endtime = 0.5;

            
            if onesec % extend the range of patterns.
                ip = patterns.intervalpat(:,realnid);
                stimid = patterns.stimid;
                stimid_p = stimid(col, :, :, :);
                for i = 1:numel(pattern)
                    pattern{i} = [pattern{i} ip{stimid_p(i)}+500];
                end
                endtime = 1;
            end
            
            if auditory
                endtime = 1;
            end
            
            
            catpat = @(x, y) cat(2, pattern{:, y, x, :});
            p = MLfit(realnid, col).params.value;
            outsidehist = []; % returns only when inout is defined.
            pcdf = 0;
            
            
            if inout == 0
                %elist = patternToEventList(fullpat(col, :,:,:,realnid));
                shspot = @(x) stairhist(x, 0:binsize:endtime, 3);

                elist = patternToEventList(pattern);
                shspot(elist(:,5));
                %stairhist(elist(:,5), 0:0.02:0.5,3);
                x = 0:0.001:endtime; % x for fitted curve
                %p = MLfit(realnid, col).params.value;
                pdf = SpotMLEPDF(x, p, endtime);
                hold on
                plot(x, pdf);
                ylabel('PDF (1/s)');
            elseif inout == 1; % inside, in this mode, curve is not plotted
                % as the curve is not plotted, it plots as firing rate
                % 6/28 by the request of reviewer, I also plotted the
                % spontaneous rate.
                %shspot = @(x) stairhist(x, 0:0.02:0.5, 0, 0, 5);
                
                % calculate here spont rate
                icount = squeeze(PatternToCount(ip, 0.2 * 1000));
                spontfr = PoissonFiringRate(sum(icount), length(icount) * 0.3);
                sfr = spontfr.fr
                
                % calculate if there is a significant reduction of the FR
                psize = size(pattern);
                ellip = {p(4) * cutradius, p(3) * cutradius, p(5), p(2), p(1)};
                [catlist, catinside] = PointsOutsideEllipse(ellip, 1:psize(3), 1:psize(2));
                bigcat = [];
                bigcatinside = [];
                for xy = catlist'
                    bigcat = [bigcat catpat(xy(1), xy(2))];
                end
                for xy = catinside'
                    bigcatinside = [bigcatinside catpat(xy(1), xy(2))];
                end
                nPanels_outside = size(catlist, 1);
                nPanels_inside = size(catinside, 1);
                
                normfact_outside = 1/(binsize * psize(4) * nPanels_outside);
                normfact_inside = 1/(binsize * psize(4) * nPanels_inside);
                
                % spontaneous
                plot([0, 1], sfr * [1 1], '--g');
                hold on;
                
                h1 = stairhist(bigcatinside/1000, 0:binsize:endtime, 0, 0, normfact_inside);
                hold on;
                [h2, outsidehist] = stairhist(bigcat/1000, 0:binsize:endtime, 0, 0, normfact_outside);
                
                outsidehist_raw = outsidehist.value / normfact_outside;
                
                % calculate poisson CDF
                minraw = min(outsidehist_raw);
                pcdf = poisscdf(minraw, sfr / normfact_outside);
                
                %h1 = shspot(bigcatinside/1000); hold on;
                %h2 = shspot(bigcat/1000);
                l = legend([h1(1), h2(1)], 'Inside RF', 'Outside RF');
                legend boxoff
                l.Position = l.Position + [0.015 0.035 0 0];
                %ylabel('FR (Hz)');
            end
            %xlabel('Time (s)')
            xlim([0 endtime]);
            box off
        end
        
        
        
        % varargin should be phase, sp_freq, orient for each cell
        function h = CRG_PSTH(obj, nID, varargin)
            f = load([obj.datafolder 'Cphasegratings_fullpat.mat'], 'patterns');
            realnid = obj.useIDs(nID);
            h = YLikeHistogram(f.patterns, realnid, varargin{:});
        end
        
        % get whether the cell is y-like
        function ylikep = getYlikeCells(obj, threshold)
            if nargin < 2
                threshold = 0.01;
            end
            cstat = load([obj.datafolder 'CphaseResult.mat'], 'cstat');
            cstat_needed = cstat.cstat(obj.useIDs);
            
            f1fs = pick(cstat_needed, 'f1fs');
            f1fs = cat(1, f1fs{:}) * 2;
            f2fs = pick(cstat_needed, 'f2fs');
            f2fs = cat(1, f2fs{:}) * 2;
            
            ffdiff = f1fs - f2fs;
            f1valid = ffdiff(:,1).value > 0 & ffdiff(:,1).sig < threshold;
            f2valid = any(ffdiff(:,2:end).value < 0 &...
                ffdiff(:,2:end).sig < threshold / 6, 2); % with bonferroni correction
            ylikep = f1valid & f2valid;
        end
        
        function [clikep, f1f0] = getClikeCells(obj)
            astat = obj.load('allOSstat', 'OSFitResult.mat');
            nNeu = obj.nNeu;
            pref_spat = pick(astat, 'pref_spat');
            preff1 = zeros(nNeu, 1);
            preff0 = zeros(nNeu, 1);
            for i = 1:nNeu
                preff1(i) = astat{i}.fourier(2, pref_spat(i)) * 2;
                preff0(i) = astat{i}.fourier(1, pref_spat(i));
            end
            sfr = pick(astat, 'interval_poi_fr');
            f1f0 = preff1 ./ (preff0 - sfr);
            clikep = abs(f1f0) < 1;
        end
        
        function plotLoomingMUBW(obj)
            lf_elec = obj.load('lf_elec', 'looming');
            emap = obj.loadEmap;
            %subplot(2, 4, d);
            %datanum = d;
            nElec = length(lf_elec);
            maxFR = 0;
            for i = 1:nElec
                %ppar(i) = lf_all_mu{i}.wminusb.value;
                ppar(i) = lf_elec{i}.color_polarity.value;
                %ppar(i) = lf_data{i}.evoked_rate2(2).value;
                %ppar(i) = diffos(lf_data{i}.max_evoked_rate(2).value, lf_data{i}.max_evoked_rate(1).value);
                sv = sum(lf_elec{i}.max_evoked_rate);
                %sv = lf_elec{i}.max_evoked_rate(2); % white
                %sv = sum(lf_data{i}.evoked_rate2);
                spar(i) = sv.value;
                esig = min(lf_elec{i}.response_significance) * 2;
                if esig > 0.01
                    ppar(i) = 500;
                    spar(i) = 0.001;
                end
                if spar(i) > maxFR
                    maxFR = spar(i);
                end
                
            end
            scatter(emap(:,2), emap(:,1), max(spar * 10, 1), ppar, 'filled');
            set(gca, 'CLim', [-0.6 0.6]);
            title(sprintf('maxFR: %3.1f', maxFR));
            colorbar
        end
        
        function emap = loadEmap(obj)
            aID = obj.load('arrayID');
            switch aID
                case 11001
                    load('emap256AN.mat');
                case 10901
                    load('emap256A.mat');
            end
            emap = emap(2:end-1, :);
        end
        
        
        
        function pattern = showPSTH(obj, nID, filename)
            if strcmp(filename, 'looming')
                loomname = getWildName([obj.datafolder 'LoomingSpots_fullpat*.mat']);
                f = load([obj.datafolder loomname], 'patterns', 'sf', 'patternVT');
            else
                f = load([obj.datafolder filename], 'patterns', 'sf', 'patternVT');
            end
            realnid = obj.useIDs(nID);
            %NeuronFullPSTH(f.patterns, f.sf, realnid, 12, f.patternVT.state, 0)
            NeuronFullPSTH(f.patterns, f.sf, realnid, 12, [], 0)
            if nargout > 0
                pattern = f.patterns.fullpat(:,:,:,:,realnid);
            end
                
        end
        
        function flag = isSBC(obj, nID)
%            realnid = obj.useIDs(nID);
            sbc = obj.load('sbcres');
            sbcp = obj.loadp('sbcres');
            if ~iscell(sbc)
                flag = -1; % no data exist. exit.
                return
            end
            
            if sbcp(nID) == 1 % bad measurement
                flag = -2;
                return
            end
            
            sres = sbc{nID};
            
            f1sig = sres.f1sig;
            sbtime = sres.angle / 2 / pi * 10 % converting to seconds
            
            supp = (sbtime < 1 | sbtime > 9);
            stim = (sbtime > 4 & sbtime < 6);
            
            
            if f1sig < 0.01
                if supp
                    flag = 1; % suppressed by contrast
                elseif stim
                    flag = 2; % stimulated by contrast
                else
                    flag = 3; % other
                end
            else
                flag = 0; % nonsignificant response
            end
        end
        
        function npattern = showRasterOTC(obj, nID, spat, withfit, mu)
            if nargin < 3
                spat = 0;
            end
            if spat == 0
                override = 1;
            else
                override = 0;
            end
            if nargin < 5
                mu = 0;
            end
            % default is override. point of giving spat?
            if mu
                load([obj.datafolder 'DriftGratings_fullpat.mat'], 'patterns_mu');
                realnid = nID;
            else
                load([obj.datafolder 'DriftGratings_fullpat.mat'], 'patterns');
                realnid = obj.useIDs(nID);
            end
            

            %OrientationRawData(patterns, realnid, spat);
            
            if nargin < 4
                withfit = 1;
            end
            
            if withfit
                if mu
                    astat = obj.load('allOSstat_mu', 'OSFitResult');
                    p = patterns_mu;
                else
                    astat = obj.load('allOSstat', 'OSFitResult');
                    p = patterns;
                end
                npattern = OrientationRawData(p, realnid, spat, astat{nID}, override);
            else
                npattern = OrientationRawData(p, realnid, spat, [], override);
            end
        end
        
        function driftid = findDriftingGratings(obj)
            load([obj.datafolder 'DriftGratings_fullpat.mat'], 'stimnum');
            driftid = stimnum;
        end
            
        
        function fracrun = fractionTimeRunning(obj, stimname)
            % find fraction of time running
            stimnum = obj.findStimulus(stimname);
            movement = obj.load('movement', 'basicinfo.mat');
            if isempty(movement)
                fracrun = 0;
                return
            end
            stimchoose = 1;
            if strcat(stimname, 'DriftGratings')
                %stimchoose = 2; % for DrittGratings, first one is gray screen
                % stimnum is achieved from actual Drifting gratings used.
                load([obj.datafolder 'DriftGratings_fullpat.mat'], 'stimnum');
                stimchoose = 1;
            end
            v = movement{stimnum(stimchoose)}.vsmooth;
            fracrun = nnz(v > 1) / length(v);
        end
        
        function fracrunc = fractionTimeRunning_cell(obj,stimname)
            fc = obj.fractionTimeRunning(stimname);
            fracrunc = repmat(fc, nnz(obj.useIDs), 1);
        end
        
        function fractionTimeRunning_filter(obj, stimname, lower, upper)
            fc = obj.fractionTimeRunning(stimname);
            if fc < lower | fc > upper
                obj.useIDs = [];
            end
        end
        
        
        
        
        % run analysis on the folder.
        function runAnalysis(obj, stimname, errflag)
            % find the stimulus number
            if nargin < 3
                errflag = 0;
            end
            stimnums = obj.findStimulus(stimname);
            DoAllAnalysis(obj.datafolder, stimnums, errflag);
            
        end
        
        function npe = getNeuronsPerElectrode(obj)
            obj.readBasicInfo;
            load([obj.datafolder 'SurfaceContour.mat'])
            emap = loademap(obj.arrayID); %gives trimmed electrodes.
            
            nElec = size(emap, 1);
            effElec = 0;
            
            for i = 1:nElec
                el = emap(i, :); % el(y, -x);
                sur = cont(mini(abs(cont(:,1) + el(2))), 2);
                if el(1) <= sur
                    effElec = effElec + 1;
                end
            end
            npe = obj.nNeu / effElec;
            %npe = effElec;
        end
        
        function elec_depth = getElectrodeDepth(obj)
            [~, elec_depth] = getDepth(obj.datafolder, 1);
        end
        
        function motionstats = getGratingsLocomotion(obj)
            % check number of neurons
            if isempty(obj.useIDs)
                motionstats = [];
                return 
            end
            load([obj.datafolder, 'DriftGratings_fullpat.mat'], 'patterns', 'patternVT', 'sf')
            load([obj.datafolder, 'OSFitResult.mat'])
            nNeu = length(obj.useIDs);

            
            mistate = patternVT.intervalstate == 1;
            sistate = patternVT.intervalstate == -1;

            for i = 1:nNeu
                n = obj.useIDs(i);
                pat = patterns.fullpat(:,:,:,:,n);
                histn = PatternToCount(pat, 200) / (sf.Duration - 0.2);
                histcell = num2cell(histn);
                
                prefspat = allOSstat{n}.pref_spat;
                
                intcount = PatternToCount(patterns.intervalpat(:,n), 200) / 0.3; % .5s - .2s
                
                %histn = PatternToFRhist(pat, Duration * 1000, Duration * 1000);
                [val_m, val_s, err_m, err_s, nm, ns] = StateDependentPSTH(histcell, patternVT.state);
                
                %this is sort of an ad-hoc thing
                
                
                
                %subplot(2, 1, 1)
                %imagesc(squeeze(val_m))
                %subplot(2, 1, 2)
                %imagesc(squeeze(val_s))
                %combined_image = [squeeze(val_s), squeeze(val_m)];
                %image_storage{n} = combined_image;
                mcurves = nanmean(squeeze(val_m), 2)';
                mcurves_err = rootmeansquare(squeeze(err_m), 2)' ./ sqrt(size(err_m, 3));
                
                motc = squeeze(val_m(:, prefspat, :));
                motc_err = squeeze(err_m(:, prefspat, :));
                
                                
                scurves = nanmean(squeeze(val_s), 2)';
                scurves_err = rootmeansquare(squeeze(err_s), 2)' ./ sqrt(size(err_s, 3));

                sotc = squeeze(val_s(:, prefspat, :));
                sotc_err = squeeze(err_s(:, prefspat, :));

                
                spont_m = nanmean(intcount(mistate));
                spont_m_err = nanstd(intcount(mistate)) / sqrt(nnz(mistate));
                
                spont_s = nanmean(intcount(sistate));
                spont_s_err = nanstd(intcount(sistate)) / sqrt(nnz(sistate));
                
                spont = nanmean(intcount);
                spont_err = nanstd(intcount) / sqrt(length(intcount));
                

                dloc.mcurves = ErrorNum(mcurves, mcurves_err);
                dloc.scurves = ErrorNum(scurves, scurves_err);
                dloc.spont_m = ErrorNum(spont_m, spont_m_err);
                dloc.spont_s = ErrorNum(spont_s, spont_s_err);
                dloc.spont = ErrorNum(spont, spont_err);
                
                dloc.motc = ErrorNum(motc, motc_err);
                dloc.sotc = ErrorNum(sotc, sotc_err);
                
                dloc.pref_m = dloc.mcurves(prefspat);
                dloc.pref_s = dloc.scurves(prefspat);
                

                motionstats(i, 1) = dloc;
            end
        end
        
        
        function fitres = loadAuditoryFitresult(obj, stims)
            % loading all of the results in cell array
            
            if nargin < 2
                stims = 0;
            end
            
            fname = [slashappend(obj.datafolder), 'AuditorySpotSummary_1'];
            if stims
                readvars = {'fitresult_stims', 'posneu1', 'posneu2'};
            else
                readvars = {'fitresult', 'posneu1', 'posneu2'};
            end
            as = load(fname, readvars{:});
            if stims
                as.fitresult = as.fitresult_stims;
            end
            fitsize = size(as.fitresult);
            fitres = cell(obj.nNeuTotal, fitsize(1), fitsize(2));
            
            for i = 1:fitsize(1)
                for j = 1:fitsize(2)
                    if ismember(i, [2 4]) & ~stims
                        fitres(as.posneu2, i, j) = as.fitresult{i, j};
                    else
                        fitres(as.posneu1, i, j) = as.fitresult{i, j};
                    end
                end
            end
            fitres = fitres(obj.useIDs, :, :);
        end
        
        
        % this is only for crude testing
        function h = plotAuditoryMap(obj, mode)
            params = obj.loadGoodAuditoryFit(mode, 0);
            x = obj.load('x', 'xy');
            y = obj.load('y', 'xy');
            
            az = params(:, 4) * 180 / pi;

            good = params(:, 1).value > 0 & az.value < 140 & az.value > 0;
            
            
            yorig = y;
            xorig = x - 1350;
            theta = 25 / 180 * pi;
            rotmat = [cos(theta), sin(theta); -sin(theta), cos(theta)];
            yxcorr = rotmat * [yorig'; xorig'];
            ycorr = yxcorr(1, :)';
            %xcorr = yxcorr(2, :)';
            h = plotEN(ycorr(good), az(good), '.', 'capsize', 0);
        end
            
        
        function params = loadGoodAuditoryFit(obj, mode, aic, stims, fitres)
            % load only good auditory fit.
            % at this point, it does not work when useID is not reset.
            % modes
            % 1: Kent distribution with old parameters, faster time scale
            % 2: Kent distribution with old parameters, slower time scale
            % 3: Kent distribution with new parameters, faster time scale
            % 4: Kent distribution with new parameters, slower time scale
            % 5:
            % 6:
            % 7:
            % 8:
            % 9: New 4-parameter simple model for auditory response
            
            % aic == 2 pass all the data.
            
            if ismember(mode, [9,10])
                nparam = 4;
            else
                nparam = 7;
            end
                        
            if nargin < 4
                stims = 0;
            end
            
            if nargin < 5
                fitres = 0; % if 1, it returns fit ressult.
            end
                
            
            fname = [slashappend(obj.datafolder), 'AuditorySpotSummary'];
            %as = load([slashappend(obj.datafolder), 'AuditorySpotSummary.mat']);
            
            if stims
                readvars = {'fitresult_stims', 'posneu1', 'posneu2'};
            else
                readvars = {'fitresult', 'posneu1', 'posneu2'};
            end
            if exist([fname '.mat'], 'file')
                as = load(fname, readvars{:});
            else
                warning('AuditorySpotSummary not found. reading _1...');
                fname = [fname '_1'];
                as = load(fname, readvars{:});
            end
            
            if stims
                as.fitresult = as.fitresult_stims;
            end
            % let's just do one.
            fits = as.fitresult{mode};
            nFit = length(fits);
            ic_accept = logical(zeros(nFit, 1));
            
            params = ErrorNum.init(obj.nNeuTotal, nparam);
            
            if isempty(fits)
                return
            end
            
            if ismember(mode, [1, 3, 5, 6, 7, 8, 9, 10]) | stims
                posneuind = (as.posneu1);
            elseif ismember(mode, [2, 4])
                posneuind = (as.posneu2);
            end
            
            %fits
            for i = 1:nFit
                if aic == 2 % passing mode. accept all
                    ic_accept(i) = 1;
                    continue;
                end
                if fits{i}.goodfit_err
                    r = fits{i};
                    if aic
                        factor = 2;
                    else
                        factor = log(85);
                    end
                    ic_model = r.nll2 + factor * nparam;
                    ic_flat = r.flatnll2 + factor * 1;
                    if ic_flat > ic_model
                        ic_accept(i) = 1;
                    end
                end
            end
            % put values into full EN array
            ic_ac_ind = find(ic_accept);
            acceptid = posneuind(ic_accept)';
            
            for i = 1:length(acceptid)
                r = fits{ic_ac_ind(i)};
                p = r.params_err;
                e = sqrt(diag(r.cov_err))';
                if ismember(mode, [1, 2]) & ~stims
                    if p(4) > pi / 2
                        p(3) = -p(3); % flip azimuth for large theta
                    end
                end
                params(acceptid(i), :) = ErrorNum(p, e);
            end
            params = params(obj.useIDs, :);
        end


        
        % this function already returns the ASDF in TTL time.
        % this function could be problematic if TTL are not sent very
        % frequently. (lower than 1Hz), as it only reads data 1s after the
        % last TTL.
        % Also, needs rework with intanSyncPulse.
        function asdf_s = segmentASDF(obj, datanumber)
            % returns asdf in the range segment ttl is available + 1s.
            obj.readBasicInfo
            load([obj.datafolder, 'asdf.mat'])
            asdf_s = ASDFChooseTime(asdf_raw,...
                obj.segttls{datanumber}(1),...
                obj.segttls{datanumber}(end) + 1000);
            
            
            asdf_s = ASDFSubsample(asdf_s, obj.useIDs);
        end
        
        function plotVT(obj, dataname)
            load([obj.datafolder, dataname], 'patternVT', 'sf');
            figure(501);
            NeuronFullPSTH(patternVT, sf, 1);
            figure(502);
            NeuronFullPSTH(patternVT, sf, 2);
        end
    end
end