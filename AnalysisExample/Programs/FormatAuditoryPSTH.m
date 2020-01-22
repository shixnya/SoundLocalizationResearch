function patterns = FormatAuditoryPSTH(asdf, ttl, stiminfo, duration, old)
if nargin < 5
    old = 0;% for compatibility
end


psthcell = fastttldivide(asdf, ttl, duration);
nNeu = asdf{end}(1);

if old
    nelev = stiminfo.nelev;
    nazim = stiminfo.nazim;
    ntrial = stiminfo.ntrial; % Yes. You need to include this.
    ninten = length(ttl) / nelev / nazim / ntrial; % number of intensities
    
    % get trial based stiminformation
    snumt = reshape(stiminfo.stiminfo.stimnums, [nelev * nazim, ntrial]) + 1;
    
    % currently, it is linear except for the neuron ID.
    psthmat_base = reshape(psthcell, [nelev * nazim, ntrial, ninten, nNeu]);
    psthmat_stimcorrected = cell(nelev * nazim, ntrial, ninten, nNeu);
    for i = 1:ntrial
        psthmat_stimcorrected(snumt(:,i), i, :, :) = psthmat_base(:, i, :, :);
    end
    %psthmat_rightdim = reshape(psthmat_stimcorrected, [nelev, nazim, ntrial, ninten, nneu]);
    %psthmat = permute(psthmat_rightdim, [4, 1, 2, 3, 5]);
    
    psthmat_rightdim = reshape(psthmat_stimcorrected, [nazim, nelev, ntrial, ninten, nNeu]);
    psthmat = permute(psthmat_rightdim, [4, 2, 1, 3, 5]);
    
    %psthmat = cell(ninten, nelev, nazim, ntrial, nneu); % match with old spots
else % new stimulus interpretation
    nmeta1 = length(stiminfo.meta1);
    nmeta2 = length(stiminfo.meta2);
    nreps = stiminfo.nreps;
    serial = stiminfo.serial;
    npattern = stiminfo.npattern;
    if isfield(stiminfo, 'ntypes')
        ntypes = stiminfo.ntypes;
    else
        ntypes = 1;
    end
    % keep track on the previous stimuluation number
    prevstim = [0 stiminfo.stimnums(1:end-1) + 1];
    %nstim = length(ttl) / npattern / nreps;
    nstim = 1; % we don't do multiple stim in one segment any more.
    if length(ttl) < npattern * nreps % incomplete stimulus
        nreps_new = floor(length(ttl) / npattern);
        nthrow = mod(length(ttl), npattern);
        answer = questdlg(...
            sprintf(['Number of TTL is less than expected.\n'...
            'It has only full %d repetitions, and %d TTL pulses'...
            'will be thrown away.\nDo you want to proceed with this?'],...
            nreps_new, nthrow));
        if strcmp(answer, 'Yes')
            psthcell = psthcell(1:(npattern * nreps_new), :);
            nreps = nreps_new;
        else
            error('You chose to stop the process.')
        end
    end
        
    
    % add one to convert zero-based (python) to one-based (matlab)
    snumt = reshape(stiminfo.stimnums(1:(npattern * nreps)), [npattern, nreps]) + 1;
    psthmat_base = reshape(psthcell, [npattern, nreps, nstim, nNeu]);
    psthmat_stimcorrected = cell(npattern, nreps, nstim, nNeu);
    
    prevstim_base = reshape(prevstim, [npattern, nreps, nstim]);
    prevstim_stimcorrected = zeros(npattern, nreps, nstim);
    for i = 1:nreps
        psthmat_stimcorrected(snumt(:,i), i, :, :) = psthmat_base(:, i, :, :);
        prevstim_stimcorrected(snumt(:,i), i, :, :) = prevstim_base(:, i, :, :);
    end
    if serial
        psthmat_rightdim = reshape(psthmat_stimcorrected, [npattern, 1, 1, nreps, nstim, nNeu]);
        prevstim_rightdim = reshape(prevstim_stimcorrected, [npattern, 1, 1, nreps, nstim]);
    else
        psthmat_rightdim = reshape(psthmat_stimcorrected, [nmeta2, nmeta1, ntypes, nreps, nstim, nNeu]);
        prevstim_rightdim = reshape(prevstim_stimcorrected, [nmeta2, nmeta1, ntypes, nreps, nstim]);
    end
    psthmat = permute(psthmat_rightdim, [3, 2, 1, 4, 6, 5]); % pushing nstim to the end...
    prevstim = permute(prevstim_rightdim, [3, 2, 1, 4, 6, 5]);
end

patterns.fullpat = psthmat;
patterns.prevstim = prevstim;
