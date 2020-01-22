function res = SpotFit_onecell(ocpattern, sf, mode)
% spot fit for one cell

if nargin < 3
    mode = 1; % mode 1: normalized by infinite gaussian
    % mode 3: normalize by bounded gaussian (use it for negative response fit)
    % mode 4: see code below (I didn't understand well.)
    % mode 5: Auditory spot (1s duration)
end

% spot fit function for saving the data later and visualization

for i = 1:size(ocpattern, 1) % color or intensity (auditory)
    ev_list = patternToEventList(ocpattern(i,:,:,:)); % ms -> s
    fulllist_one = ev_list(:, [2 3 5]);
    
    nSpike = size(fulllist_one,1);
    res(1,i).nSpike = nSpike;
    
    if nSpike < 20
        warning('Spike number too small');
        res(1,i).params = ErrorNum(nan(1,10), nan(1,10));
        res(1,i).nllf = nan;
        res(1,i).cov = nan(10);
        continue
    end
    
    maxx = sf.nSteps1;
    maxy = sf.nSteps2;
    
    inputdata.data = fulllist_one;
    inputdata.maxx = maxx;
    inputdata.maxy = maxy;
    inputdata.mode = mode;
    
    % do both calculation fallback
    means_one = mean(fulllist_one);
    stds_one = std(fulllist_one);
    
    % means_one: mean of x, y, time
    % stds_one: std of x, y, time
    initp = [means_one stds_one, 0, 0.1];
    bound1 = [1 1 0.02 0.5 0.5 0.01 -2*pi 0];
    bound2 = [maxy maxx 0.4 10 10 0.25 2*pi 1];
    
    stepsize = abs(initp / 100);
    stepsize(7) = 0.01;
    stepbounds = [1:length(initp);stepsize;bound1;bound2]';

    % x y sx sy? or y x sy sx? 4/5
    
    % initp is [x, y, t, sx, sy, st, ang, frac]
    % initp2 is [x, y, sx, sy, ang, t, st1, st2, x????, frac]
    order = [1 2 4 5 7 3 6 6 1 8]; % 9th parameter is dummy.
    initp2 = initp(order);
    initp2(9) = 0.1; % right side bias
    initp2(10) = 0.9; % signal component
    
    stepbounds2 = stepbounds(order, :);
    stepbounds2(9, :) = [9, initp2(9) / 100, 0 10];
    stepbounds2(:, 1) = 1:10; % make sure param ID is set
    if ismember(mode, [1])
        stepbounds2(10, 3) = 0;
    else
        stepbounds2(10, 3) = -1; % negative value allows negative response.
    end
    stepbounds2(10, 4) = 1.0;
    
    
    if ismember(mode, [5, 6, 7])
        [initp2, stepbounds2] = AuditorySetInits(inputdata);
        %initp2(6) = 0.3; % mean of time is larger for auditory
        %initp2([7, 8]) = 0.3; % std time should also be larger
    else
        initp2(6) = 0.1;
    end
    
    if mode ~=4
        if mode == 6
            stepbounds2(1:5, 2) = 0; % turning off modulation of spatial part
        elseif mode == 7
            initp2(1) = 0;
            stepbounds2(1, 2) = 0.01; % the first argument becomes ILD offset
            
        end
        [o1, o2, o3, o4] = fminuit('minfunc_spot_asymtime', initp2, 'b', '-c',...
            'set pri -1;seek;scan;minimize;improve;minos', '-s', stepbounds2, inputdata);
        %minc = @(x) minfunc_spot_asymtime(x, inputdata);
        %[o1, o3, ~, ~, ~, ~, o4] = fmincon(minc,...
        %    initp2, [], [], [], [], stepbounds2(:, 3), stepbounds2(:, 4));
        %o4 = inv(o4 / 2);
        %o2 = sqrt(diag(o4));
    else
        % also play with the input
        counts = sum(PatternToCount(ocpattern(i,:,:,:)), 4); % getting average response
        counts = squeeze(counts);
        
        
        [~, minpos] = minmd(counts);
        initp3 = [minpos 2, 2, 0, 0.1 0.05 0.05, 0.0, -0.01];
        
        % for this fit, I set the minimum radius to be 1.
        stepbounds2([3, 4], 3) = 1;
        stepbounds2(9, 2) = 0; % prohibit baseline movement.
        stepbounds2([3, 4], 4) = 5;
        
        
        
        [o1, o2, o3, o4] = fminuit('minfunc_spot_asymtime', initp3, 'b', '-c',...
            'set pri -1;scan 1; scan 2; scan 1; minimize;improve;hess', '-s', stepbounds2, inputdata);
    end
    
    res(1,i).params = ErrorNum(o1, o2);
    res(1,i).nllf = o3;
    res(1,i).cov = o4;
end





