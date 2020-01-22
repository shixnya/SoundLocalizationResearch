function res = MinuitFunctionFit(func, x, y, initp, err, extrad, bound1, bound2, mode)
% function res = MinuitFunctionFit(func, x, y, initp, err, extrad, bound1, bound2, mode)
%
%    func - Function handle func(p, x)
%    x - (any) 2nd argument of func
%    y - (same as x) data to fit to
%    initp - (any) initial parameters. 1st artment of func
%    err - (same as y) error of the data
%    extrad - (1, 1) extra degree of freedom to subtract
%    bound1 - lower bound of the parameters
%    bound2 - upper bound of the parameters
%    mode - fit mode. 0: poisson likelihood (default), 1: gaussian chi2
%
% Returns:
%    res.goodfit_err: whether the fit was good or not.
%    res.goodfit_noweight
%
% Description:
%    do function fit using fminuit.
%    save results with statistics.
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
% Created: during 2017 paper
% Modified: 9/24/2018: this help is added



if nargin < 6
    extrad = 0;
end
if nargin < 7
    bound1 = false;
    bound2 = false;
end
if nargin < 9
    mode = 0; % default is poisson. gaussian if 1, 
end

goodind = find(err);
res.dof = length(goodind) - length(initp) - extrad;

if res.dof <= 0
    warning('Degree of freedomj <= 0');
    res.goodfit_noweight = 0;
    res.goodfit_err = 0;
    return
end

if sum(y) == 0
    warning('No spikes found. Quitting.')
    res.goodfit_noweight = 0;
    res.goodfit_err = 0;
    return
end



passdata.func = func;
passdata.x = x;
passdata.y_obs = y;
passdata.y_err = err;
if mode == 0
    passdata.minfunc = @minfunc_2nllf_poisson;
elseif mode == 1
    passdata.minfunc = @minfunc_chi2;
else
    error('Please spedcify valid fit mode (0: poisson, 1:chi2)');
end

mincom = 'set pri -1; scan; minimize; simplex; improve; minos';

% params, p_err, nll2, err_matrix

if any(bound1)
    % fit with bounded parameters
    stepbounds = [1:length(initp);bound1;bound2]';
    stepbounds(isinf(stepbounds)) = 0;% bound1 = bound2 = 0 is considered unbounded. (MINUIT manual page 14)
    
    [params, params_err, nll2, cov_err] = fminuit('minuit_passfunc', initp,...
        'b', '-c', mincom, '-s', stepbounds, passdata);
else
    [params, params_err, nll2, cov_err] = fminuit('minuit_passfunc', initp,...
        'b', '-c', mincom, passdata);
end


if any(isnan(params)) || isinf(nll2) || isnan(nll2) || any(any(isnan(cov_err))) || any(diag(cov_err) < 0)
    res.goodfit_err = 0;
    res.goodfit_noweight = 0;
else
    res.params_err = params;
    res.cov_err = cov_err;
    
    if isnumeric(x)
        res.chi2_err = sum((func(params, x(goodind)) - y(goodind)).^2 ./ err(goodind).^2) / res.dof;
        res.chi2sig_err = redchi2sig(res.chi2_err, res.dof);
    else
        res.chi2_err = nan;
        res.chi2sig_err = nan;
    end
    
    res.passdata = passdata;
    res.nll2 = nll2;
    res.nll2red = nll2 / res.dof;
    res.nll2sig = redchi2sig(res.nll2red, res.dof);
    
    res.goodfit_err = 1;
    res.goodfit_noweight = 0; % just not doing it.
end