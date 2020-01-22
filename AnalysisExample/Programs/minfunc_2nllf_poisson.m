function ret = minfunc_2nllf_poisson(y_obs, y_err, y_model)
% Negative log likelihood function for poisson statistics
% primarily used for minuit likelihood fit.


if any(y_model < 0) % punishing negative FR
    ret = inf;
    return
end

if any(y_err == 0)
    error('There are zeros in the y_err. Please replace it with the minimum allowed value of error.');
end

effective_duration = y_obs ./ y_err.^2; % (Hz/Hz^2 = s)

e0 = effective_duration == 0;
effective_duration(e0) = 1 ./ y_err(e0);

yc_obs = y_obs .* effective_duration; % return to count statistics
yc_model = y_model .* effective_duration;
yc_model(isnan(yc_obs)) = 0;

second = yc_obs .* log((yc_model)); % second term
second(yc_obs==0) = 0; % overwrite NaN produced by log(0).
second(isnan(yc_obs)) = 0;

% including third term (this is not necessary for minimization, but helpful
% for significance test by making it likelihood ratio.
third = yc_obs .* (log(yc_obs)-1);
third(yc_obs==0) = 0; % again, overwrite NaN with zero.
third(isnan(yc_obs)) = 0;

ret = 2 * sum(sum(((yc_model) - second + third)));
% this factor of 2 is due to adjusting differences between chi2 and
% likelihood statistics. Although it does not change the fit parameters, it
% adjusts the error matrix.
% One can remove it and set UP=0.5 instead of doing this. (refer to MINUIT
% manual ver 94.1 page 39.)
