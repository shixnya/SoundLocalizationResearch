function res = vonMisesFit1DEN(theta, yen, initp)
% function res = vonMisesFit1DEN(theta, yen, initp)
%
%    theta - (n,1) vector that contains angles
%    yen - (n,1) ErrorNum vector of y-value
%    initp - (3, 1) initial parameters for von Mises.
%                   [amplitude, mu, kappa];
%
% Returns:
%    res: fit result from MinuitFunctionFit.
%
% Description:
%    Performs minuit function fit with von Mises distribution.
%    If you don't give initp, it will estimate them automatically.
%
% Example:
%    
%
% Requires:
%    ErrorNum.m
%    MinuitFunctionFit.m (including fminuit)
%    
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/24/2018
% Modified: 

dlen = length(theta);

if nargin < 3
    % kappa is diffficult to estimate, so assume 1.
    initp = [sum(yen.value) / dlen * 2 * pi, theta(mini(-yen.value)), 1];
end

res = MinuitFunctionFit(@vonMises, theta, yen.value, initp, yen.err);
