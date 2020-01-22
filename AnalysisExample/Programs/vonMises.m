function y = vonMises(p, x)
% function val = vonMises(p, x);
%
%    p - (3,1) Amplitude (A), center (mu), concentration (kappa) of von
%    Mises distribution
%    x - (vector) x, angle
%
% Returns:
%    y - (vector) same dimension as x, value of von Mises distribution
%
% Description:
%    An implementation of von Mises distribution.
%    von Mises distribution is a circular version of Gaussian.
%    It takes two arguments and 
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
% Created: 9/21/2018
% Modified: 9/21/2018

y = p(1) * exp(p(3) * cos(x - p(2))) / (2 * pi * besseli(0, p(3)));
