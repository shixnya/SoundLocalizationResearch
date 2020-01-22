function [freq, strength] = fftstrength(signal, fs)
% function [freq, strength] = fftstrength(signal, fs)
% a function to produce a plottable fft results.

nsamp = length(signal);
duration = nsamp / fs;
fincr = 1/duration;

freq = (1:nsamp/2) * fincr;

ffts = fft(signal);
strength = abs(ffts(2:end/2+1, :));

