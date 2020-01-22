function [expandedname] = getWildName(wildname)
%GETWILDNAME Summary of this function goes here
%   Detailed explanation goes here
a = dir(wildname);
if isempty(a)
    error('No file matches with the given pattern: %s', wildname)
else
    expandedname = a.name;
end

