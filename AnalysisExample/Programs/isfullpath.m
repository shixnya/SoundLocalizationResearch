function isfull = isfullpath(dirname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ispc
    isfull = dirname(2) == ':';
else
    isfull = dirname(1) == filesep;
end


end

