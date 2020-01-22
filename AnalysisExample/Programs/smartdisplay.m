function [resstring, sigstring] = smartdisplay(value, err, precise)
% a function to make a string of values in a format
% value +/- err, using the leading term of the error.
if nargin < 3
    precise = 0;
end

digit = floor(log10(err));
if precise
    digit = digit - 1;
end

if isinf(digit) % means zero error. special case
    resstring = sprintf(['%s ' char(177) ' 0'], num2str(value));
elseif digit < 0 % floating number, should be most cases
    format = sprintf('%%.%df', abs(digit));
    resstring = sprintf([format ' ' char(177) ' ' format], value, err);
else
    format = '%.0f';
    resstring = sprintf([format ' ' char(177) ' ' format], value, err);
end
    
if nargout > 1
    sigval = sigmaToPval(value, err, 0);
    if sigval < 0.0001;
        sigstring = 'p < 0.0001';
    else
        lord = abs(floor(log10(sigval)));
        lord = max(lord, 2); % ensure 2 digits.

        format = sprintf('%%.%df', lord);
        sigstring = sprintf(['p = ' format], sigval);
    end
end