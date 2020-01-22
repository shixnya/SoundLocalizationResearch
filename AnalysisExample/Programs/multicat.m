function strings = multicat(string1, string2)
% function strings = multicat(string, cellofstrings)
% concatnate string with cells of strings and return cells of string

if iscellstr(string2)
    strings = cell(size(string2));
    for i = 1:numel(string2)
        strings{i} = [string1 string2{i}];
    end
elseif iscellstr(string1) % the opposite operation
    strings = cell(size(string1));
    for i = 1:numel(string1)
        strings{i} = [string1{i} string2];
    end
end
