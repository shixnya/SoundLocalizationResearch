function brokenpath = readpath(fpath)


n = 1;
brokenpath = {''};
for i = 1:length(fpath)
    if fpath(i) == '/' | fpath(i) == '\'
        brokenpath{n} = [brokenpath{n} fpath(i)];
        n = n + 1;
        brokenpath{n} = '';
    else
        brokenpath{n} = [brokenpath{n} fpath(i)];
    end
end

