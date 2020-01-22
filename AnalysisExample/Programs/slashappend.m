function slashed = slashappend(somepath)

if somepath(end) ~= filesep
    slashed = [somepath filesep];
else
    slashed = somepath;
end