function pathstr = dirname(file)
    pathstr = fileparts(file);
    if isempty(pathstr), pathstr = '.'; end
end
