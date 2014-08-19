function filename = basename(file)
    [~, name, ext] = fileparts(file);
    filename = [name ext];
end