function ext = get_extension(filename)
%   GET_EXTENSION   Gets the extension to a filename
%     [EXT] = GET_EXTENSION(FILENAME)
% 
%   Basically uses fileparts
%   
%   Created by Zarrar Shehzad on 2012-10-01.

[~,~,ext] = fileparts(filename);

end %  function