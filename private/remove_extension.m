function filename = remove_extension(filename)
%   REMPOVE_EXTENSION   Removes the extension to a filename
%     [FILENAME] = REMOVE_EXTENSION(FILENAME)
% 
%   Basically uses fileparts.
%   
%   Created by Zarrar Shehzad on 2012-10-01.

[pathstr,name,~] = fileparts(filename);
if isempty(pathstr), pathstr = '.'; end
filename = [pathstr filesep name];

end %  function
