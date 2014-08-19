function check_output_file(filename, logg)
% Logs the existence of given output file
%   check_output_file(filename, logg)
%   
%   Created by Zarrar Shehzad on 2012-09-08.
%

if ~exist(filename, 'file')
    logg.write('...specified file ''%s'' does not currently exist: Creating file...\n', filename);
else
    logg.write('...specified file ''%s'' already exists: Overwriting old file...\n', filename);
end

end %  function
