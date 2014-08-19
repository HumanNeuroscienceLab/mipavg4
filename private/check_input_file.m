function check_input_file(filename, logg)
%   CHECK_INPUT_FILE   Checks if the filename exists
%     check_input_file(filename, logg)
% 
%   Checks if the filename exists, otherwise throws an error
%   
%   Created by Zarrar Shehzad on 2012-10-03.

if ~exist(filename, 'file')
    logg.error('Input filename %s could not be found\n', filename);
    error('Input file not found. Program cannot proceed');
end

end %  function