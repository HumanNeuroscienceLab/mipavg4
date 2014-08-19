function recodes = read_recode(RCDfiles, logg)
%   READ_RECODE   Read codes from file allowing one to recode events
%     [RECODES] = READ_RECODE(RCDFILES, LOGG)
% 
%   Input:
%       RCDfiles    - cell array of filenames and in each file there should
%                     be two columns, one for the original event codes and 
%                     another for the new event codes
%       logg        - instance of logger object
%
%   Output:
%       recodes     - cell array of length RCDfiles with each element being a
%                     2 column matrix of [orig-code new-code]
%   
%   Created by Zarrar Shehzad on 2012-09-12.

    nRCDfiles   = length(RCDfiles);
    recodes     = cell(1, nRCDfiles);
    for i=1:nRCDfiles
        recodes{i} = dlmread(RCDfiles{i});
        logg.write('\nRecode file specified: %s\n', RCDfiles{i});
        logg.write('Number of recodes specified = %d\n', size(recodes{i},1));
    end
    
end %  function