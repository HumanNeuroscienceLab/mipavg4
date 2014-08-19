function ALLEEG = recode_eeg(ALLEEG, ALLrecodes, logg)
%   RECODE_EEG   Changes the event codes based on recode matrix
%     [ALLEEG] = RECODE_EEG(ALLEEG, ALLRECODES, LOGG)
% 
%   Inputs:
%       ALLEEG      - struct array of EEGLab data structures
%       ALLrecodes  - cell array of recode matrices, each matrix should be 
%                     2 columns and have as many rows as there are events in
%                     the associated EEG object
%       logg        - instance of logger object
%   
%   Created by Zarrar Shehzad on 2012-09-12.

n = length(ALLEEG);

if n ~= length(ALLrecodes)
    error(['\nMismatch between number of EEG files (%d) and recode files (%d)\n' ...
           'Program cannot proceed'], n, length(ALLrecodes));
end

for i=1:n
    EEG = ALLEEG(i); recodes = ALLrecodes{i};
    
    norig = length(EEG.event); nrecodes = size(recodes, 1);
    if norig ~= nrecodes
        logg.error('Mismatch between number of original and recodes: %d %d\n', norig, nrecodes);
        error('Program cannot proceed');
    end
    
    if ischar([EEG.event.type])
        logg.error('The input EEG.event must be numeric')
        error('Program cannot proceed')
    end
    
    logg.write('\nSeq\t Old-EEG\t Old-Recode\t New-Recode\n');
    
    for j=1:norig
        logg.write('%4d\t %4d\t %4d\t %4d\n', j, EEG.event(j).type, ....
                    recodes(j,1), recodes(j,2));
        if EEG.event(j).type ~= recodes(j,1)
            logg.error('Mismatch in EEG event code and user recode\n');
            error('Program cannot proceed');
        else
            EEG.event(j).type = recodes(j,2);
        end
    end
    
    ALLEEG(i) = EEG;
end

end %  function