function EEG = fieldtrip2econnectome(timelock, name, Fs)
%     FIELDTRIP2ECONNECTOME   Converts average time-locked object to econnectome
%         [EEG] = FIELDTRIP2ECONNECTOME(TIMELOCK, NAME, FS)
% 
%     Inputs:
%       timelock    - output from ft_timelockanalysis
%       name        - name for given trials
%       Fs          - frequency that data was sampled
%     
%     Outputs:
%       EEG         - econnectome data structure
%     
%     Created by Zarrar Shehzad on 2012-09-12.

    EEG = [];
    EEG.name        = name;
    EEG.type        = 'EEG';
    EEG.nbchan      = length(timelock.label);
    EEG.points      = length(timelock.time);
    EEG.srate       = Fs;
    EEG.labeltype   = 'standard';
    EEG.labels = {};
    for i = 1:length(timelock.label)
       EEG.labels{i} = timelock.label{i};
    end
    EEG.data = timelock.avg;
    EEG.unit = 'uV';

end %  function