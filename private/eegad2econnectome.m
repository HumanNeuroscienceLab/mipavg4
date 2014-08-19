function EEG = fieldtrip2econnectome(hdr, avg)
%     EEGAD2ECONNECTOME   Converts average ERP to econnectome
%         [EEG] = eegad2econnectome(avg)
% 
%     Inputs:
%       hdr         - not exactly a header but an EEGAD friendly object 
%                     returned by EEGRead2
%       avg         - average ERP
%     
%     Outputs:
%       EEG         - econnectome data structure
%     
%     Created by Zarrar Shehzad on 2012-10-04.

    EEG = [];
    EEG.name        = hdr.expName;
    EEG.type        = 'EEG';
    EEG.nbchan      = hdr.nChannels;
    EEG.points      = hdr.nPoints;
    EEG.srate       = 1000./hdr.sampling;
    EEG.labeltype   = 'standard';
    EEG.labels = {};
    for i = 1:hdr.nChannels
       EEG.labels{i} = hdr.chanNames{i};
    end
    EEG.data = avg;
    EEG.unit = 'uV';

end %  function
