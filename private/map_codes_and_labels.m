function ALLEEG = map_codes_and_labels(ALLEEG, codes, labels, epoch, logg)
%   [ALLEEG] = map_codes_and_labels(ALLEEG, codes, labels, epoch, logg)
% 
%   Replaces event type values based on mapping in codes, which results in
%   changing [ALLEEG(i).event.type]. Corresponding labels are also added
%   to [ALLEEG(i).urevent.label] and [ALLEEG(i).event.label].
%
%   Inputs:
%       ALLEEG      - struct array of EEG data
%                     contains EEG.event with trial information
%       codes       - a cell array where each element is a 2-column matrix
%                     1st column = input codes, corresponding to [EEG.event.type]
%                     2nd column = output codes, replacing input values
%       labels      - similar to codes, cell array with each element a 2-column matrix
%                     1st column = input labels
%                     2nd column = output labels
%       epoch       - 2 element vector with start & end time of trial in seconds
%                     this is only used to check that all codes have complete epochs
%       logg        - instance of logger object
%   
%   Outputs:
%       ALLEEG      - struct array of EEG data
%
%   Created by Zarrar Shehzad on 2012-09-06

ALLcodes = check_against_ALLEEG(codes);
ALLlabels = check_against_ALLEEG(labels);
[ALLEEG.eventcodes] = deal([]); [ALLEEG.eventlabels] = deal({});
clear codes labels;
n = length(ALLEEG);

for i=1:n
    logg.write('\nProcessing EEG run #%i\n', i);
    
    EEG     = ALLEEG(i);
    codes   = ALLcodes{i};
    labels  = ALLlabels{i};
    
    if size(codes,1) ~= size(labels,1)
       logg.error('Size mismatch between codes and labels\n');
       error('Program cannot proceed');
    end
    
    raw_codes = [EEG.event.type];
    [EEG.event.type] = deal(0);
    
    %% Map codes and labels
    logg.write('\nMapping codes and labels (displaying info for inputs)\n');
    % Add label field
    [EEG.urevent.label] = deal('');
    [EEG.event.label]   = deal('');
    for j=1:size(codes, 1)
        % Find input codes in EEGLAB event data structure
        search = raw_codes == codes(j,1);
        if any(search)
            logg.write('\tFound %d events for %s trials (code %d)\n', ...
                        sum(search), labels{j,1}, codes(j,1));
            % Add input labels to original event structure
            [EEG.urevent([EEG.event(search).urevent]).label] = deal(labels{j,1});
            % Change event codes (ie type) based on user specified output codes
            [EEG.event(search).type] = deal(codes(j,2));
            % Add output labels to event structure
            [EEG.event(search).label] = deal(labels{j,2});
        else
            logg.write('\tNO events found for %s trials (code %d)\n', ...
                        labels{j,1}, codes(j,1));
        end
    end
    
    %% Store unique codes and labels
    logg.write('\nStoring unique codes with labels\n');
    [~, IC]  = unique(codes(:,2));
    [~, IL] = unique(labels(:,2));
    IC = sort(IC); IL = sort(IL);
    EEG.eventcodes  = codes(IC,2);
    EEG.eventlabels = labels(IL,2);
    
    %% Check that everything is good
    logg.write('\nChecking codes and labels\n');
    
    % Ensure that unique codes and labels match
    ncodes  = length(EEG.eventcodes);
    nlabels = length(EEG.eventlabels);
    if ncodes ~= nlabels
        logg.error(['# of unique output codes (%d) not equal to # of ' ...
                    'unique output labels (%d). Check the 3rd and 4th ' ...
                    'columns of your SGC file.'], ncodes, nlabels);
        error('Program cannot proceed');
    end
    if any(IC ~= IL)
        logg.error(['Mismatch between output codes and labels. Check the ' ...
                    '3rd and 4th columns of your SGC file.']);
        error('Program cannot proceed');
    end
    
    % Check that at least one output code/label exists
    for j=1:ncodes
        if ~any([EEG.event.type] == EEG.eventcodes(j))
            logg.write('NOTE: No trials were found with output event code %d',  ...
                        EEG.eventcodes(j));
        end
        if ~any(ismember({EEG.event.label}, EEG.eventlabels{j}))
            logg.write('NOTE: No trials were found with output event label %d', ...
                        EEG.eventlabels{j});
        end
    end
    
    % Make sure that all codes have a complete epoch and do not extend past 
    % either end of data array. If so, remove event.
    epoch_pts = epoch .* EEG.srate;
    badevents = [EEG.event.latency] + epoch_pts(1) < 0 | ...
                [EEG.event.latency] + epoch_pts(2) > EEG.pnts;
    if any(badevents)
        logg.write('');
        for e = find(badevents)
            logg.write(['NOTE: For trial #%d, removing event code %d / ' ...
                        'label %s with latency of %d, event epoch exceeds data' ...
                        ' limits.'], e, EEG.event(e).type, EEG.event(e).label,  ...
                        EEG.event(e).latency);
        end
        EEG.event = EEG.event(~badevents);
    end
    
    %% Report unique codes/labels
    logg.write('\nNumber of events based on output codes:\n');
    for j=1:ncodes
        n = sum([EEG.event.type] == EEG.eventcodes(j));
        logg.write('\t%d for code %d and label %s\n', n, ...
                   EEG.eventcodes(j), EEG.eventlabels{j});
    end
        
    ALLEEG(i) = EEG;
end

% Check that at least one output code/label exists ACROSS RUNS
[ucodes, I]  = unique([ALLEEG.eventcodes]);
ulabels = [ALLEEG.eventlabels]; ulabels = ulabels(I);
allevents = [ALLEEG.event];
allcodes = [allevents.type];
alllabels = {allevents.label};
for j=1:length(ucodes)
    if ~any(allcodes == ucodes(j))
        logg.error('No trials were found with output event code %d',  ...
                    ucodes(j));
        error('Program cannot proceed');
    end
    if ~any(ismember(alllabels, ulabels{j}))
        logg.error('No trials were found with output event label %d', ...
                    ulabels{j});
        error('Program cannot proceed');
    end
end


function variable = check_against_ALLEEG(variable)
%   [variable] = check_against_ALLEEG(variable)
%
%   Checks that variable length matches that of ALLEEG.
%   Will repeat variable if its length is 1 in order to match its length with ALLEEG.
%
    vname=@(x) inputname(1);
    nVAR = length(variable); nEEG = length(ALLEEG);
    if nVAR ~= nEEG;
        if nVAR == 1
          variable = repmat(variable, 1, nEEG); 
        else
          logg.error('Length of %s must match ALLEEG or be 1\n', vname(variable));
          error('Program cannot continue');
        end
    end
end

end %  function