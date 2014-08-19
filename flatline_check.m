function [r, rvalue] = flat3_check(single_trial, npts, channel, wSize, wProgression, wEnvelope)

% initialize reject flag and reject value variable
r = 0;
rvalue = 0;

% create a window starting from 1 to window size
start_val = 1;
end_val = wSize;

%set default if no custom envelope
if wEnvelope ==0
    wEnvelope = 80;   %use 1000 for non-converted data
end;    



% loop through entire trial, moving the window by wProgression points

for n=1:((npts/wProgression)-(wSize/wProgression - 1)) % calculates the # of times the window moves
    
    % find the min and max within the window
    min_of_trial = min(single_trial(channel, start_val:end_val)');
    max_of_trial = max(single_trial(channel, start_val:end_val)');
      if((max_of_trial-min_of_trial) < wEnvelope)     
        r = 2;
        rvalue = channel;
        break;
        
    else    % otherwise increment the window by the progression
        start_val = start_val + wProgression;
        end_val = end_val + wProgression;
    end
    
end





%% Notes:

% window size 20, moving by 5 -> works best for files with flatlines on
% part of the epochs

% window size 50, move by 5 -> works for files with flatlines for entire
% epoch

