function ft_event = eeglab2fieldtrip_for_event(eeglab_event)
%   [ft_event] = event_eeglab2fieldtrip_for_event(eeglab_event)
%   
%   Converts the EEGlab event structure to Fieldtrip event structure.
%   This code is mostly copied from the read_eeglabevent.m script.
%   
%   Created by Zarrar Shehzad on 2012-09-07.
%

ft_event = [];

for index = 1:length(eeglab_event)

  if isfield(eeglab_event,'code')
    type = eeglab_event(index).code;
  elseif isfield(eeglab_event,'value')
    type = eeglab_event(index).value;
  else
    type = 'trigger';
  end;

  % events can have a numeric or a string value
  if isfield(eeglab_event,'type')
    value  = eeglab_event(index).type;
  else
    value = 'default';
  end;

  % this is the sample number of the concatenated data to which the event corresponds
  sample = eeglab_event(index).latency;

  % a non-zero offset only applies to trial-events, i.e. in case the data is
  % segmented and each data segment needs to be represented as event. In
  % that case the offset corresponds to the baseline duration (times -1).
  offset = 0;

  if isfield(eeglab_event, 'duration')
    duration = eeglab_event(index).duration;
  else
    duration = 0;
  end;
  
  % add the current event in fieldtrip format
  ft_event(index).type     = type;     % this is usually a string, e.g. 'trigger' or 'trial'
  ft_event(index).value    = value;    % in case of a trigger, this is the value
  ft_event(index).sample   = sample;   % this is the sample in the datafile at which the event happens
  ft_event(index).offset   = offset;   % some events should be represented with a shifted time-axix, e.g. a trial with a baseline period
  ft_event(index).duration = duration; % some events have a duration, such as a trial

end

end %  function


