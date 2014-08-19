function [EEG] = fieldtrip2eeglab(data, contEEG, eventcode, eventlabel, logg)
%   function [EEG] = fieldtrip2eeglab(data, contEEG, eventcode, eventlabel, logg)
% TODO
%

trl = data.cfg.trl(data.cfg.trl(:,4) == eventcode,:);

if size(trl, 1) ~= size(data.trial,2)
    error('Mismatch between size of trl (%d) and size of trial (%d)',  ...
            size(trl,1), size(data.trial,2));
end

EEG = eeg_emptyset;

EEG.setname         = eventlabel;
EEG.comments        = sprintf('preprocessed with fieldtrip\n\nParent dataset: %s\n%s', ...
                                contEEG.setname, contEEG.comments);
EEG.nbchan          = size(data.trial{1},1);
EEG.trials          = size(data.trial,2);
EEG.pnts            = size(data.trial{1},2);
EEG.srate           = data.fsample;
EEG.xmin            = data.time{1}(1);
EEG.xmax            = data.time{1}(end);
EEG.times           = data.time{1};
for i=1:size(data.trial,2)
    EEG.data(:,:,i) = single(data.trial{i});
end
EEG.chanlocs        = contEEG.chanlocs;
EEG.urchanlocs      = contEEG.urchanlocs;
EEG.chaninfo        = contEEG.chaninfo;
EEG.ref             = contEEG.ref;

% Include events from continuous EEG, 
% which are found ANYTIME during each trial epoch
oldevent            = contEEG.event;
[oldevent.epoch]    = deal([]);
newevent            = struct(oldevent(1));
for i=1:size(data.trial,2)
    oldinds = find([oldevent.latency] >= trl(i,1) & [oldevent.latency] <= trl(i,2));
    if isempty(oldinds)
        logg.warn('MATLAB:fieldtrip2eeglab:NoEvent', ...
                  'Did not find any event(s) for trial #%d', i);
    else
        n                           = length(newevent);
        newinds                     = (n+1):(n+length(oldinds));
        newevent(newinds)           = oldevent(oldinds);
        [newevent(newinds).epoch]   = deal(i);
        onset_latency               = trl(i,1) - trl(i,3);
        [newevent(newinds).latency] = split([oldevent(oldinds).latency] - onset_latency ...
                                            + -EEG.xmin*EEG.srate + 1 + EEG.pnts*(i-1));
    end
end
EEG.event       = newevent(2:end);
EEG.epoch       = [];
EEG             = eeg_checkset(EEG, 'eventconsistency');    % populates epoch structure
EEG.urevent     = contEEG.urevent;

EEG.history     = contEEG.history;
EEG.saved       = 'no';


function varargout = split(a)
    for k = 1:nargout
      varargout{k} = a(k);
    end
end

end