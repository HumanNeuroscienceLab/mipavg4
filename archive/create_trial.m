function cfg = create_trial(cfg, Fs, fidLOG)
%   cfg = create_trial(cfg)
% 
%   TODO
%   
%   Created by Zarrar Shehzad on 2012-09-07.
%
% cfg.event
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

if nargin < 3; fidLOG = 1; end

if isnumeric(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = num2strcell(cfg.trialdef.eventvalue); 
end

sel_events = cfg.event(ismember({cfg.event.type}, cfg.trialdef.eventtype));
sel_events = sel_events(ismember({sel_events.value}, cfg.trialdef.eventvalue));

if length(sel_events) == 0
    fprintf('\nCould not find any trials with type %s and value %s\n',  ...
            cfg.trialdef.eventtype, char(cfg.trialdef.eventvalue));
    error('Program cannot continue');
end

prestim_pts = round(cfg.trialdef.prestim * Fs);
tmp_prestim = prestim_pts/Fs;
if (cfg.trialdef.prestim ~= tmp_prestim)
    fprintf('Adjusted prestim from %d to %d secs', cfg.trialdef.prestim, tmp_prestim);
    cfg.trialdef.prestim = tmp_prestim;
end

poststim_pts = round(cfg.trialdef.poststim * Fs);
tmp_poststim = poststim_pts/Fs;
if (cfg.trialdef.poststim ~= tmp_poststim)
    fprintf('Adjusted poststim from %d to %d secs', cfg.trialdef.poststim, tmp_poststim);
    cfg.trialdef.poststim = tmp_poststim;
end

cfg.trl = [[sel_events.sample] - prestim_pts; [sel_events.sample] + poststim_pts; ...
            repmat(-1*prestim_pts, 1, length(sel_events))]';

end %  function

