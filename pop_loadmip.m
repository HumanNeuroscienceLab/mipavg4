% pop_loadmip() - load header/data created by MIP for Windows.
%
% Usage:
%   >> EEG = pop_loadmip; % pop-up window mode
%   >> EEG = pop_loadmip( filename, 'key', 'val', ...);
%
% Inputs:
%   filename       - file name
%
% Optional inputs:
%   'headeronly'   - ['yes'|'no'] Read only the header. Default is no.
%   'conversion'   - Override the microvolt conversion factor in the
%                       header with specified value (optional).
%   'labels'       - An optional cell array with labels for each electrode.
%                    If there are 34 electrodes, it will use a default set 
%                    of labels. For all other setups, the electrodes will 
%                    be numbered.
%   'labels_file'  - A file with a list of labels (aka a montage file).
%
% Outputs:
%   EEG            - EEGLAB data structure
%   command        - commands run, useful to add to ALLCMDS
%
% Author: Zarrar Shehzad, 2012

function [EEG, command] = pop_loadmip(filename, varargin) 

ADunits = 16384; adunitsPerUvolt = 1;

cap34 = {
    {'HEOG', -18, 0.768, 0.95, 0.309, -0.0349, 18, -2, 1,'EOG'},...
    {'VEOG', 18, 0.768, 0.95, -0.309, -0.0349, -18, -2, 1,'EOG'},...
    {'FPZ', 0, 0.511, 0.999, -0, -0.0349, -0, -2, 1, 'EEG'},...
    {'FZ', 0, 0.256, 0.719, -0, 0.695, -0, 44, 1, 'EEG'},...
    {'FCz', 0, 0.128, 0.391, -0, 0.921, -0, 67, 1, 'EEG'},...
    {'Cz', 90, 0, 3.75e-33,	-6.12e-17, 1, -90, 90, 1, 'EEG'},...
    {'CPZ', 180, 0.128, -0.391,	-4.79e-17, 0.921, -180, 67, 1, 'EEG'},...
    {'PZ', 180, 0.256, -0.719, -8.81e-17, 0.695, -180, 44, 1, 'EEG'},...
    {'FP1', -18, 0.511, 0.95, 0.309, -0.0349, 18, -2, 1, 'EEG'},...
    {'FP2', 18, 0.511, 0.95, -0.309, -0.0349, -18, -2, 1, 'EEG'},...
    {'F3', -39, 0.333, 0.673, 0.545, 0.5, 39, 30, 1, 'EEG'},...
    {'F4', 39, 0.333, 0.673, -0.545, 0.5, -39, 30, 1 ,'EEG'},...
    {'F7', -54, 0.511, 0.587, 0.809, -0.0349, 54, -2, 1,'EEG'},...
    {'F8', 54, 0.511, 0.587, -0.809, -0.0349, -54, -2, 1,'EEG'},...
    {'FC3', -62, 0.278, 0.36, 0.676, 0.643, 62, 40, 1,'EEG'},...
    {'FC4', 62, 0.278, 0.36, -0.676, 0.643, -62, 40, 1,'EEG'},...
    {'FT7', -72, 0.511, 0.309, 0.95, -0.0349, 72, -2, 1,'EEG'},...
    {'FT8', 72, 0.511, 0.309, -0.95, -0.0349, -72, -2, 1,'EEG'},...
    {'C3', -90, 0.256, 4.4e-17, 0.719, 0.695, 90, 44, 1,'EEG'},...
    {'C4', 90, 0.256, 4.4e-17, -0.719, 0.695, -90, 44, 1,'EEG'},...
    {'T7', -90, 0.511, 6.12e-17, 0.999, -0.0349, 90, -2, 1,'EEG'},...
    {'T8', 90, 0.511, 6.12e-17, -0.999, -0.0349, -90, -2, 1,'EEG'},...
    {'CP3', -118, 0.278, -0.36, 0.676, 0.643, 118, 40, 1, 'EEG'},...
    {'CP4', 118, 0.278, -0.36, -0.676, 0.643, -118, 40, 1, 'EEG'},...
    {'TP7', -108, 0.511, -0.309, 0.95, -0.0349,	108, -2, 1, 'EEG'},...
    {'TP8', 108, 0.511, -0.309, -0.95, -0.0349,	-108, -2, 1, 'EEG'},...
    {'P3', -141, 0.333, -0.673, 0.545, 0.5,	141, 30, 1, 'EEG'},...
    {'P4', 141, 0.333, -0.673, -0.545, 0.5,	-141, 30, 1, 'EEG'},...
    {'P7', -126, 0.511, -0.587, 0.809, -0.0349,	126, -2, 1, 'EEG'},...
    {'P8', 126, 0.511, -0.587, -0.809, -0.0349,	-126, -2, 1, 'EEG'},...
    {'O1', -162, 0.511, -0.95, 0.309, -0.0349,	162, -2, 1, 'EEG'},...
    {'Oz', 180, 0.511, -0.999, -1.22e-16, -0.0349, -180, -2, 1, 'EEG'},...
    {'O2', 162, 0.511, -0.95, -0.309, -0.0349, -162, -2, 1, 'EEG'},...
    {'Blank', 0, 0, 0, 0, 0, 0, 0, 0,'Blank'},...
};


%% Check optional arguments and set any defaults
g = finputcheck(varargin, ...
                {'headeronly'   'string'    {'yes', 'no'}   'no'; ...
                 'conversion'   'integer'   []              []  ; ...
                 'labels'       'cell'      []              {}  ;
                 'labels_file'  'string'    []              ''  });
if ischar(g); error(g); end


%% Load data

if ~exist(filename, 'file'); error('EEG file "%s" doesn''t exist'); end

EEG = eeg_emptyset;
[~, name, ~] = fileparts(filename);
EEG.setname = name;
EEG.comments = sprintf('loaded (%s) with mipavg', filename);

if strcmp(g.headeronly, 'yes')
    eeg = MipRead2(filename, 0);
else
    eeg = MipRead2(filename, 1);
end

if ~strfind(eeg.DAPinfo.DAPname,'4200a/526')
    fprintf('\nUnexpected DAP processor - %s - check A/D units and calibrations\n', ...
            eeg.DAPinfo.DAPName);
end

EEG.subject = deblank(eeg.runInfo.subject);
EEG.nbchan = eeg.nChannels;
EEG.trials = 1;
EEG.pnts = eeg.nPoints;
EEG.srate = eeg.design.freq;
EEG.xmin = 0;
EEG.xmax = EEG.pnts ./ EEG.srate;

if strcmp(g.headeronly, 'yes')
    % save the trigger data
    trigChan = squeeze(eeg.data);    
else
    % get conversion data
    uvcf(1:eeg.nChannels) = (eeg.calInfo.pp(1:eeg.nChannels) ./ eeg.calInfo.ppv) .* ADunits;
    uvcf(1:eeg.nChannels) = uvcf(1:eeg.nChannels) ./ eeg.calInfo.uVolts;
    
    % convert to uVolts
    if ~isempty(g.conversion); uvcf(:) = g.conversion; end
    for k = 1:eeg.nChannels
        % the division by 4 downshifts from 16-bits to 14-bits as per A/D
        divisor = 4 * uvcf(k) * (1/adunitsPerUvolt);
        eeg.data(k,:) = eeg.data(k,:) ./ divisor;
    end
    
    % save the EEG data
    EEG.data = eeg.data(1:eeg.nChannels,:);
    
    % save the trigger data
    trigChan = squeeze(eeg.data(end,:));
end


%% Deal with channel info

% Set default channel names
if ~isempty(g.labels_file)
    if ~exist(g.labels_file, 'file')
        error('Montage/labels file (%s) doesn''t exist', g.labels_file);
    end
    fidMON = fopen(g.labels_file);
    [g.labels, ~] = getstrs(fidMON);
    fclose(fidMON);
elseif isempty(g.labels)
    if eeg.nChannels == 34
        for k=1:eeg.nChannels; g.labels{k} = cap34{k}{1}; end 
    else
        for k=1:eeg.nChannels; g.labels{k} = sprintf('chan%03d', k); end 
    end
elseif length(g.labels) ~= eeg.nChannels
    error('Number of labels does not match number of channels');
end

EEG.chanlocs = struct('theta', cell(1,EEG.nbchan), 'radius', cell(1,EEG.nbchan), ...
                  'labels', cell(1,EEG.nbchan), 'sph_theta', cell(1,EEG.nbchan), ...
                  'sph_phi', cell(1,EEG.nbchan), 'sph_radius', cell(1,EEG.nbchan), ...
                  'X', cell(1,EEG.nbchan), 'Y', cell(1,EEG.nbchan), 'Z', cell(1,EEG.nbchan), ...
                  'type', cell(1,EEG.nbchan));

if eeg.nChannels == 34
    for k = 1:eeg.nChannels
        EEG.chanlocs(k).theta = cap34{k}{2};
        EEG.chanlocs(k).radius = cap34{k}{3};
        EEG.chanlocs(k).labels = char(g.labels{k});
        EEG.chanlocs(k).sph_theta = cap34{k}{7};
        EEG.chanlocs(k).sph_phi = cap34{k}{8};
        EEG.chanlocs(k).sph_radius = cap34{k}{9};
        EEG.chanlocs(k).X = cap34{k}{4};
        EEG.chanlocs(k).Y = cap34{k}{5};
        EEG.chanlocs(k).Z = cap34{k}{6};
        EEG.chanlocs(k).type = char(cap34{k}{10});
    end
else
    nrow = ceil(sqrt(eeg.nChannels));
    nrow_center = (nrow + 1) / 2;
    ncol = ceil(eeg.nChannels/nrow);
    ncol_center = (ncol + 1) / 2;
    
    for k = 1:eeg.nChannels
        EEG.chanlocs(k).theta = 0;
        EEG.chanlocs(k).radius = 0;
        EEG.chanlocs(k).labels = char(g.labels{k});
        EEG.chanlocs(k).sph_theta = 0;
        EEG.chanlocs(k).sph_phi = 0;
        EEG.chanlocs(k).sph_radius = 0;
        nr = rem(k,nrow);
        if (nr == 0); nr = nrow; end
        EEG.chanlocs(k).X = nr - nrow_center;
        nc = mod(ceil(k/nrow),ncol);
        if (nc == 0); nc = ncol; end
        EEG.chanlocs(k).Y = nc - ncol_center;
        EEG.chanlocs(k).Z = 1;
    end
end


%% Process event codes

points = find(trigChan > 0 & trigChan < 251);
codes = trigChan(points);
nEvents = length(codes);

EEG.event = struct('type', cell(1,nEvents), 'latency', cell(1,nEvents));

for e=1:nEvents
   EEG.event(e).type = codes(e);
   EEG.event(e).latency = points(e);
end


%% Check data structure and save command that was called

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

if ~isempty(varargin)
    command = sprintf('EEG = pop_loadcnt(''%s'' %s);', filename, vararg2str(varargin)); 
else
    command = sprintf('EEG = pop_loadcnt(''%s'');', filename); 
end;


return;

% for testing can use:
% filename = '/mnt/nfs/psych/sponty01/raw/eeg/subject101_sponty_task_FIX.EEG';

end