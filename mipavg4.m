function mipavg4(varargin)
%     MIPAVG4   Creates signal averaged ERPs from continuous EEG data
%         mipavg4(varargin)
% 
%     INPUTS (in any order):
%       EEG data file   - Can be a *.eeg, *.cnt, or *.set file (required).
%       SGC file        - A *.sgc (segmentation control) file specifies the 
%                         parsing of the raw EEG file into single trial epochs 
%                         (required).
%       ARF file        - A *.arf (artifact rejection) file specifies the 
%                         criteria and settings for rejecting epochs due to 
%                         detected artifacts (optional).
%       LOG file        - A *.log (log) file to which all the command output 
%                         is written (default: OUTprefix_mipavg.log).
%       RCD file        - A *.rcd (recode) file that substitutes new codes for
%                         the original event codes provided in the EEG data
%                         file (optional).
%       -b path1;path2  - Specify the paths to EEGLAB (path1) and Fieldtrip 
%                         (path2), overriding the defaults.
%       -p OUTprefix    - Specify the output prefix (default: MipAvg).
%       -o OUTtype      - Specify the type of file to output, which can be 
%                         either: raw_codes, eegad, eeglab, fieldtrip, econnectome, 
%                         continuous_eeglab, or continuous_fieldtrip. This 
%                         option can be called multiple times, once for each
%                         output type (required at least once).
%       -r reref;ignore - Channels to use as the reference, which can be 
%                         optionally followed by a semicolon with channels to 
%                         ignore when re-referencing. The channels can be 
%                         anything accepted by ft_channelselection. Multiple 
%                         channels can be specified and should be seperated by 
%                         commas (optional).
%       -f low,high     - Apply band-pass filtering prior to signal averaging. 
%                         The low-frequency cutoff and the high-frequncy cutoff
%                         must be given. Either cutoff can be set to zero, in 
%                         which case that filtering is not done (optional).
%       -n              - Apply a 60, 120, 180 Hz notch filter on the data prior
%                         to signal averaging (optional).
%       -c conv-factor  - Override the microvolt conversion factor contained in
%                         the EEG (MIP) file's header with the specified value 
%                         (optional).
%       MON file        - A *.mon (montage) file provides a label for the
%                         channels/electrodes.
%       -q              - Quiet the output written to the command window. The 
%                         usual stuff is still written to the log file.
%       -s              - Suppress logging output to a file
%     
%     Created by Zarrar Shehzad on 2012-09-07.


printhead = @(logg) logg.write('\n\n=================\n');
printfoot = @(logg) logg.write('\n=================\n\n');


%%==============================================================================
%%                                              Initialize file and option flags
%%==============================================================================

nEEGfiles = 0; nSGCfiles = 0; nMONfiles = 0; nLOGfiles = 0; 
nOUTprefixes = 0; nARFfiles = 0; nRCDfiles = 0;

EEGfiles = {}; SGCfiles = {}; RCDfiles = {};
OUTprefix = []; ARFfile = [];

typeEEG = 0; typeCNT = 0; typeSET = 0; 
output = []; output_opts = {'raw_codes', 'eegad', 'eeglab', 'fieldtrip', ...
                            'econnectome', 'continuous_eeglab', ....
                            'continuous_fieldtrip'};
for i=1:length(output_opts); output.(output_opts{i}) = 0; end

lowpass_value = 0; highpass_value = 0; notchfilter_option = 0;
refchannel = []; conversion_factor = 0; verbose = 1; to_logfile = 1;
output_raw_codes = 0;

adunitsPerUvolt = 1; 

basepath = '/usr/local/packages/MATLABPackages';
eeglabpath = [basepath filesep 'eeglab'];
fieldtrippath = [basepath filesep 'fieldtrip'];
eeglabset = 0; fieldtripset = 0;

% old stuff...not sure to keep or not
display_option = 0;

% Arguments for pop_loadXXX functions
load_args = {};


%%==============================================================================
%%                                                          Parse user arguments
%%==============================================================================

for j = 1:length(varargin)
    arg_string = lower(char(varargin{j}));
    [~, ~, arg_ext] = fileparts(arg_string);
    if strcmp(arg_ext, '.eeg')
        EEGfiles{end+1} = varargin{j};
        typeEEG = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif strcmp(arg_ext, '.cnt')
        EEGfiles{end+1} = varargin{j};
        typeCNT = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif strcmp(arg_ext, '.set')
        EEGfiles{end+1} = varargin{j};
        typeSET = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif strcmp(arg_ext, '.sgc')
        SGCfiles{end+1} = varargin{j};
        nSGCfiles = nSGCfiles + 1;
    elseif strcmp(arg_ext, '.log')
        LOGfile = varargin{j};
        nLOGfiles = nLOGfiles + 1;
    elseif strcmp(arg_ext, '.arf')
        ARFfile = varargin{j};
        nARFfiles = nARFfiles + 1;
    elseif strcmp(arg_ext, '.rcd')
        RCDfiles{end+1} = varargin{j};
        nRCDfiles = nRCDfiles + 1;
    elseif strfind(arg_string, '-b')
        opt = strtrim(varargin{j}(strfind(arg_string, '-b')+2:end));
        opt = regexp(opt, ';\ *', 'split');
        if length(opt) ~= 2, error('Unrecognized argument for -b'); end
        [eeglabpath fieldtrippath] = opt{:};        
    elseif strfind(arg_string, '-p')
        OUTprefix = strtrim(varargin{j}(strfind(arg_string, '-p')+2:end));
        nOUTprefixes = nOUTprefixes + 1;
    elseif strfind(arg_string, '-o')
        o = lower(strtrim(varargin{j}(strfind(arg_string,'-o')+2:end)));
        if isempty(intersect(o, output_opts))
            error('Unrecognized argument ''%s'' for -o', o); 
        end
        output.(o) = 1;
    % Preprocessing
    elseif strfind(arg_string, '-r')
        opt = strtrim(arg_string(strfind(arg_string,'-r')+2:end));
        opt = regexp(opt, ';\ *', 'split');
        refchannel = regexp(opt{1}, ',\ *', 'split');
        if isempty(refchannel{1}), refchannel = 'all'; end
        if length(opt) > 1
            refignore = regexp(opt{2}, ',\ *', 'split');
        else
            refignore = [];
        end
    elseif strfind(arg_string, '-f')
        filter_values = arg_string(strfind(arg_string,'-f')+2:end);
        filter_values = str2double(filter_values);
        if length(filter_values) ~= 2; error('Unrecognized argument for -f'); end
        lowpass_value = filter_values(1); highpass_value = filter_values(2);
    elseif strfind(arg_string, '-n')
        notchfilter_option = 1;
%     elseif(strfind(arg_string,'-d'))
%         display_option = 1;
%         temp = char(varargin{j});
%         display_channel = fix(str2double(temp(strfind(temp,'-d')+2:end)));
%         if(isempty(display_channel)) display_channel = 1; end
%         nDisplayChannels = length(display_channel);
    % Legacy options for MIP files
    elseif strfind(arg_string, '-c')
        temp = arg_string;
        conversion_factor = str2double(temp(strfind(temp,'-c')+2:end));
        if isempty(conversion_factor), conversion_factor = 1; end
        load_args{end+1} = 'conversion';
        load_args{end+1} = conversion_factor;
        clear temp conversion_factor;
    elseif strcmp(arg_ext, '.mon')
        load_args{end+1} = 'labels_file';
        load_args{end+1} = varargin{j};
        nMONfiles = nMONfiles + 1;
    % Logging
    elseif strfind(arg_string, '-q')
        verbose = 0;
    elseif strfind(arg_string, '-s')
        to_logfile = 0;
    else
        error('Invalid argument "%s"', arg_string);
    end
end


%%==============================================================================
%%      Check that we have all the necessary files and set defaults where needed
%%==============================================================================

if sum([typeEEG typeCNT typeSET]) > 1
    error('Cannot use any combination of .cnt, .eeg, or .set files');
end

if ~nEEGfiles, error('No EEG file was specified. Program cannot proceed'); end

if ~nSGCfiles
    error('No SGC file was specified. Program cannot proceed');
elseif (nSGCfiles > 1 && nSGCfiles ~= nEEGfiles)
    error('Insufficient SGC files specified to match EEG files. Program cannot proceed');
end

if nMONfiles > 1
  error('More than one montage file was specified. Program cannot proceed');
end

if ~nOUTprefixes
    fprintf('\nNo output prefix was specified. Using MipAvg.\n');
    OUTprefix = 'MipAvg';
elseif nOUTprefixes > 1
    error('More than one output prefix was specified. Program cannot proceed');
end

if ~any(struct2array(output))
    error('No output was specified. Program cannot proceed');
end

if ~nLOGfiles && to_logfile
    fprintf('\nNo LOG file was specified. Using %s.log.\n', OUTprefix);
    LOGfile = [OUTprefix '.log'];
elseif nLOGfiles > 1
    error('More than one LOG file specified. Program cannot proceed');
end

if ~nARFfiles
    fprintf('\nNo ARF file was specified. No artifact rejection will be applied.\n');
elseif nARFfiles > 1
    error('More than one ARF file specified. Program cannot proceed');
end

if ~exist(eeglabpath, 'dir')
    error('Cannot find required eeglab directory %s', eeglabpath);
end
if ~exist(fieldtrippath, 'dir')
    error('Cannot find required fieldtrip directory %s', fieldtrippath);
end

is_writable = can_write(dirname(OUTprefix));
if ~is_writable
    error('Cannot write to the output directory %s', dirname(OUTprefix));
end
clear is_writable;


%%==============================================================================
%%                                                            Initialize logging
%%==============================================================================

logg = logger;

% Log to standard output if verbose
if verbose, logg = logg.to_standard_output; end

% Log to file
if to_logfile, logg = logg.to_file(LOGfile); end

% First messages!
logg.write('\nMIPAVG4 Program executed on %s\n', datestr(now));
opts = sprintf('''%s'', ', varargin{:});
logg.write('\nCalled:\nmipavg4(%s)\n\n', opts(1:end-2));

% Report on options
if (conversion_factor)
    logg.write('\nCalibration override specified. Microvolt conversion factor set equal to %5.2f\n', conversion_factor);
end


%%==============================================================================
%%                                                         Set output file paths
%%==============================================================================

EEGLABcontinuousfile = [OUTprefix '_continuous.set'];
FIELDTRIPcontinuousfile = [OUTprefix '_ft_continuous.mat'];
EEGLABprefix = [OUTprefix '_trials'];
FIELDTRIPprefix = [OUTprefix '_ft_trials'];
ECONNECTOMEprefix = [OUTprefix '_econnectome'];
AVGfile = [OUTprefix '_eegad.avg'];
HDRfile = [OUTprefix '_eegad.hdr'];


%%==============================================================================
%%                                               Add necessary functions to path
%%==============================================================================

if exist('pop_importdata.m','file') ~= 2
    if exist([eeglabpath '/functions/popfunc'],'dir') == 7
        logg.write('\nLoading required EEGLAB functions\n');
        addpath([eeglabpath '/functions/guifunc']);
        addpath([eeglabpath '/functions/popfunc']);
        addpath([eeglabpath '/functions/adminfunc']);
        addpath([eeglabpath '/functions/sigprocfunc']);
        eeglabset = 1;
    else
        error(['Required eeglab functions cannot be found - ' ...
               'check your paths (%s)!'], eeglabpath);
    end
end

if exist('ft_definetrial', 'file') ~= 2
    if exist(fieldtrippath, 'dir') == 7
        logg.write('\nLoading required fieldtrip functions\n');
        addpath(fieldtrippath);
        ft_defaults;
        fieldtripset = 1;
    else
        error(['Required fieldtrip functions cannot be found - ' ...
               'check your paths (%s)!'],  fieldtrippath);
    end
end


%%==============================================================================
%%                                                             Check input files
%%==============================================================================

for fname=[EEGfiles, SGCfiles, ARFfile, RCDfiles]
    check_input_file(fname{1}, logg);
end


%%==============================================================================
%%                                                              Read in EEG data
%%==============================================================================

printhead(logg);
logg.write('Reading in EEG files\n');

ALLEEG = [];
for i=1:nEEGfiles
    logg.write('\nFile #%d: %s\n', i, EEGfiles{i});
    if typeEEG
        EEG = logg.callfun('pop_loadmip', EEGfiles{i}, load_args{:});
    elseif typeCNT
        EEG = logg.callfun('pop_loadcnt', EEGfiles{i}, 'dataformat', 'int32');
    elseif typeSET
        EEG = logg.callfun('pop_loadset', EEGfiles{i});
    end
    [ALLEEG EEG ~] = eeg_store(ALLEEG, EEG);
end

printfoot(logg);



%%==============================================================================
%%                                                      Recode events (optional)
%%==============================================================================

% Output raw event codes that can be used for creating RCD file
if output.raw_codes
    logg.write('Saving raw event codes\n');
    for i=1:nEEGfiles
        fname               = [OUTprefix '_rawcodes_file' num2str(i) '.txt'];
        rawevents           = [[EEG.event.type]', [EEG.event.latency]'];
        dlmwrite(fname, rawevents, 'delimiter', '\t');
        logg.write('...%s', fname);
        clear fname raw_event;
    end
end

if nRCDfiles
    printhead(logg);
    logg.write('Recoding section\n');
    
    ALLrecodes  = read_recode(RCDfiles, logg);
    ALLEEG      = recode_eeg(ALLEEG, ALLrecodes, logg);
    
    printfoot(logg);
end


%%==============================================================================
%%            Read in SGC file with epoch/baseline limits and trial codes/labels
%%==============================================================================

printhead(logg);
logg.write('Reading in segmentation control file(s)\n')

[epoch baseline first_pt_msec ALLcodes ALLlabels] = read_sgc(SGCfiles, logg);

printfoot(logg);


%%==============================================================================
%%                                Map codes and labels from SGC file to EEG file
%%==============================================================================

printhead(logg);
logg.write('Mapping codes and labels from SGC to EEG');

ALLEEG = map_codes_and_labels(ALLEEG, ALLcodes, ALLlabels, epoch, logg);

printfoot(logg);


%%==============================================================================
%%                      Merge data, setup channel info, and save continuous data
%%==============================================================================

printhead(logg);

% Merge data together
logg.write('\nMerging data together\n');
if nEEGfiles > 1, 
    EEG = logg.callfun('pop_mergeset', ALLEEG, 1:nEEGfiles);
    % also merge unique event codes/labels
    EEG.eventcodes      = [ALLEEG.eventcodes];
    EEG.eventlabels     = [ALLEEG.eventlabels];
    % only keep unique codes/labels from merger
    [EEG.eventcodes I]  = unique(EEG.eventcodes);
    EEG.eventlabels     = EEG.eventlabels(I);
    clear I;
else
    EEG = ALLEEG(1); 
end
% clear ALLEEG;

% Ghetto fix of channel labels for cnt files
logg.write('\nFixing channel labels\n');
old_vals = {'CB1', 'CB2', 'HEO', 'VEO'};
new_vals = {'I1', 'I2', 'HEOG', 'VEOG'};
for i=1:length(old_vals)
    search = ismember({EEG.chanlocs.labels}, old_vals{i});
    if any(search); 
        EEG.chanlocs(search).labels = new_vals{i};
        logg.write('\tchanging %s => %s\n', old_vals{i}, new_vals{i});
    end
    clear search;
end
clear old_vals new_vals;

% Get channel locations
logg.write('\nLooking up standard channel locations\n');
chanfile = [eeglabpath '/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'];
EEG = logg.callfun('pop_chanedit', EEG, 'lookup', chanfile);

% Save the EEGLab continuous file (will be removed later if user-specified)
logg.write('\nSaving continuous EEGLab file\n');
check_output_file(EEGLABcontinuousfile, logg);
logg.callfun('pop_saveset', EEG, 'filename', EEGLABcontinuousfile);

% Save the fieldtrip continuous file
if output.continuous_fieldtrip
    logg.write('\nSaving continuous fieldtrip file\n');
    check_output_file(FIELDTRIPcontinuousfile, logg);
    data            = logg.callfun('eeglab2fieldtrip', EEG, 'preprocessing');
    data.hdr        = ft_read_header(EEGLABcontinuousfile);
    data.cfg.event  = ft_read_event(EEGLABcontinuousfile);
    save(FIELDTRIPcontinuousfile, 'data');
end

printfoot(logg);


%%==============================================================================
%%                                                    Epoching and Preprocessing
%%==============================================================================

printhead(logg);

% Read continuous data into memory and apply filtering
logg.write('\nReading in continuous data\n');
cfg = [];
cfg.dataset                 = EEGLABcontinuousfile;
cfg.continuous              = 'true';
if lowpass_value && highpass_value
    cfg.bpfilter            = 'yes';
    cfg.bpfreq              = [lowpass_value highpass_value];
elseif lowpass_value
    cfg.hpfilter            = 'yes';
    cfg.hpfreq              = lowpass_value;
elseif highpass_value
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = highpass_value;
end
if notchfilter_option
    cfg.dftfilter           = 'yes';
    cfg.dftfreq             = [60, 120, 180];
end
if lowpass_value || highpass_value || notchfilter_option
    logg.write('...will filter the data as well\n');
end
continuous_data = logg.callfun('ft_preprocessing', cfg);

% Re-Reference
if ~isempty(refchannel)
    logg.write('\nRe-referencing data\n');
    cfg = [];
    cfg.reref               = 'yes';
    cfg.refchannel          = refchannel;
    if ~isempty(refignore)
        % re-reference desired channels
        cfg.channel         =  cellfun(@(x) ['-' x], refignore, 'UniformOutput', false);
        cfg.channel         = {'all', cfg.channel{:}};
        data_chansA         = logg.callfun('ft_preprocessing', cfg, continuous_data);
        % ignore desired channels
        cfg.reref           = 'no';
        cfg.channel         = refignore;
        data_chansB         = logg.callfun('ft_preprocessing', cfg, continuous_data);
        % combine back together
        continuous_data     = logg.callfun('ft_appenddata', cfg, data_chansA, data_chansB);
        clear data_chansA data_chansB;
    else
        continuous_data     = logg.callfun('ft_preprocessing', cfg, continuous_data);
    end
end

% Define trials
logg.write('\nDefining trials\n');
cfg = [];
cfg.dataset                 = EEGLABcontinuousfile;
cfg.continuous              = 'true';
cfg.trialfun                = 'mip_trialfun';
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = EEG.eventcodes;
cfg.trialdef.prestim        = -1 * epoch(1);
cfg.trialdef.poststim       = epoch(2);
cfg = logg.callfun('ft_definetrial', cfg);
tmpevent = cfg.event;
tmptrl   = cfg.trl;

% Split data into trials
logg.write('\nSplitting continuous data into trials\n');
data = ft_redefinetrial(cfg, continuous_data);
clear continuous_data;

% Check first time point for any rounding
if data.time{1}(1) ~= epoch(1)
    logg.warn('MATLAB:ft_redefinetrial:Rounding', ...
              'First time point was adjusted from %.4f to %.4f.', ...
              epoch(1), data.time{1}(1));
end

% Remove baseline
logg.write('\nRemoving baseline\n');
cfg = [];
cfg.demean                  = 'yes';
cfg.baseline                = baseline;
data = logg.callfun('ft_preprocessing', cfg, data);

% Add back event and trl structure
data.cfg.event  = tmpevent;
data.cfg.trl    = tmptrl;
clear tmpevent tmptrl;

% Remove continuous EEGLAB file
if ~output.continuous_eeglab
    logg.write('\nRemoving continuous set file: %s\n', EEGLABcontinuousfile);
    [pathstr name] = fileparts(EEGLABcontinuousfile);
    delete(EEGLABcontinuousfile);
    if isempty(pathstr), pathstr = '.'; end
    delete([pathstr filesep name '.fdt']);
end

printfoot(logg);


%%==============================================================================
%%                                                            Artifact Rejection
%%==============================================================================

if nARFfiles
    printhead(logg);
    logg.write('\nArtifact Rejection\n');
    
    % Read in artifact rejection settings
    ar_settings = read_arf(ARFfile, data.label, epoch, logg);

    % Loop through and apply each user-defined artifact rejection criteria
    cfg = [];
    cfg.continuous = 'no';
    for i=1:length(ar_settings)
        switch ar_settings(i).method
            case 'flat'
                cfg.trl                             = data.cfg.trl;
                cfg.artfctdef.clip.channel          = ar_settings(i).channels;
                cfg.artfctdef.clip.prestim          = ar_settings(i).prestim;
                cfg.artfctdef.clip.poststim         = ar_settings(i).poststim;
                cfg.artfctdef.clip.thresh           = ar_settings(i).criteria;
                cfg = logg.callfun('run_ft_artifact', cfg, data, 'clip');
                logg.write('*** %d trials rejected with flat ***\n', size(cfg.artfctdef.clip.artifact, 1));
            case 'ppa'
                cfg.trl                             = data.cfg.trl;
                % adjust starting sample based on prestim
                if ar_settings(i).prestim < -1*epoch(1)
                    cfg.trl(:,1) = cfg.trl(:,1) + (-1*epoch(1) - ar_settings(i).prestim) * data.fsample;
                end
                % adjust ending sample based on poststim
                if ar_settings(i).poststim < epoch(2)
                    cfg.trl(:,2) = cfg.trl(:,2) + (ar_settings(i).poststim - epoch(2)) * data.fsample;
                end
                cfg.artfctdef.threshold.channel     = ar_settings(i).channels;
                cfg.artfctdef.threshold.range       = ar_settings(i).criteria;
                cfg.artfctdef.threshold.bpfilter    = 'no';
                for j=1:length(ar_settings(i).opts)
                    if strcmpi(ar_settings(i).opts{j}, 'yes-freq')
                        cfg.artfctdef.threshold.bpfilter = 'yes';
                    end
                end
                % TODO: are the other filter settings good?
                cfg = logg.callfun('run_ft_artifact', cfg, data, 'threshold');
                logg.write('*** %d trials rejected with ppa ***\n', ...
                            size(cfg.artfctdef.threshold.artifact, 1));
            case 'zthr'
                cfg.trl                             = data.cfg.trl;
                % adjust starting sample based on prestim
                if ar_settings(i).prestim < -1*epoch(1)
                    cfg.trl(:,1) = cfg.trl(:,1) + (-1*epoch(1) - ar_settings(i).prestim) * data.fsample;
                end
                % adjust ending sample based on poststim
                if ar_settings(i).poststim < epoch(2)
                    cfg.trl(:,2) = cfg.trl(:,2) + (ar_settings(i).poststim - epoch(2)) * data.fsample;
                end
                cfg.artfctdef.eog.channel           = ar_settings(i).channels;
                cfg.artfctdef.eog.cutoff            = ar_settings(i).criteria; % z-value
                cfg.artfctdef.eog.bpfilter          = 'yes';
                cfg.artfctdef.eog.bpfilttype        = 'fir';
                cfg.artfctdef.eog.hilbert           = 'yes';
                cfg.artfctdef.eog.interactive       = 'no';
                if  (-1*epoch(1) - ar_settings(i).prestim) < 0.1 || (ar_settings(i).poststim - epoch(2)) < 0.1
                    cfg.artfctdef.eog.trlpadding    = -0.1; % Padding added for each trial in 
                                                            %   checking for an artifact.                
                end
                cfg.artfctdef.eog.fltpadding        = 0.1;  % Padding added before filtering
                                                            %   and removed after filtering.
                cfg.artfctdef.eog.artpadding        = 0.1;  % Padding added around period of time
                                                            %   determined to be an artifact,
                for j=1:length(ar_settings(i).opts)
                    if strcmpi(ar_settings(i).opts{j}, 'no-freq')
                        cfg.artfctdef.eog.bpfilter  = 'no';
                        cfg.artfctdef.eog.bpfilttype= 'fir';
                        cfg.artfctdef.eog.hilbert   = 'no';
                    elseif strcmpi(ar_settings(i).opts{j}, 'yes-interactive')
                        cfg.artfctdef.eog.interactive = 'yes';
                    end
                end
                % TODO: does there need to be an option for the filt freq range?
                cfg = logg.callfun('run_ft_artifact', cfg, data, 'eog');
                logg.write('*** %d trials rejected with zthr ***\n', ...
                            size(cfg.artfctdef.eog.artifact, 1));
            otherwise
                logg.error('Unrecognized artifact rejection approach: %s\n', ...
                            ar_settings(i).method);
                error('Program cannot proceed');
        end
    end

    % Remove offending trials
    cfg.artfctdef.reject    = 'complete';
    tmpevent = data.cfg.event;
    data = logg.callfun('ft_rejectartifact', cfg, data);    
    data.cfg.event = tmpevent;
    clear tmpevent;
    
    printfoot(logg);
else
    data.cfg.trlold = data.cfg.trl;
end

% Indices of good trials
% cleaninds = arrayfun(@(x) find(x == data.cfg.trlold(:,1)), data.cfg.trl(:,1));

% Adjust time-axis
printhead(logg);
data = set_first_time(data, first_pt_msec/1000, logg);
printfoot(logg);


%%==============================================================================
%%                                       Write Output (EEGAD, EEGLab, Fieldtrip)
%%==============================================================================

printhead(logg);
logg.write('\nFor each unique event code, save single trial and averaged output\n');

ntpts   = round(diff(epoch) * data.hdr.Fs) + 1;
nchans  = data.hdr.nChans;
ncodes  = length(EEG.eventcodes);

if output.eegad
    % Create average ERP matrix for EEGAD
    avg     = zeros(nchans+1, ntpts, ncodes);

    % AVG file (EEGAD Output)
    logg.write('\nCreating AVG file: %s\n', AVGfile);
    check_output_file(AVGfile, logg);
    fidAVG = fopen(AVGfile,'w');
end

% Split trials based on unique event code
for i=1:ncodes
    code = EEG.eventcodes(i); label = EEG.eventlabels{i};
    logg.write('\nProcessing code %d or label %s\n', code, label);
    
    %% Fieldtrip
    % Select trials with unique event code
    trials = data.cfg.trl(:,4) == code;
    tdata = ft_selectdata(data, 'rpt', trials);
    if output.fieldtrip
        logg.write('\n\tSaving fieldtrip output\n')
        save([FIELDTRIPprefix '_' label '.mat'], 'tdata');
    end
    
    
    %% EEGLAB
    if output.eeglab
        % Convert from fieldtrip to eeglab and then save
        logg.write('\n\tSaving EEGLAB output\n');
        EEGtrial = fieldtrip2eeglab(tdata, EEG, code, label, logg);
        logg.callfun('pop_saveset', EEGtrial, 'filename', [EEGLABprefix '_' label '.set']);
    end
    
    
    %% Average ERP
    
    % Average single trials
    logg.write('\n\tAveraging trials\n');
    cfg = [];
    timelock = logg.callfun('ft_timelockanalysis', cfg, tdata);
    
    % Baseline correct
    logg.write('\n\tAdjusting based on the baseline\n');
    cfg = [];
    cfg.baseline = baseline;
    timelock = logg.callfun('ft_timelockbaseline', cfg, timelock);
    
    
    %% Econnectome
    if output.econnectome
        logg.write('\n\tSaving econnectome output\n');
        tmp = EEG;
        EEG = fieldtrip2econnectome(timelock, label, data.hdr.Fs);
        save([ECONNECTOMEprefix '_' label '.mat'], 'EEG');
        EEG = tmp; clear tmp;
    end
    
    
    %% EEGAD
    if output.eegad
        logg.write('\n\tSaving EEGAD output\n');
        
        % Retain averaged evoked reponse
        avg(1:nchans,:,i) = timelock.avg;
    
        % Add trigger channel info
        % (basically the event code at the onset)
        onset = -epoch(1) * data.hdr.Fs + 1;
        avg(nchans+1,onset,i) = code;
        
        % Save for EEGAD
        fwrite(fidAVG, avg(:,:,i)', 'float32');
    end
end

%% Log Rejected Trials
logg.write('\n Code\t Label\t Accepted #\t Total #\t Percent Averaged\n');
for j = 1:length(EEG.eventcodes)
    code = EEG.eventcodes(j); label = EEG.eventlabels{j};
    total = sum(data.cfg.trlold(:,4)==code);
    accepted = sum(data.cfg.trl(:,4)==code);
    logg.write(' %d\t %s\t %d\t\t %d\t\t %4.2f\n', code, label, ...
                accepted, total, accepted/total*100);
end

%% HDR File (EEGAD Output)
if output.eegad
    logg.write('\nCreating AVG HDR file: %s\n', HDRfile);
    check_output_file(HDRfile, logg);
    fidHDR = fopen(HDRfile, 'w');
    
    % Line 1: expName is the experiment name
    fprintf(fidHDR,'%s\n', OUTprefix); 
    % Line 2: expDate is the experiment date string
    fprintf(fidHDR,'%s\n', datestr(now, 'dd-mmm-yyyy HH:MM:SS'));
    % Line 3: nChannels is the number of data channels (electrodes).
    fprintf(fidHDR,'%d\n', nchans+1);
    % Line 4: nPoints is the number of data points collected per channel.
    fprintf(fidHDR,'%d\n', ntpts);
    % Line 5: sampling is the sampling rate of the data points in ms/point.
    fprintf(fidHDR,'%d\n', 1000/data.hdr.Fs);
    % Line 6: uvunits is the microvolt conversion factor in raw data units/microvolt.
    fprintf(fidHDR,'%d\n', adunitsPerUvolt);
    % Line 7: onset is the stimulus onset from the first data point in ms.
    fprintf(fidHDR,'%d\n', -first_pt_msec);
    % Line 8+: Event Code Labels
    for k = 1:length(EEG.eventcodes)
        label = char(EEG.eventlabels{k});
        fprintf(fidHDR,'%s\n',label);
    end
    % Line X+: Channel Names
    for k = 1:nchans, fprintf(fidHDR,'%s\n', data.label{k}); end
    fprintf(fidHDR,'%s\n', 'trigger');
end

printfoot(logg);


%%==============================================================================
%%                                                    Remove functions from path
%%==============================================================================

if eeglabset
    rmpath([eeglabpath '/functions/popfunc']);
    rmpath([eeglabpath '/functions/adminfunc']);
    rmpath([eeglabpath '/functions/guifunc']);
end

if fieldtripset
    warning('off', 'MATLAB:rmpath:DirNotFound');
    rmpath(genpath(fieldtrippath));
    warning('on', 'MATLAB:rmpath:DirNotFound');
end


end
