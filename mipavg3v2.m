function mipavg3v2(varargin)
%%
% Creates bin averages of raw EEG file acquired by MIPWIN
% mipavg2 released on 17-March-2008
% mipavg2 modified to include residual option 23-January-2009
% mipavg3 modified to fix minor problems on 10-April-2009
% mipavg3 modified to fix minor error with 'p' option and to add codes labels to output files on 3-December-2009
% mipavg3 modified to write eeglab .set files instead of .mat files

% varargin = {'SampleEEG.eeg', 'SampleMON.mon', 'SampleSGC.sgc', 'SampleAVG.avg', '-d1,2,13',   'SampleARF.arf',  'SampleLOG.log'};

% Gregory McCarthy
% 2001-2011

tic;

%% Initialize file and option flags

arg_pointer = zeros(length(varargin),1);

nEEGfiles = 0; nSGCfiles = 0; nMONfiles = 0; nLOGfiles = 0; nAVGfiles = 0; nARFfiles = 0; nRCDfiles = 0;
display_option = 0; filter_option = 0; conversion_factor = 0; conversion_option = 0; verbose = 0;
reject = 0; notchfilter_option = 0; single_trial_option = 0; patient_stem = ''; eeglabpath = 0;

ADunits = 16384; adunitsPerUvolt = 1; nCalCols = 4; plotPauseSec = .5; sgc_version = 'sgc_v2.0';
epoch_pts = zeros(2,1); epoch_ms = zeros(2,1); base_pts = zeros(2,1); base_ms = zeros(2,1);
reject_code = {'ppa','flat','rms'}; scalp_montage = 'Cap34_ChinRef';

%% Initialize 3D locations of scalp electrodes for 34 channel montage

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

%% Parse the variable argument list

for j = 1:length(varargin)
    arg_string = lower(char(varargin{j}));
    if(strfind(arg_string,'.eeg'))
        arg_pointer(j) = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif(strfind(arg_string,'.sgc'))
        arg_pointer(j) = 2;
        nSGCfiles = nSGCfiles + 1;
    elseif(strfind(arg_string,'.mon'))
        arg_pointer(j) = 3;
        nMONfiles = nMONfiles + 1;
    elseif(strfind(arg_string,'.log'))
        arg_pointer(j) = 4;
        nLOGfiles = nLOGfiles + 1;
    elseif(strfind(arg_string,'.avg') | strfind(arg_string,'.gav'))
        arg_pointer(j) = 5;
        nAVGfiles = nAVGfiles + 1;
    elseif(strfind(arg_string,'.arf'))
        arg_pointer(j) = 6;
        nARFfiles = nARFfiles + 1;
    elseif(strfind(arg_string,'.rcd'))
        arg_pointer(j) = 7;
        nRCDfiles = nRCDfiles + 1;
    elseif(strfind(arg_string,'-d'))
        display_option = 1;
        temp = char(varargin{j});
%         display_channel = fix(str2double(temp(strfind(temp,'-d')+2:end)));
        display_channel = fix(str2num(temp(strfind(temp,'-d')+2:end)));
        if(isempty(display_channel)) display_channel = 1; end
        nDisplayChannels = length(display_channel);
    elseif(strfind(arg_string,'-f'))
        filter_option = 1;
        temp = arg_string;
%         filter_values = str2double(temp(strfind(temp,'-f')+2:end));
        filter_values = str2num(temp(strfind(temp,'-f')+2:end));
        lowpass_value = filter_values(1); highpass_value = filter_values(2);
    elseif(strfind(arg_string,'-n'))
        notchfilter_option = 1;
    elseif(strfind(arg_string,'-s'))
        single_trial_option = 1;
        temp = char(varargin{j});
        patient_stem = strcat(temp(3:end),'_');
    elseif(strfind(arg_string,'-c'))
        conversion_option = 1;
        temp = arg_string;
        conversion_factor = str2double(temp(strfind(temp,'-c')+2:end));
        if(isempty(conversion_factor)) conversion_factor = 1; end
    elseif(strfind(arg_string,'-v'))
        verbose = 1;
    else
        fprintf('\n Invalid argument: %s\n',arg_string);
    end
end

indexEEG = find(arg_pointer == 1);
indexSGC = find(arg_pointer == 2);
indexMON = find(arg_pointer == 3);
indexLOG = find(arg_pointer == 4);
indexAVG = find(arg_pointer == 5);
indexARF = find(arg_pointer == 6);
indexRCD = find(arg_pointer == 7);

%scaling values for -d optionadd
if (conversion_option ==1)  %-c option
    Ymax = 500;
    Ymin = -Ymax;
else
    Ymax = 500;
    Ymin = -Ymax;
end

%% Check to make sure we have the necessary files to continue

if(~nEEGfiles) error('No EEG file was specified. Program cannot proceed'); end

if (~nSGCfiles)
    error('No SGC file was specified. Program cannot proceed');
elseif (nSGCfiles > 1 && nSGCfiles ~= nEEGfiles)
    error('Insufficient SGC files specified to match EEG files. Program cannot proceed');
end

if(~nMONfiles)
    fprintf('\nNo MON file was specified. Channels will be numerically labeled\n');
    montage_name = '';
elseif (nMONfiles == 1)
    if (~exist(varargin{indexMON}))
        error('Specified montage file does not exist. Program cannot proceed');
    else
        fidMON = fopen(varargin{indexMON});
        [montageLabels,montageCount] = getstrs(fidMON);
        fclose(fidMON);
    end
    [pathstr, montage_name] = fileparts(varargin{indexMON});
elseif (nMONfiles > 1)
    error('More than one montage file was specified. Program cannot proceed');
end

if(~nAVGfiles)
    fprintf('\nNo AVG or GAV file was specified. Creating MipAvg.avg in working directory\n');
    fidAVG = fopen('MipAvg.avg','wb');
    fidHDR = fopen('MipAvg.hdr','w');
elseif (nAVGfiles == 1)
    if (~exist(varargin{indexAVG},'file'))
        fprintf('Specified AVG or GAV file does not currently exist: Creating file\n');
    else
        fprintf('Specified AVG or GAV file already exists: Overwriting old file\n');
    end
    fidAVG = fopen(varargin{indexAVG},'w');
    hdrFile = varargin{indexAVG};
    if(strfind(hdrFile,'.avg')) hdrFile(strfind(hdrFile,'.avg'):end) = '.hdr'; end
    if(strfind(hdrFile,'.gav')) hdrFile(strfind(hdrFile,'.gav'):end) = '.hdr'; end
    fidHDR = fopen(hdrFile,'w');
elseif (nAVGfiles > 1)
    error('More than one AVG or GAV file was specified. Program cannot proceed');
end

if(~nLOGfiles)
    fprintf('\nNo LOG file was specified. Creating MipAvg.log in working directory\n');
    fidLOG = fopen('MipAvg.log','w');
elseif(nLOGfiles == 1)
    fidLOG = fopen(varargin{indexLOG},'w');
else
    error('More than one LOG file specified. Program cannot proceed');
end

if(~nARFfiles) %% IF no arf file, use variance on channel specified in scg file
    fprintf('\nNo ARF file was specified. No artifact rejection will be applied.\n');
elseif(nARFfiles > 1)
    error('More than one ARF file specified. Program cannot proceed');
end

%% Initialize the log file

fprintf(fidLOG,'\n MIPAVG3 Program executed on %s\n\n EEG files specified:\n',datestr(now));
for j = 1:length(indexEEG)
    fprintf(fidLOG,' %s\n',varargin{indexEEG(j)});
end

fprintf(fidLOG,'\n AVG or GAV file specified:\n');
if (nAVGfiles)
    fprintf(fidLOG,' %s\n',varargin{indexAVG});
else
    fprintf(fidLOG,' No AVG or GAV file was specified. MipAVG.avg will be created in working directory\n');
end

if(nRCDfiles)
    fidRCD = fopen(varargin{indexRCD},'r');
    nrecode = 0;
    while ~feof(fidRCD)
        rcd_string = strtrim(fgetl(fidRCD));
        if(isempty(rcd_string)) break; end
        nrecode = nrecode + 1;
    end
    frewind(fidRCD);
    recode = zeros(nrecode,2);
    for n = 1:nrecode
        rcd_string = fgetl(fidRCD);
        [old_code,new_code] = strread(rcd_string,'%d %d');
        recode(n,1) = old_code;
        recode(n,2) = new_code;
    end
    fclose(fidRCD);
    fprintf(fidLOG,'\nRecode file specified: %s\nNumber of recodes specified = %d\n',char(varargin{indexRCD}),nrecode);
end

if (nMONfiles)
    fprintf(fidLOG,'\n MON files specified:\n');
    fprintf(fidLOG,' %s\n\n',varargin{indexMON});
    for j = 1:length(montageLabels)
        fprintf(fidLOG,' %3d. %s\n',j,montageLabels{j});
    end
else
    fprintf(fidLOG,'\n No MON file was specified. Channels will be numerically labeled.\n');
end
if (conversion_option)
    fprintf(fidLOG,'\n Calibration override specified. Microvolt conversion factor set equal to %5.2f\n',conversion_factor);
end

%% Loop for EEG files: First make sure each exists, extract critical information from the header

nChannels = zeros(nEEGfiles,1); freq = zeros(nEEGfiles,1);
trigPoints = cell(nEEGfiles,1); trigCodes = cell(nEEGfiles,1);
total_points = zeros(nEEGfiles,1);

for j = 1:length(indexEEG)
    
    if(~exist(varargin{indexEEG(j)}))
        fprintf('\nEEG file does not exist --> %s\n',varargin{indexEEG(j)});
        error('Program cannot proceed');
    end
    
    eeg = MipRead2(varargin{indexEEG(j)},0);
    if(~strfind(eeg.DAPinfo.DAPname,'4200a/526'))
        fprintf(fidLOG,'\nUnexpected DAP processor - %s - check A/D units and calibrations\n',eeg.DAPinfo.DAPName);
        fprintf('\nUnexpected DAP processor - %s - check A/D units and calibrations\n',eeg.DAPinfo.DAPName);
    end
    
    nChannels(j) = eeg.nChannels; freq(j) = eeg.design.freq;
    trigChan = squeeze(eeg.data);
    total_points(j) = length(trigChan);
    points = find(trigChan > 0 & trigChan < 251);
    codes = trigChan(points);
    trigCodes{j} = squeeze(codes); trigPoints{j} = squeeze(points);
    
end

%% Make a list of the original codes and print if in verbose mode

old_code = 0;
for j = 1:length(indexEEG)
    old_code = length(trigCodes{j}) + old_code;
end
if(verbose)
    fprintf(fidLOG,'\nNumber of original codes stored in EEG file = %d',old_code);
    code_num = 0;
    for j = 1:length(indexEEG)
        for n = 1:length(trigCodes{j})
            code_num = code_num + 1;
            fprintf(fidLOG,'\n%3d %3d',code_num,trigCodes{j}(n));
        end
    end
    fprintf(fidLOG,'\n\n');
end

%% If recoding is requested, then do it here - check for errors

if(nRCDfiles)
    fprintf(fidLOG,'\nRecoding section\n');
    if(old_code ~= nrecode)
        fprintf(fidLOG,'\nMismatch between number of original and recodes: %d %d',old_code,nrecode);
        fclose('all');
        error('Mismatch between number of original and recodes - Program cannot proceed');
    else
        new_code = 0;
        fprintf(fidLOG,'\n Seq Old  Old  New\n');
        for j = 1:length(indexEEG)
            for n = 1:length(trigCodes{j})
                new_code = new_code + 1;
                fprintf(fidLOG,'\n%4d %4d %4d %4d',new_code,trigCodes{j}(n),recode(new_code,1),recode(new_code,2));
                if(trigCodes{j}(n) ~= recode(new_code,1))
                    fprintf(fidLOG,'\nError in recoding - code mismatch\n');
                    fclose('all');
                    error('Error in recoding - code mismatch. Program cannot proceed');
                else
                    trigCodes{j}(n) = recode(new_code,2);
                end
            end
        end
        fprintf(fidLOG,'\n\n');
    end
end

%% Make sure that all of EEG files have the same number of channels and same sampling rate

if(length(unique(nChannels)) ~= 1 || length(unique(freq)) ~= 1)
    fprintf('\nAn inconsistency was detected in the acquisition parameters of the EEG files\n');
    error('Program cannot proceed');
end

%% Loop for SGC files: First make sure each exists, extract critical information, test epoch limits

epoch_limits = cell(nSGCfiles,1); base_limits = cell(nSGCfiles,1); nBinsSGC = zeros(nSGCfiles,1);
codeMapSGC = cell(nSGCfiles,1); labelsSGC = cell(nSGCfiles,1); nChan = nChannels(1);

for j = 1:length(indexSGC)
    if(~exist(varargin{indexSGC(j)}))
        fprintf('\n SGC file does not exist --> %s\n',varargin{indexSGC(j)});
        error('Program cannot proceed');
    else
        fid = fopen(varargin{indexSGC(j)});
        version = lower(strtrim(fscanf(fid,'%s',1)));
        if(~strcmp(version,sgc_version))
            fprintf('\n SGC file is not current version\n');
            error('Program cannot proceed');
        end
        epoch = fscanf(fid,'%d%d',2);
        base  = fscanf(fid,'%d%d',2);
        first_pt_msec = fscanf(fid,'%d',1);
        nBinsSGC(j) = fscanf(fid,'%d',1);
        codes = fscanf(fid,'%d%d\n',[2,nBinsSGC(j)])';
        [labels,labelsCount] = getstrs(fid);
        labels = deblank(labels);
        for k = labelsCount:-1:1
            if(isempty(labels{k}))
                labels(k) = [];
            else
                break;
            end
        end
        epoch_limits{j} = epoch; base_limits{j} = base; codesSGC{j} = codes; labelsSGC{j} = labels;
    end
    fclose(fid);
end

clear epoch base codes labels;

% Make sure that the epoch limits are the same in all SGC files specified

fprintf(fidLOG,'\n SGC files specified:\n');
for j = 1:length(indexSGC)
    fprintf(fidLOG,' %s\n',varargin{indexSGC(j)});
end
if(length(unique(cat(1,epoch_limits{:}))) ~= 2)
    fprintf('\nAn inconsistency was detected in the epoch limits within the SGC files\n');
    error('Program cannot proceed');
end

msecpt = 1000/freq(1);
epoch_ms = epoch_limits{1};
base_ms = base_limits{1};
epoch_pts(1) = round(epoch_ms(1) / msecpt);
epoch_pts(2) = round(epoch_ms(2) / msecpt);
base_pts(1) = round(base_ms(1) / msecpt);
base_pts(2) = round(base_ms(2) / msecpt);
npts = epoch_pts(2) - epoch_pts(1) + 1;

fprintf(fidLOG,'\nSegmentation Control File\n\nThe designated epoch is specified relative to the stimulus code\n');
fprintf(fidLOG,'Epoch begin and end (msec): %d %d\n',epoch_ms(1),epoch_ms(2));
fprintf(fidLOG,'Epoch begin and end (pts) : %d %d\n',epoch_pts(1),epoch_pts(2));
fprintf(fidLOG,'Number of points in epoch : %d\n',npts);
fprintf(fidLOG,'\nBaseline epoch specified relative to the stimulus code\n');
fprintf(fidLOG,'\nBase begin and end (msec) : %d %d\n',base_ms(1),base_ms(2));
fprintf(fidLOG,'Base begin and end (pts)  : %d %d\n',base_pts(1),base_pts(2));

if(epoch_pts(1) >= epoch_pts(2) || base_pts(1) >= base_pts(2) || base_pts(1) < epoch_pts(1) || base_pts(2) > epoch_pts(2))
    fprintf(fidLOG,'\nInconsistency in designated epochs - program cannot continue\n');
    fclose('all');
    error('Inconsistency in designated epochs - program cannot continue');
end
base_ms(1) = base_ms(1) - epoch_ms(1);
base_ms(2) = base_ms(2) - epoch_ms(1);
base_pts(1) = base_pts(1) - epoch_pts(1) + 1;
base_pts(2) = base_pts(2) - epoch_pts(1) + 1;
fprintf(fidLOG,'\nBaseline epoch specified relative to the designated epoch\n');
fprintf(fidLOG,'\nBase begin and end (msec) : %d %d\n',base_ms(1),base_ms(2));
fprintf(fidLOG,'Base begin and end (pts)  : %d %d\n',base_pts(1),base_pts(2));

%% Get critical information from ARF file

rms_check = 0; flat_check = 0; ppa_check = 0; any_check = 0;

if(nARFfiles)
    
    rejCrit = zeros(eeg.nChannels, 12);
    fid = fopen(varargin{indexARF(1)});
    
    while ~feof(fid)
        arf_string = strtrim(fgetl(fid));
        [channel,typeCrit,value1,value2,value3] = strread(arf_string,'%d %s %d %d %d');
        if(isempty(value2))
            value2 = 1;
        else
            value2 = round(value2 / msecpt) - epoch_pts(1) + 1;
        end
        if(isempty(value3))
            value3 = npts;
        else
            value3 = round(value3 / msecpt) - epoch_pts(1) + 1;
        end
        if (value2 > npts || value3 > npts)
            fprintf(fidLOG,'\nARF temporal range exceeds raw epoch limits (pts): %d %d\n',value2,value3);
            fclose('all');
            error('Error in ARF specified artifact rejection range - Program cannot proceed');
        end
        
        if(channel <= eeg.nChannels && ~isempty(value1))
            
            switch lower(char(typeCrit))
                case 'ppa'
                    rejCrit(channel,1) = 1;
                    rejCrit(channel,2) = value1;
                    rejCrit(channel,3) = value2;
                    rejCrit(channel,4) = value3;
                    ppa_check = 1; any_check = 1;
                    
                case 'flat'
                    rejCrit(channel,5) = 1;
                    rejCrit(channel,6) = value1;
                    rejCrit(channel,7) = value2;
                    rejCrit(channel,8) = value3;
                    flat_check = 1; any_check = 1;
                    
                case 'rms'
                    rejCrit(channel,9) = 1;
                    rejCrit(channel,10) = value1;
                    rejCrit(channel,11) = value2;
                    rejCrit(channel,12) = value3;
                    rms_check = 1; any_check = 1;
                    
                otherwise
                    fprintf(fidLOG,'\n Unknown specifier in ARF file: %s\n',typeCrit);
            end
        end
    end
    fclose(fid);
end

if(any_check)
    fprintf(fidLOG,'\nARF file specified: %s\n',varargin{indexARF});
    fprintf(fidLOG,'\nThe following flags are set:\n rms_check: %d\n ppa_check: %d\n flat_check: %d\n',rms_check,ppa_check,flat_check);
    
    % print artifact rejection criteria to log file
    
    fprintf(fidLOG,'\nArtifact rejection option set for \nchan\tfunc\tthreshold\tstart\tend \n');
    for j=1:nChan
        if (sum(rejCrit(j,:)) >0)
            if (rejCrit(j,1)>0)
                fprintf(fidLOG,'\n%4d\t ppa\t %d\t%d\t%d\t',j,rejCrit(j,2),rejCrit(j,3),rejCrit(j,4));
            end;
            if (rejCrit(j,5)>0)
                fprintf(fidLOG,'\n%4d\t flat\t %d\t%d\t%d\t',j,rejCrit(j,6),rejCrit(j,7),rejCrit(j,8));
            end;
            if (rejCrit(j,9)>0)
                fprintf(fidLOG,'\n%4d\t rms\t %d\t%d\t%d\t',j,rejCrit(j,10),rejCrit(j,11),rejCrit(j,12));
            end;
        end
    end
    
    fprintf(fidLOG,'\n');
    
end
if (~nARFfiles) fprintf(fidLOG,'\nNo ARF file was specified.\n'); end

%% All SGC and EEG files exist and have no critical inconsistencies
% Make sure that all codes have a complete epoch and do not extend past either end of data array. If so, set trigger code to zero.

for j = 1:nEEGfiles
    for k = 1:length(trigPoints{j})
        if(trigPoints{j}(k) + epoch_pts(1) < 0 || trigPoints{j}(k) + epoch_pts(2) > total_points(j))
            fprintf(fidLOG,'\nTrigger code (not remapped) %d at Trigger point %d for EEGfile %d exceeds total points limits\n',trigCodes{j}(k),trigPoints{j}(k),j);
            trigCodes{j}(k) = 0;
        else
            break;
        end
    end
    for k = length(trigPoints{j}):-1:1
        if(trigPoints{j}(k) + epoch_pts(1) < 0 || trigPoints{j}(k) + epoch_pts(2) > total_points(j))
            fprintf(fidLOG,'\nTrigger code (not remapped) %d at Trigger point %d for EEGfile %d exceeds total points limits\n',trigCodes{j}(k),trigPoints{j}(k),j);
            trigCodes{j}(k) = 0;
        else
            break;
        end
    end
end

%% Remap the trigCodes for each EEG file according to the code mapping in the corresponding SGC file

for j = 1:nEEGfiles
    codes = trigCodes{j};
    points = trigPoints{j};
    mapped_codes = zeros(1,length(codes));
    if(nSGCfiles == 1) mapSGC = codesSGC{1}; else mapSGC = codesSGC{j}; end
    [nlines,nmp] = size(mapSGC);
    %   for k = 1:length(mapSGC)
    for k = 1:nlines
        mapped_codes(find(codes == mapSGC(k,1))) = mapSGC(k,2);
    end
    trigPoints{j} = points(find(mapped_codes>0));
    trigCodes{j} = mapped_codes(find(mapped_codes>0));
    lastTime = 0;
    fprintf(fidLOG,'\n EEG file - > %s\n Code  Code_Time  Event_Time  Elapsed_Time\n\n',varargin{indexEEG(j)});
    for k = 1:length(trigPoints{j})
        time1 = trigPoints{j}(k) * msecpt;
        time2 = (trigPoints{j}(k) + epoch_pts(1)) * msecpt;
        elpsdTime = time2 - lastTime;
        lastTime = time2;
        fprintf(fidLOG,'%3d %10d %10d %10d\n',trigCodes{j}(k),time1,time2,elpsdTime);
    end
end

clear points codes mapped_codes mapSGC;

%% Count the number of unique trigger codes following remapping

unique_codes = unique(cat(2,trigCodes{:}));
code_count = zeros(length(unique_codes),1);
for j = 1:length(unique_codes)
    code_count(j) = length(find(cat(2,trigCodes{:}) == unique_codes(j)));
end

fprintf(fidLOG,'\n Number of unique codes -> %d\n',length(unique_codes));
for j = 1:length(unique_codes)
    fprintf(fidLOG,' Code %d : %d\n',unique_codes(j),code_count(j));
end
code_actual = code_count;

%% Create memory arrays for each unique bin, we can assume that the epoch limits and freq are consistent

avg = zeros(nChan+1,npts,length(unique_codes));

%% If single trial_option is set, open binary files to write single trial data

if(single_trial_option)
    fprintf(fidLOG,' \n\n Single Trial Option specified:\n\n');
    fid_single_trial = zeros(length(unique_codes),1);
    for k = 1:length(unique_codes)
        fid_single_trial(k) = fopen(sprintf('single_trial_code%03d.seg',unique_codes(k)),'w+');
        fprintf(fidLOG,'Code %03d stored in %s\n',unique_codes(k),sprintf('single_trial_code%03d.seg',unique_codes(k)));
    end
end

%% If no montage file was supplied, make sure to create a numeric channel list

if(~exist('montageLabels','var'))
    montageLabels = cell(nChan+1,1);
    for k = 1:nChan
        montageLabels{k} = sprintf('chan%03d',k);
    end
    montageLabels{nChan+1} = 'trigger';
else
    if(length(montageLabels) > nChan+1)
        fprintf(fidLOG,' \n\n **** WARNING **** Too many channels were specified in MON file.\n\n');
    elseif(length(montageLabels) < nChan+1)
        fprintf(fidLOG,' \n\n **** WARNING **** Too few channels were specified in MON file. Compensating!\n\n');
        for k = length(montageLabels)+1:nChan
            montageLabels{k} = sprintf('chan%03d',k);
        end
        montageLabels{nChan+1} = 'trigger';
    end
end

%% Create a display if display option is set

if(display_option)
    display_channel(find(display_channel > nChan+1)) = 1;
    nHandles = zeros(nDisplayChannels,1);
    figureHandle = figure;
    nrows = fix(sqrt(nDisplayChannels));
    ncols = round((nDisplayChannels/nrows));
    Ytemp = zeros(npts,1);
    Xtemp = [epoch_ms(1):msecpt:epoch_ms(2)];
    Xmin = epoch_ms(1); Xmax = epoch_ms(2);
    for k = 1:nDisplayChannels
        nHandles(k) = subplot(nrows,ncols,k);
        plot(Xtemp,Ytemp);
        if(display_channel(k) == nChan+1)
            axis([Xmin, Xmax, 0, max(unique_codes)]);
        else
            axis([Xmin, Xmax, Ymin, Ymax]);
        end
        set(gca,'YTick',[Ymin, 0 Ymax],'Xtick',[0,Xmax],'FontSize',7);
        title(montageLabels{display_channel(k)},'FontSize',10);
    end
    clear Xtemp Ytemp;
end

%% Begin looping through the EEG files for averaging


for j = 1:length(indexEEG)
    
    % get conversion data
    
    fprintf(fidLOG,'\n EEG file - > %s\n microVolt Conversion Factors:\n\n',varargin{indexEEG(j)});
    eeg = MipRead2(varargin{indexEEG(j)},1);
    uvcf(1:eeg.nChannels) = (eeg.calInfo.pp(1:eeg.nChannels) ./ eeg.calInfo.ppv) .* ADunits;
    uvcf(1:eeg.nChannels) = uvcf(1:eeg.nChannels) ./ eeg.calInfo.uVolts;
    maxlines = fix(eeg.nChannels / nCalCols);
    remainingCals = eeg.nChannels - (maxlines * nCalCols);
    for k = 1:maxlines
        fprintf(fidLOG,' Channel %3i = %.3f\t Channel %3i = %.3f\t Channel %3i = %.3f\t Channel %3i = %.3f\n',...
            k,uvcf(k),k+maxlines,uvcf(k+maxlines),k+(maxlines*2),uvcf(k+maxlines*2),k+(maxlines*3),uvcf(k+maxlines*3));
    end
    for k = 1:remainingCals
        offset = (maxlines * nCalCols) + k;
        fprintf(fidLOG,'\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t Channel %3i = %.3f\n',offset,uvcf(offset));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conversion to uVolts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(conversion_option) uvcf(:) = conversion_factor; end
    
    for k = 1:eeg.nChannels
        divisor = 4 * uvcf(k) * (1/adunitsPerUvolt);        %the division by 4 downshifts from 16-bits to 14-bits as per A/D
        eeg.data(k,:) = eeg.data(k,:) ./ divisor;
    end
    
    if(notchfilter_option)
        fprintf(fidLOG,'\nApplying notch filter\n');
        for k = 1:eeg.nChannels
            eeg.data(k,:) = notch60(squeeze(eeg.data(k,:)),freq(j));
        end
    end
    
    % Filter the EEG data if option is set
    
    if (filter_option)
        fprintf(fidLOG,'\nApplying bandpass filter (low high): %6.2f f6.2f\n',lowpass_value, highpass_value);
        for k = 1:eeg.nChannels
            eeg.data(k,:) = bandpass_EEG(squeeze(eeg.data(k,:)),lowpass_value, highpass_value,freq(j));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % artifact rejection & display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Sum the data by bin to the avg array
    % Because the display and artifact rejection options are compute intensive, branch to avoid wasting time
    
    for m = 1:length(trigCodes{j}(:))
        
        reject=0;
        first_pt = trigPoints{j}(m) + epoch_pts(1);
        last_pt =  trigPoints{j}(m) + epoch_pts(2);
        nbin = find(unique_codes == trigCodes{j}(m));
        single_trial = eeg.data(1:eeg.nChannels+1,first_pt:last_pt);
        base = mean(single_trial');
        single_trial = single_trial - repmat(base',[1,npts]);
        
        % execute next code block if no artifact rejection selected
        
        if(~any_check)
            
            avg(:,:,nbin) = avg(:,:,nbin) + single_trial;
            if(single_trial_option)
                fwrite(fid_single_trial(nbin),single_trial,'double');
            end
            % compute min/max and variance of each single trials over entire epoch for quality check
            
            ppa_all = max(single_trial') - min(single_trial');
            rms_all = var(single_trial');
            ppa_max_chan = find(ppa_all == max(ppa_all));
            rms_max_chan = find(rms_all == max(rms_all));
            
            fprintf(fidLOG,' %3d\t %3d %6.0f\t %3d %6.0f\n',m, rms_max_chan, rms_all(rms_max_chan), ppa_max_chan, ppa_all(ppa_max_chan));
            
            if(display_option)
                rms_vals = zeros(nDisplayChannels,1);
                for nplot=1:nDisplayChannels
                    set(get(nHandles(nplot),'Children'),'YData',squeeze(single_trial(display_channel(nplot),1:npts)));
                    rms_vals(nplot) = rms_all(display_channel(nplot));
                end
                set(figureHandle,'Name',[char(labelsSGC{1}(nbin)),' TotalEpochVar: ',sprintf('%5.0f ',rms_vals)]);
                pause;
            end
        end
        
        % execute next code block if any artifact rejection option was selected
        
        if(any_check)
            
            % If no filtering was previously specified for the EEG, then apply the notch filter before testing for artifacts
            
            if(~notchfilter_option && ~filter_option)
                for p = 1:eeg.nChannels
                    single_trial(p,:) = notch60(squeeze(single_trial(p,:)),freq(j));
                end
            end
            
            % Test for artifact rejection - stop after detecting first rejection event
            
            channr = 0;
            reject = 0;
            rejValue = 0;
            
            %ARF file indicated rms analysis on channels specified
            if(rms_check)
                for chans = 1:eeg.nChannels
                    if (rejCrit(chans,9))
                        rmsReject = var(squeeze(single_trial(chans,rejCrit(chans,11):rejCrit(chans,12)))');
                        if (rmsReject >= rejCrit(chans,10));
                            reject = 3;
                            channr = chans;
                            rejValue = rmsReject;
                        end
                    end
                    if(reject) break; end
                end
            end
            
            % ARF file indicates analysis by ppa for this channel
            if(ppa_check && ~reject)
                for chans = 1:eeg.nChannels
                    if (rejCrit(chans,1))
                        [min_ppa,max_ppa] = minmax(squeeze(single_trial(chans,rejCrit(chans,3):rejCrit(chans,4)))');
                        ppaReject = max_ppa - min_ppa;
                        if(ppaReject >= rejCrit(chans,2))
                            reject = 1;
                            channr = chans;
                            rejValue = ppaReject;
                        end
                    end
                    if(reject) break; end
                end
            end
            
            % ARF file indicates analysis for flatline for this channel
            if(flat_check && ~reject)
                for chans = 1:eeg.nChannels
                    if(rejCrit(chans, 5))
                        flatW = rejCrit(chans, 4);
                        flatProg = rejCrit(chans, 5);
                        flatEnvelope= rejCrit(chans, 6);
                        [reject, rejValue] = flatline_check(single_trial, npts, chans, flatW, flatProg, flatEnvelope);
                        
                        if(reject)
                            channr = chans;
                            break;
                        end
                    end
                end
            end
            
            if (reject)
                rejectType = char(reject_code{reject});
                code_actual(nbin) = code_actual(nbin) - 1;
                fprintf(fidLOG,'Trial: %3d Reject Channel: %3d\t %s\t%6.0f\n',m, channr, rejectType, rejValue);
            else
                if(~notchfilter_option && ~filter_option)
                    single_trial = eeg.data(1:eeg.nChannels+1,first_pt:last_pt);
                    base = mean(single_trial');
                    single_trial = single_trial - repmat(base',[1,npts]);
                end
                avg(:,:,nbin) = avg(:,:,nbin) + single_trial;
                if(single_trial_option)
                    fwrite(fid_single_trial(nbin),single_trial,'double');
                end
                
            end
            
            if(display_option)
                for nplot=1:nDisplayChannels
                    set(get(nHandles(nplot),'Children'),'YData',single_trial(display_channel(nplot),1:npts));
                end
                if (reject)
                    set(figureHandle, 'Name',[char(labelsSGC{1}(nbin)),sprintf('  Reject %s on Channel %03d = %6.0f',rejectType,channr,rejValue)]);
                else
                    set(figureHandle, 'Name',[char(labelsSGC{1}(nbin)),sprintf('  Accept')]);
                end
                pause;
            end
        end
    end
end

%% Average and subtract prestimulus baseline

for k = 1:length(unique_codes)
    avg(:,:,k) = avg(:,:,k) ./ code_actual(k);
end

baseline = zeros(length(unique_codes),nChan+1);

for k = 1:length(unique_codes)
    baseline(k,:) = mean(avg(:,base_pts(1):base_pts(2),k)');
end
for k = 1:length(unique_codes)
    for m = 1:nChan
        avg(m,:,k) = avg(m,:,k) - baseline(k,m);
    end
    % Divide the trigger channel by 10 to keep it close to EEG range
    avg(nChan+1,:,k) = avg(nChan+1,:,k) / 10;
end

%% Write avg file

for k = 1:length(unique_codes)
    fwrite(fidAVG,avg(:,:,k)','float32');
end

%% Write HDR file (description below is taken from EEGAD documentation for HDR)

% Line 1: expName is the experiment name
% Line 2: expDate is the experiment date string
% Line 3: nChannels is the number of data channels (electrodes).
% Line 4: nPoints is the number of data points collected per channel.
% Line 5: sampling is the sampling rate of the data points in ms/point.
% Line 6: uvunits is the microvolt conversion factor in raw data units/microvolt.
% Line 7: onset is the stimulus onset from the first data point in ms.
% binNames is a cell array of bin names.
% chanNames is a cell array of channel names.
% coords is a matrix of electrode coordinates with size [nChannels 2] or [nChannels 3].

labelsCount = length(labelsSGC{end});
if (length(unique_codes) ~= labelsCount)
    fprintf(fidLOG,' \n\n **** WARNING **** Number of unique codes does not match number of labels. Compensating!\n\n');
end
if (length(unique_codes) > labelsCount)
    for k = labelsCount+1:unique_codes
        labelsSGC{indexSGC(end)}(k) = sprintf('label%03d',k);
    end
end

fprintf(fidHDR,'%s\n',eeg.runInfo.comments(1,1:end));
fprintf(fidHDR,'%s\n',eeg.runInfo.time);
%fprintf(fidHDR,'%d\n%d\n%d\n%d\n',nChan+1,npts,msecpt,adunitsPerUvolt,(base_pts(2)-1)*msecpt);
fprintf(fidHDR,'%d\n%d\n%d\n%d\n',nChan+1,npts,msecpt,adunitsPerUvolt,-first_pt_msec);
for k = 1:length(unique_codes)
    label = char(labelsSGC{end}(k));
    fprintf(fidHDR,'%s\n',label);
end
for k = 1:nChan+1, fprintf(fidHDR,'%s\n',montageLabels{k}); end

fprintf(fidLOG,'\n Label Total Number Total Accepted Percent Averaged\n');
for j = 1:length(unique_codes)
    fprintf(fidLOG,' Code %d : (%s) %d %d %4.2f\n',unique_codes(j),char(labelsSGC{end}(j)),code_count(j),...
        code_actual(j),code_actual(j)/code_count(j));
end

fclose('all');

%% if single_trial option (-s) is set, open single trial files,
% read contents, and create a eeglab .set file for each condition
% This next section was substantially modified by GM on 12-Dec-2009

if(~single_trial_option)
    fclose('all');
    toc
    return
end

[sx,sy] = size(single_trial);
single_trial_count = sx*sy;

% Before we go further, make sure that the required eeglab functions can be found

if(exist('pop_importdata.m','file') ~= 2)
    if(exist('/usr/local/packages/MATLABPackages/eeglab/functions/popfunc','dir') == 7)
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/popfunc');
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/adminfunc');
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/guifunc');
        eeglabpath = 1;
    else
        fprintf('\nRequired eeglab functions cannot be found - check you paths!\n');
        fclose('all');
        return
    end
end

% if we are using the standard 34 channel scalp montage, we drop the last channel

if strcmpi(montage_name,scalp_montage) nChan = nChan - 1; end

n_secpt = msecpt/1000;

for k = 1:length(unique_codes)
    
    infile = sprintf('single_trial_code%03d.seg',unique_codes(k));
    fid_temp = fopen(infile,'r');
    single_trials = zeros(sx,sy,code_actual(k));
    label = strtrim(char(labelsSGC{end}(k)));
    p = strfind(label,' ');
    for q = 1:length(p) label(p(q)) = '_'; end
    outname = sprintf('%s%s_code%03d',patient_stem,label,unique_codes(k));
    [outpath,outname] = fileparts(outname);
    if(strcmp(outpath,''))
        outpath = pwd;
    end
    
    for n = 1:code_actual(k)
        eofstat = feof(fid_temp);
        if(~eofstat)
            single_trial = reshape(fread(fid_temp,single_trial_count,'double'),sx,sy);
        else
            fprintf('\nPremature end of file encountered after %d single trials in %s \n',...
                n,sprintf('single_trial_code%03d.seg',unique_codes(k)));
        end
        single_trials(:,:,n) = single_trial(:,:);
    end
    
    chanloc = struct('theta',cell(1,nChan),'radius',cell(1,nChan),'labels',cell(1,nChan),...
        'sph_theta',cell(1,nChan),'sph_phi',cell(1,nChan),'sph_radius',cell(1,nChan),...
        'X',cell(1,nChan),'Y',cell(1,nChan),'Z',cell(1,nChan),'type',cell(1,nChan));
    
    
    nrow = ceil(sqrt(nChan));
    nrow_center = (nrow + 1) / 2;
    ncol = ceil(nChan/nrow);
    ncol_center = (ncol + 1) / 2;
    
    if ~(strcmpi(montage_name,scalp_montage))
        for n = 1:nChan
            chanloc(n).theta = 0;
            chanloc(n).radius = 0;
            chanloc(n).labels = char(montageLabels{n});
            chanloc(n).sph_theta = 0;
            chanloc(n).sph_phi = 0;
            chanloc(n).sph_radius = 0;
            nr = rem(n,nrow);
            if (nr == 0) nr = nrow; end
            chanloc(n).X = nr - nrow_center;
            nc = mod(ceil(n/nrow),ncol);
            if (nc == 0) nc = ncol; end
            chanloc(n).Y = nc - ncol_center;
            chanloc(n).Z = 1;
        end
    else
        for n = 1:nChan
            chanloc(n).theta = cap34{n}{2};
            chanloc(n).radius = cap34{n}{3};
            chanloc(n).labels = char(cap34{n}{1});
            chanloc(n).sph_theta = cap34{n}{7};
            chanloc(n).sph_phi = cap34{n}{8};
            chanloc(n).sph_radius = cap34{n}{9};
            chanloc(n).X = cap34{n}{4};
            chanloc(n).Y = cap34{n}{5};
            chanloc(n).Z = cap34{n}{6};
            chanloc(n).type = char(cap34{n}{10});
        end
    end
    
    ntrials = size(single_trials,3);
    
    EEGOUT = pop_importdata('data',single_trials(1:nChan,:,:),'dataformat','array','srate',1/n_secpt,...
        'nbchan',nChan,'xmin',epoch_ms(1)/1000,'pnts',npts,'setname',outname,...
        'chanlocs',chanloc);
    
    EEGOUT.event = struct('type',cell(1,ntrials),'epoch',cell(1,ntrials));
    EEGOUT.avg = avg(1:nChan,:,k);
    
    for n = 1:ntrials
        EEGOUT.event(n).type = label;
        EEGOUT.event(n).epoch = n;
    end
    
    pop_saveset(EEGOUT,'filepath',outpath,'filename',strcat(outname,'.set'),'savemode','onefile');
    
    fclose(fid_temp);
    delete(sprintf('single_trial_code%03d.seg',unique_codes(k)));
    
end
if eeglabpath
    rmpath('/usr/local/packages/MATLABPackages/eeglab/functions/popfunc');
    rmpath('/usr/local/packages/MATLABPackages/eeglab/functions/adminfunc');
    rmpath('/usr/local/packages/MATLABPackages/eeglab/functions/guifunc');
end
fclose('all');
return
end

% Gregory McCarthy, 08/01/02.  First working version of original MIPAVG
% See MIPAVG code for revision history of MIPAVG
% Gregory McCarthy 09/30/07. Begin major revision of MIPAVG2
%
% Assume data was collected with a DAP 4200 and at 14-bit resolution
% Change output of mipavg from integer to float, and change output scale to microvolts
% Fixed bug that caused fatal error when prestimulus epoch was equal to poststimulus epoch
% Added the variance to the banner display with plot option
%
% Gregory McCarthy 03/11/2008. Continue major revision

%***********************************************************
% BEGIN MipRead2 FUNCTION
%***********************************************************

function MIP=MipRead2(varargin)
%MipRead2 Read header and (optionally) data from file created by MIP for Windows.
%
%   MIP=MipRead2(file,data_flag);
%
%   file is the name of the MIP EEG data file.
%   data_flag = 0, then only read header
%   data_flag = 1, then also read data
%
%   MIP is a structure with many fields including the following:
%     nChannels is the number of data channels (electrodes),
%       not including the digital channel.
%     nPoints is the number of data points collected per channel.
%     data is the EEG data in a [nChannels+1 nPoints] matrix.
%       The last row contains the digital channel.

% CVS ID and authorship of this code
% CVSId = '$Id: MIPRead.m,v 1.4 2005/02/21 23:32:46 michelich Exp $';
% CVSRevision = '$Revision: 1.4 $';
% CVSDate = '$Date: 2005/02/21 23:32:46 $';
% CVSRCSFile = '$RCSfile: MIPRead.m,v $';

%TODO: Fix header string decoding

% Defaults and consants
MIPhSz=8088;                   % Header size
MIPsig1 = 10706671;            % Header signature for valid Windows header
MIPsig2 = 10706672;            % Header signature for converted DOS header
MIPver=11;                     % Header version
MIPdSig=1152444886;            % Design parameters signature
MIPcSig=1051667;               % Calibration parameters signature
MIPmaxChans=128;               % Max possible channels
MIPmaxBins=32;                 % Max possible bins (for montage info)

if (isempty(varargin))
    error('No arguments to MipRead2 - program cannot proceed');
elseif (length(varargin) == 1)
    data_flag=1;
else
    data_flag = varargin{2};
end
file = varargin{1};
% Calculate data file size
EEGfid=fopen(file);
if EEGfid==-1
    error(['Couldn''t open data file "' file '"'])
end
fseek(EEGfid,0,'eof');
fsize=ftell(EEGfid);
fseek(EEGfid,0,'bof');

% Read header info (0:8087=8088 bytes)
MIP.sig=fread(EEGfid,1,'int32');                                        % 0
fread(EEGfid,4,'uchar');                                                % 4
if MIP.sig~=MIPsig1 && MIP.sig~=MIPsig2
    error(['"' file '" is not a modern or converted MIP EEG file!'])
end
MIP.ver=fread(EEGfid,1,'float64');                                      % 8
if MIP.ver~=MIPver
    error(['Unrecognized header version in "' file '"!'])
end
MIP.hdrLen=fread(EEGfid,1,'int32');                                     % 16
MIP.dataStart=fread(EEGfid,1,'int32');                                  % 20
MIP.dataFooter=fread(EEGfid,1,'uint32');                                % 24
fread(EEGfid,4,'uchar');                                                % 28
% MIP.DAPinfo (32:263=232 bytes)
MIP.DAPinfo.minChanSep=fread(EEGfid,1,'float64');                       % 32
MIP.DAPinfo.it=fread(EEGfid,1,'float64');                               % 40
MIP.DAPinfo.ot=fread(EEGfid,1,'float64');                               % 48
MIP.DAPinfo.incr=fread(EEGfid,1,'float64');                             % 56
MIP.DAPinfo.oma=fread(EEGfid,1,'float64');                              % 64
MIP.DAPinfo.sioa=fread(EEGfid,1,'float64');                             % 72
MIP.DAPinfo.DAPname=sprintf('%s',char(fread(EEGfid,[1 20],'uchar')));   % 80
MIP.DAPinfo.sysName=sprintf('%s',char(fread(EEGfid,[1 20],'uchar')));   % 100
MIP.DAPinfo.maxVolt=fread(EEGfid,1,'float64');                          % 120
MIP.DAPinfo.minVolt=fread(EEGfid,1,'float64');                          % 128
MIP.DAPinfo.voltGrads=fread(EEGfid,1,'int32');                          % 136
MIP.DAPinfo.DAPchannels=fread(EEGfid,1,'int32');                        % 140
MIP.DAPinfo.timeVal=fread(EEGfid,1,'float64');                          % 144
MIP.DAPinfo.digGrads=fread(EEGfid,1,'int32');                           % 152
MIP.DAPinfo.ic=fread(EEGfid,1,'int32');                                 % 156
MIP.DAPinfo.dc=fread(EEGfid,1,'int32');                                 % 160
MIP.DAPinfo.nc=fread(EEGfid,1,'int32');                                 % 164
MIP.DAPinfo.dtaddr=fread(EEGfid,1,'int32');                             % 168
MIP.DAPinfo.nataddr=fread(EEGfid,1,'int32');                            % 172
MIP.DAPinfo.clockNum=fread(EEGfid,1,'int32');                           % 176
MIP.DAPinfo.andMask=fread(EEGfid,1,'int32');                            % 180
MIP.DAPinfo.orMask=fread(EEGfid,1,'int32');                             % 184
MIP.DAPinfo.xorMask=fread(EEGfid,1,'int32');                            % 188
MIP.DAPinfo.aScale=fread(EEGfid,1,'float64');                           % 192
MIP.DAPinfo.dScale=fread(EEGfid,1,'float64');                           % 200
MIP.DAPinfo.digTxt=fread(EEGfid,1,'uchar');                             % 208
fread(EEGfid,7,'uchar');                                                % 209
MIP.DAPinfo.calUVolt=fread(EEGfid,1,'float64');                         % 216
MIP.DAPinfo.testFreq=fread(EEGfid,1,'float64');                         % 224
MIP.DAPinfo.calFreq=fread(EEGfid,1,'float64');                          % 232
MIP.DAPinfo.natType=sprintf('%s',char(fread(EEGfid,[1 10],'uchar')));   % 240
fread(EEGfid,2,'uchar');                                                % 250
MIP.DAPinfo.readCodes=fread(EEGfid,1,'int32');                          % 252
MIP.DAPinfo.compressCodes=fread(EEGfid,1,'int32');                      % 256
fread(EEGfid,4,'uchar');                                                % 260
% MIP.design (264:2951=2688 bytes)
MIP.design.sig=fread(EEGfid,1,'int32');                                 % 264
if MIP.design.sig~=MIPdSig
    error(['Unrecognized design header in "' file '"!'])
end
MIP.design.seconds=fread(EEGfid,1,'int32');                             % 268
MIP.design.expand=fread(EEGfid,1,'int16');                              % 272
fread(EEGfid,6,'uchar');                                                % 274
MIP.design.freq=fread(EEGfid,1,'float64');                              % 280
MIP.design.idleFreq=fread(EEGfid,1,'float64');                          % 288
MIP.design.oneFreq=fread(EEGfid,1,'float64');                           % 296
MIP.design.logName=sprintf('%s',char(fread(EEGfid,[1 9],'uchar')));     % 304
MIP.design.logExt=sprintf('%s',char(fread(EEGfid,[1 4],'uchar')));      % 313
fread(EEGfid,3,'uchar');                                                % 317
MIP.design.channels=fread(EEGfid,[1 MIPmaxChans],'int32');              % 320
MIP.design.prePts=fread(EEGfid,1,'int32');                              % 832
MIP.design.postPts=fread(EEGfid,1,'int32');                             % 836
MIP.design.binCodes=fread(EEGfid,[1 MIPmaxBins],'uint32');              % 840
MIP.design.binNums=fread(EEGfid,[1 MIPmaxBins],'uint32');               % 968
for n=1:MIPmaxChans                                                     % 1096
    MIP.design.montLabels{n}=sprintf('%s',char(fread(EEGfid,[1 10],'uchar')));
end
MIP.design.montLabels=char(MIP.design.montLabels);
MIP.design.montChans=fread(EEGfid,[1 MIPmaxChans],'int32');             % 2376
MIP.design.montName=sprintf('%s',char(fread(EEGfid,[1 61],'uchar')));   % 2888
MIP.design.avgMode=char(fread(EEGfid,1,'uchar'));                       % 2949
fread(EEGfid,2,'uchar');                                                % 2950
% MIP.runInfo (2952:3399=448 bytes)
MIP.runInfo.time=sprintf('%s',char(fread(EEGfid,[1 26],'uchar')));      % 2952
MIP.runInfo.op=sprintf('%s',char(fread(EEGfid,[1 60],'uchar')));        % 2978
MIP.runInfo.subject=sprintf('%s',char(fread(EEGfid,[1 60],'uchar')));   % 3038
MIP.runInfo.control=sprintf('%s',char(fread(EEGfid,[1 60],'uchar')));   % 3098
for n=1:4                                                               % 3158
    MIP.runInfo.comments{n}=sprintf('%s',char(fread(EEGfid,[1 60],'uchar')));
end
MIP.runInfo.comments=char(MIP.runInfo.comments);
fread(EEGfid,2,'uchar');                                                % 3398
% MIP.calInfo (3400:6039=2640 bytes)
MIP.calInfo.sig=fread(EEGfid,1,'int32');                                % 3400
if MIP.calInfo.sig~=MIPcSig
    error(['Unrecognized calibration header in "' file '"!'])
end
MIP.calInfo.time=sprintf('%s',char(fread(EEGfid,[1 26],'uchar')));      % 3404
fread(EEGfid,2,'uchar');                                                % 3430
MIP.calInfo.nChannels=fread(EEGfid,1,'int32');                          % 3432
MIP.calInfo.channels=fread(EEGfid,[1 MIPmaxChans],'int32');             % 3436
fread(EEGfid,4,'uchar');                                                % 3440
MIP.calInfo.base=fread(EEGfid,[1 MIPmaxChans],'float64');               % 3952
MIP.calInfo.pp=fread(EEGfid,[1 MIPmaxChans],'float64');                 % 4976
MIP.calInfo.uVolts=fread(EEGfid,1,'float64');                           % 6000
MIP.calInfo.voltGrads=fread(EEGfid,1,'int32');                          % 6008
fread(EEGfid,4,'uchar');                                                % 6012
MIP.calInfo.ppv=fread(EEGfid,1,'float64');                              % 6016
MIP.calInfo.freq=fread(EEGfid,1,'float64');                             % 6024
MIP.calInfo.avgVals=fread(EEGfid,1,'int32');                            % 6032
% No need to read final header padding
%fread(EEGfid,4,'uchar');                                                % 6036
% unused space (6040:8087=2048 bytes))
%fread(EEGfid,2048,'uchar');                                             % 6040

% Check parameters
if (isempty(MIP.hdrLen) | MIP.hdrLen~=MIPhSz)
    fclose(EEGfid);
    error(['Invalid header length in "' file '"']);
end
if (isempty(MIP.dataStart) | MIP.dataStart~=MIP.hdrLen)
    fclose(EEGfid);
    error(['Invalid pointer to data start in "' file '"']);
end
if (isempty(MIP.dataFooter) | MIP.dataFooter<MIP.hdrLen | MIP.dataFooter>fsize)
    %    fclose(EEGfid);
    %    error(['Invalid data length in "' file '"']);
    %    Greg concluded that this error is not always fatal - so program may proceed
    fprintf('\nError: Invalid data length in %s file. \nResults may be invalid! \n',file);
    fprintf('MIP.dataFooter = %d\nMIP.hdrLen = %d\nFile size (fsize) = %d\n',MIP.dataFooter,MIP.hdrLen,fsize);
end

% Calculate number of channels and points
MIP.nChannels=length(find(MIP.design.channels));
MIP.nPoints=(MIP.dataFooter-MIP.dataStart)./(MIP.nChannels+1)./2;

% If data_flag = 0, skip EEG data and read only the digital trigger channel
if(~data_flag)
    MIP.data=[];
    fseek(EEGfid,MIP.dataStart + MIP.nChannels*2,'bof');
    MIP.data=fread(EEGfid,MIP.nPoints,'1*int16',MIP.nChannels*2);
else
    % If data_flag = 1, return EEG
    MIP.data=[];
    fseek(EEGfid,MIP.dataStart,'bof');
    MIP.data=fread(EEGfid,[MIP.nChannels+1 MIP.nPoints],'int16');
end
fclose(EEGfid);
return
end

% Modification History:
%
% $Log: MIPRead.m,v $
% Revision 1.4  2005/02/21 23:32:46  michelich
% Changes by Charles Michelich & Francis Favorini:
% Remove garbage after null terminated strings.
% Cast uchars to chars when using sprintf to address MATLAB 7 warning.
%
% Revision 1.3  2005/02/03 16:58:19  michelich
% Changes based on M-lint:
% Make unused CVS variables comments for increased efficiency.
% Remove unnecessary semicolon after function declarations.
% Remove unnecessary commas after try, catch, and else statements.
%
% Revision 1.2  2003/10/22 15:54:34  gadde
% Make some CVS info accessible through variables
%
% Revision 1.1  2002/10/08 23:46:45  michelich
% Initial CVS Import
%
%
% Pre CVS History Entries:
% Francis Favorini, 06/15/01.

%***********************
% Function bandpass_eeg
%***********************

function y = bandpass_EEG(data, lowfreq, highfreq, sample_rate)

% Bandpass EEG data with a Butterworth filter
% newdata = bandpass_EEG(data, lowfreq, highfreq, sample_rate);
% if lowfreq = 0 -> lowpass filter only
% if highfreq = 0 -> highpass filter only

% Written by Martin McKeown

nyquist = 0.5*sample_rate;
filterorder = 5;

if highfreq == 0
    % high pass filter
    [b,a] = butter(filterorder, lowfreq/nyquist, 'high');
elseif lowfreq == 0
    % low pass filter
    [b,a] = butter(filterorder, highfreq/nyquist);
else
    [b,a] = butter(filterorder, [lowfreq highfreq]/nyquist);
end

y=zeros(size(data));
for i = 1:size(data,1)
    y(i,:) = filtfilt(b,a,data(i,:));
end
return
end

%***********************
% Function notch60
%***********************

function y = notch60(x,fs)
%
% function y = notch60(x,fs)
%
%  where x is a vector or matrix with cols >> rows, and fs is sample
%  frequency
%
% Uses a 5th order butterworth bandstop filter
%

if nargin < 2
    help notch60
    return
end

f=60;  % notch frequency
Wn = fs/2;
notch = f/Wn;
filter_order = 5;

[B,A] = butter(filter_order,[notch - 1/Wn, notch+ 1/Wn ],'stop');

y = zeros(size(x));

for i=1:size(x,1)
    y(i,:) = filtfilt(B,A,x(i,:));
end
return
end

%***********************
% Function getstrs
%***********************

function [strs,count]=getstrs(fid,n)
%GETSTRS Read lines from text file and put into cell array.
%
%       [strs,count]=GETSTRS(fid);
%       [strs,count]=GETSTRS(fid,n);
%
%       fid is the file identifier for an open file.
%       n optionally specifies the max number of lines to read.
%
%       strs is a cell array of strings.
%       count is the number of lines actually read.
%

% CVS ID and authorship of this code
% CVSId = '$Id: getstrs.m,v 1.4 2005/02/03 20:17:45 michelich Exp $';
% CVSRevision = '$Revision: 1.4 $';
% CVSDate = '$Date: 2005/02/03 20:17:45 $';
% CVSRCSFile = '$RCSfile: getstrs.m,v $';

if nargin<2, n=0; end
l=0;
strs={};
while 1
    s=fgetl(fid);
    if ~ischar(s), break; end                            % No more strings
    l=l+1;
    strs{l}=s;
    if l==n, break; end
end
strs=strs';
count=l;
return
end

% Modification History:
%
% $Log: getstrs.m,v $
% Revision 1.4  2005/02/03 20:17:45  michelich
% M-Lint: Replace deprecated isstr with ischar.
%
% Revision 1.3  2005/02/03 16:58:33  michelich
% Changes based on M-lint:
% Make unused CVS variables comments for increased efficiency.
% Remove unnecessary semicolon after function declarations.
% Remove unnecessary commas after try, catch, and else statements.
%
% Revision 1.2  2003/10/22 15:54:35  gadde
% Make some CVS info accessible through variables
%
% Revision 1.1  2002/08/27 22:24:15  michelich
% Initial CVS Import
%
%
% Pre CVS History Entries:
% Charles Michelich, 2001/01/23. Changed case of function name to all lowercase
% Francis Favorini,  1997/06/12. Returns a cell array of strings.
% Francis Favorini,  1997/01/14. n is optional.
% Francis Favorini,  1996/10/02.
