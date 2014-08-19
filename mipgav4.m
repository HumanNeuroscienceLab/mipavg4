function mipgav4(varargin)
%     MIPGAV4   Computes the grand average ERP across subjects
%         mipgav4(varargin)
% 
%     TODO
%     
%     Created by Zarrar Shehzad on 2012-10-02.


typeAVG = 0; typeSET = 0; typeMAT = 0;
input_files = {}; input_hdr = ''; log_file = '';
nInputs = 0; nLOGfiles = 0; nOUTprefixes = 0;
diff_conditions = {};
force = 0; verbose = 1; to_logfile = 1;
output = []; output_opts = {'eegad', 'econnectome'};
for i=1:length(output_opts); output.(output_opts{i}) = 0; end


%%==============================================================================
%%                                                          Parse user arguments
%%==============================================================================

for j = 1:length(varargin)
    arg_string = lower(char(varargin{j}));
    [~, ~, arg_ext] = fileparts(arg_string);
    % Files
    if strcmp(arg_ext, '.avg')
        input_files{end+1} = varargin{j};
        typeAVG = 1;
        nInputs = nInputs + 1;
    elseif strcmp(arg_ext, '.set')
        input_files{end+1} = varargin{j};
        typeSET = 1;
        nInputs = nInputs + 1;
    elseif strcmp(arg_ext, '.mat')
        input_files{end+1} = varargin{j};
        typeMAT = 1;
        nInputs = nInputs + 1;
    elseif strcmp(arg_ext, '.hdr')
        input_hdr = varargin{j};
    elseif strcmp(arg_ext, '.log')
        log_file = varargin{j};
        nLOGfiles = nLOGfiles + 1;
    % Options
    elseif strfind(arg_string, '-d')
        arg_string = char(varargin{j});
        opt = strtrim(arg_string(strfind(lower(arg_string),'-d')+2:end));
        % Get the different conditions
        diff_conditions = regexp(opt, ';\ *', 'split');
        % Parse 'Label1 - Label2' into {'Label1', 'Label2'}
        diff_conditions = regexp(diff_conditions, '(\w+)\ *-\ *(\w+)', ...
                                    'tokens');
        diff_conditions = cellfun(@(x) x{1}, diff_conditions, ...
                                    'UniformOutput', false);        
    elseif strfind(arg_string, '-f')
        force = 1;
    % Output Options
    elseif strfind(arg_string, '-p')
        out_prefix = strtrim(varargin{j}(strfind(arg_string, '-p')+2:end));
        nOUTprefixes = nOUTprefixes + 1;
    elseif strfind(arg_string, '-o')
        o = lower(strtrim(varargin{j}(strfind(arg_string,'-o')+2:end)));
        if isempty(intersect(o, output_opts))
            error('Unrecognized argument ''%s'' for -o', o); 
        end
        output.(o) = 1;
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
%%     Check that we have all the necessary things and set defaults where needed
%%==============================================================================

if sum([typeAVG typeSET typeMAT]) > 1
    error('Cannot use any combination of .avg, .set, or .mat files');
end

if isempty(input_hdr), input_hdr = [remove_extension(input_files{1}) '.hdr']; end

if ~nOUTprefixes
    fprintf('\nNo output prefix was specified. Using MipGav.\n');
    out_prefix = 'MipGav';
elseif nOUTprefixes > 1
    error('More than one output prefix was specified. Program cannot proceed');
end

if ~any(struct2array(output))
    error('No output was specified. Program cannot proceed');
end

if ~nLOGfiles && to_logfile
    fprintf('\nNo LOG file was specified. Using %s.log.\n', out_prefix);
    log_file = [out_prefix '.log'];
elseif nLOGfiles > 1
    error('More than one LOG file specified. Program cannot proceed');
end

is_writable = can_write(dirname(out_prefix));
if ~is_writable
    error('Cannot write to the output directory %s', dirname(OUTprefix));
end
clear is_writable;

out_econn_prefix    = [out_prefix '_econn'];
diff_econn_prefix   = [out_prefix '_diff_econn'];
out_gav             = [out_prefix '_eegad.gav'];
out_hdr             = [out_prefix '_eegad.hdr'];
diff_gav            = [out_prefix '_diff_eegad.gav'];
diff_hdr            = [out_prefix '_diff_eegad.hdr'];


%%==============================================================================
%%                                                            Initialize logging
%%==============================================================================

logg = logger;

% Log to standard output if verbose
if verbose, logg = logg.to_standard_output; end

% Log to file
if to_logfile, logg = logg.to_file(log_file); end

% First messages!
logg.write('\nMIPGAV4 Program executed on %s\n', datestr(now));
opts = sprintf('''%s'', ', varargin{:});
logg.write('\nCalled:\nmipgav4(%s)\n\n', opts(1:end-2));


%%==============================================================================
%%                                                      Check inputs and outputs
%%==============================================================================

logg.write('\nChecking inputs and outputs\n');

% Check that all inputs exist
for input_file = [input_files input_hdr]
    check_input_file(input_file{1}, logg);
end

% Check if output file(s) exists
if output.eegad
    if exist(out_gav, 'file') && ~force
        error('Output %s already exists. Use -f if you want to overwrite.', out_gav);
    end
    if exist(diff_gav, 'file') && ~force
        error('Output %s already exists. Use -f if you want to overwrite.', diff_gav)
    end
end


%%==============================================================================
%%                                                                 Grand Average
%%==============================================================================

logg.write('\nReading in data and averaging');

% Read in data and average
hdr = EEGRead2(input_files{1});
avg = hdr.data;
for i = 2:length(input_files)
    tmp = EEGRead2(input_files{i});
    avg = avg + tmp.data;
end
avg = avg ./ length(input_files);

logg.write('\nSaving');

if output.eegad
    % Save
    logg.write('\nSaving average to %s', out_gav);
    fidAVG = fopen(out_gav, 'wb');
    for i = 1:hdr.nBins
        for j = 1:hdr.nChannels
            fwrite(fidAVG, avg(i,j,:), 'float32');
        end
    end
    fclose(fidAVG);

    % Copy over header
    logg.write('\nCopying over header file %s => %s\n', input_hdr, out_hdr);
    copyfile(input_hdr, out_hdr);
end

if output.econnectome
    for i = 1:hdr.nBins
        condition = hdr.binNames{i};
        
        % Convert
        EEG = eegad2econnectome(hdr, squeeze(avg(i,:,:))');
        
        % Save
        out_econn = [out_econn_prefix '_' condition '.mat'];
        if exist(out_econn, 'file') && ~force
            error('Output %s already exists. Use -f if you want to overwrite.', out_econn);
        end
        save(out_econn, 'EEG');
    end    
end


%%==============================================================================
%%                                                  Difference of Grand Averages
%%==============================================================================

if ~isempty(diff_conditions)
    logg.write('\nComputing difference between averages');
    
    % Setup
    diff = zeros(length(diff_conditions), hdr.nChannels, hdr.nPoints);

    % Compute differences
    for i = 1:length(diff_conditions)
        condition = diff_conditions{i};
    
        ind1 = find(strcmp(hdr.binNames, condition{1}));
        ind2 = find(strcmp(hdr.binNames, condition{2}));
    
        if isempty(ind1)
            error('Could not find difference bin name %s', condition{1});
        elseif isempty(ind2)
            error('Could not find difference bin name %s', condition{2});
        end
    
        diff(i,:,:) = avg(ind1,:,:) - avg(ind2,:,:);
    end

    if output.eegad
        % Save
        logg.write('\nSaving difference of averages to %s', diff_gav);
        fidDIFF = fopen(diff_gav, 'wb');
        for i = 1:length(diff_conditions)
            for j = 1:hdr.nChannels
                fwrite(fidDIFF, diff(i,j,:), 'float32');
            end
        end
        fclose(fidDIFF);
    
        % Write to header
        logg.write('\nWriting difference header file to %s\n', diff_hdr);
        fidHDR = fopen(diff_hdr, 'w');
        fprintf(fidHDR,'%s\n', hdr.expName); 
        fprintf(fidHDR,'%s\n', hdr.expDate);
        fprintf(fidHDR,'%d\n', hdr.nChannels);
        fprintf(fidHDR,'%d\n', hdr.nPoints);
        fprintf(fidHDR,'%d\n', hdr.sampling);
        fprintf(fidHDR,'%d\n', hdr.uvunits); % Microvolt conversion factor in raw data units/microvolt
        fprintf(fidHDR,'%d\n', hdr.onset);
        % Each difference on one line (read in as bin name)
        for k = 1:length(diff_conditions)
            label = sprintf('%s - %s', diff_conditions{k}{:});
            fprintf(fidHDR,'%s\n',label);
        end
        % Each channel name on one line
        for k = 1:hdr.nChannels, fprintf(fidHDR,'%s\n', hdr.chanNames{k}); end
    end
    
    if output.econnectome
        hdr.binNames = cellfun(@(x) sprintf('%s - %s', x{:}), diff_conditions, 'UniformOutput', false);
        hdr.nBins = length(hdr.binNames);
        
        for i = 1:hdr.nBins
            condition = hdr.binNames{i};
            
            % Convert
            EEG = eegad2econnectome(hdr, squeeze(diff(i,:,:))');
            
            % Save
            diff_econn = [diff_econn_prefix '_diff_' condition '.mat'];
            if exist(diff_econn, 'file') && ~force
                error('Output %s already exists. Use -f if you want to overwrite.', diff_econn);
            end
            save(diff_econn, 'EEG');
        end
        
    end
end

end %  function
