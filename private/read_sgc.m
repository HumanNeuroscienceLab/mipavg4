function [epoch baseline first_pt_msec ALLcodes ALLlabels] = read_sgc(SGCfiles, logg)
%   read_sgc - Read in epoch options from segmentation control file
%   [epoch baseline first_pt_msec ALLcodes ALLlabels] = read_sgc(SGCfiles, logg)
%   
%   Inputs:
%     SGCfiles      - a cell array of segmentation control files
%     logg          - instance of logger object
%   
%   Outputs:
%     epoch         - [onset offset] in secs for each trial
%     baseline      - [onset offset] in secs of baseline window
%     first_pt_msec - 1st point in trial in msecs
%     ALLcodes      - cell array with each element consisting of 2 columns: 
%                     [input-codes output-codes]
%     ALLlabels     - cell array with each element consisting of 2 columns: 
%                     [input-labels output-labels]
%   
%   Created by Zarrar Shehzad on 2012-09-06.


%% Short Functions
% returns the number of unique values across a cell array
number_of_unique = @(cell_array) length(unique(cat(1,cell_array{:})));


%%==============================================================================
%%                                                                    Initialize
%%==============================================================================

sgc_version   = 'sgc_v3.0';
nSGCfiles     = length(SGCfiles);
                                        % For each element:
ALLepochs     = cell(nSGCfiles, 1);     % onset/offset in ms of epoch
ALLbaselines  = cell(nSGCfiles, 1);     % onset/offset in ms of baseline window
ALLcodes      = cell(nSGCfiles, 1);     % 2 columns with code & recode
ALLlabels     = cell(nSGCfiles, 1);     % 2 columns with label & relabel

logg.write('\n%d SGC file(s)\n', nSGCfiles);


%%==============================================================================
%%                                                      Read data from SGC files
%%==============================================================================

for j = 1:nSGCfiles
    logg.write('File #%d: %s\n', j, SGCfiles{j});
  
    if ~exist(SGCfiles{j}, 'file')
        logg.error('SGC file does not exist --> %s\n', SGCfiles{j});
        error('Program cannot proceed');
    else
        fid = fopen(SGCfiles{j});
        
        % 1st line: version
        version = strtrim(fscanf(fid,'%s',1));
        if(~strcmpi(version, sgc_version))
            logg.error('SGC file is not current version\n');
            error('Program cannot proceed');
        end
        % 2nd line: onset/offset in ms of epoch
        epoch = fscanf(fid,'%d%d',2);
        % 3rd line: onset/offset in ms of baseline window
        baseline  = fscanf(fid,'%d%d',2);
        % 4th line: 1st epoch point in ms for EEGAD
        first_pt_msec = fscanf(fid,'%d',1);
        % 5+ lines: input codes/labels and output codes/labels
        in_codes = []; in_labels = []; out_codes = []; out_labels = [];
        while ~feof(fid)
            codes_and_labels = textscan(fid, '%d %s %d %s\n');
            in_codes = [in_codes; codes_and_labels{1}];
            in_labels = [in_labels; codes_and_labels{2}];
            out_codes = [out_codes; codes_and_labels{3}];
            out_labels = [out_labels; codes_and_labels{4}];
        end
        
        fclose(fid);
        
        % check that number of unique codes and labels are the same
        if length(unique(out_codes)) ~= length(unique(out_labels))
            logg.error(['The number of unique output codes and labels do not match.\n ' ...
                        'Check the the 3rd and 4th columns of your sgc file\n']);
            error('Program cannot proceed');
        end
        
        ALLepochs{j} = epoch; ALLbaselines{j} = baseline;
        ALLcodes{j} = [in_codes, out_codes]; ALLlabels{j} = [in_labels, out_labels];
    end
end

clearvars epoch baseline codes_and_labels in_codes in_labels out_codes out_labels;


%%==============================================================================
%%                                      Check epoch/baseline limits and keep one
%%==============================================================================

% Make sure that the epoch and baseline limits are the same in all SGC files
if number_of_unique(ALLepochs) ~= 2    % function below
    logg.error('An inconsistency was detected in the epoch limits within the SGC files\n');
    error('Program cannot proceed');
end
if number_of_unique(ALLbaselines) ~= 2 % function below
    logg.error('An inconsistency was detected in the baseline limits within the SGC files\n');
    error('Program cannot proceed');
end

% Only keep one set of epoch/baseline info
% and save in seconds
epoch = ALLepochs{1}./1000;
baseline = ALLbaselines{1}./1000;
clear ALLepochs ALLbaselines;

% Check epoch/baseline timing
if(epoch(1) >= epoch(2) || baseline(1) >= baseline(2) || baseline(1) < epoch(1) || baseline(2) > epoch(2))
    logg.error('Inconsistency in designated epochs\n');
    error('Program cannot continue');
end

logg.write('\n');
logg.write('Epoch begin and end (msec): %d %d\n',epoch(1),epoch(2));
logg.write('Baseline begin and end (msec) : %d %d\n',baseline(1),baseline(2));

end % function
