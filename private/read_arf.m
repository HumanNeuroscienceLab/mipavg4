function ar_settings = read_arf(ARFfile, data_channels, standard_epoch, logg)
%   read_arf - Reads in artifact rejection settings from *.arf file
%   [ar_setting] = read_arf(ARFfile, data_channels, standard_epoch, logg)
% 
%   Inputs:
%       ARFfile         - filename with artifact rejection settings 
%       data_channels   - cell array with names of all the channels
%       standard_epoch  - [onset offset] in secs, default epoch
%       logg            - instance of logger object
%   
%   Outputs:
%       ar_settings     - a struct array with fields: method, channels, 
%                         criteria, prestim, poststim, and opts
%   
%   Created by Zarrar Shehzad on 2012-09-09.
    
    %% Short Functions
    str_is_number = @(x) all(isstrprop(x, 'digit'));
    
    
    %% Check/set input args
    
    if ~exist(ARFfile, 'file')
        logg.error('File "%s" does not exist\n', ARFfile);
        error('Program cannot proceed');
    end
    
    logg.write('\nReading in ARF file: %s\n', ARFfile);
    
    fid = fopen(ARFfile);
    
    valid_methods = {'flat', 'ppa', 'zthr'};
    ar_settings = struct('method', [], 'channels', [], 'criteria', [], ...
                         'prestim', [], 'poststim', [], 'opts', []);
    
    logg.write(['\nArtifact rejection option set for ...\n' ...
                '\nchan\tfunc\tcrit\tstart\tend \n']);
    
    
    %% Read in ARF file
    i = 0;
    while ~feof(fid)
        i = i + 1;
        arf_string = strtrim(fgetl(fid));
        if isempty(arf_string), continue; end
        [channel method criteria epoch_start epoch_end opt1 opt2] = ... 
            strread_with_textscan(arf_string, '%s %s %d %d %d %s %s');
        
        % Channel
        if isempty(channel); channel = 'all'; end
        if str_is_number(channel); channel = str2num(channel); end
        channels = ft_channelselection(channel, data_channels);
        if isnumeric(channel); channel = num2str(channel); end
        if isempty(channels)
            logg.error('Line #%d: Selected channels (%s) not found', i, channel);
            error('Program cannot proceed');
        end
        
        % Method
        method = lower(method);
        if ~any(ismember(valid_methods, method))
            logg.error('Line #%d: Unrecognized method "%s"\n', i, method);
            error('Program cannot proceed');
        end
        
        % Epoch
        if isempty(epoch_start); epoch_start = standard_epoch(1); 
        else epoch_start = double(epoch_start)/1000; end
        if isempty(epoch_end); epoch_end = standard_epoch(2); 
        else epoch_end = double(epoch_end)/1000; end
        
        if epoch_start < standard_epoch(1) || epoch_end > standard_epoch(2)
            logg.error('ARF temporal range (%.3f to %.3f) exceeds raw epoch limits (%.3f to %.3f)\n', ...
                        epoch_start, epoch_end, standard_epoch(1), standard_epoch(2));
            error('Program cannot proceed');
        end
        
        % Additional options
        opts = {};
        if ~isempty(opt1); opts{end+1} = opt1; end
        if ~isempty(opt2); opts{end+1} = opt2; end
        
        % Save
        ar_settings(i).channels = channels;
        ar_settings(i).method   = method;
        ar_settings(i).criteria = double(criteria);
        ar_settings(i).prestim  = -1 * epoch_start;
        ar_settings(i).poststim = epoch_end;
        ar_settings(i).opts     = opts;
        
        % Log
        logg.write('%s\t%s\t%d\t%.3f\t%.3f\n', channel, method, criteria,  ...
                                               epoch_start, epoch_end);
    end 
    logg.write('\n');
    fclose(fid);
    
    
    %% Function
    function [varargout] = strread_with_textscan(string, format, varargin)
    %   strread_with_textscan - Implements deprecated strread using textscan
    %   [varargout] = strread_with_textscan(string, format, varargin)
    % 
    %   Inputs:
    %       string      - Input string to be parsed
    %       format      - Format (like that for sprintf) to parse string
    %       varargin    - See help for textscan
    %
    %   Outputs:
    %       varargout   - Depending on the number of elements indicated in format
    %                     will return that number of elements. 
    %
    %   Usage:
    %       >> [a,b,c,...] = strread('string', 'format')
    %   
    %   Example:
    %       >> s = sprintf('a,1,2\nb,3,4\n');
    %       >> [a,b,c] = strread(s,'%s%d%d','delimiter',',')
    %
    %   Created by Zarrar Shehzad on 2012-09-09.
        varargout = textscan(string, format, varargin{:});
        for ii=1:length(varargout); 
            if ~isempty(varargout{ii}) && iscell(varargout{ii})
                varargout{ii} = varargout{ii}{1}; 
            end
        end
    end
end
