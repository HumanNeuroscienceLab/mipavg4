classdef logger
    properties
        fids
    end
    
    methods(Static)
        function mfprintf(fids, varargin)
            for i=1:length(fids), fprintf(fids(i), varargin{:}); end
        end
    end
    
    methods
        %% constructor
        function obj = logger
            obj.fids = [];
        end
        
        function obj = to_standard_output(obj)
            obj.fids(end+1) = 1;
        end
        
        function obj = to_standard_error(obj)
            obj.fids(end+1) = 2;
        end
        
        function obj = to_file(obj, filename)
            [fid msg] = fopen(filename, 'w');
            if fid == -1
                error('Failed to open log file "%s" for writing: %s', filename, msg);
            end
            obj.fids(end+1) = fid;
        end
        
        %% log message
        function write(obj, varargin)
            obj.mfprintf(obj.fids, varargin{:});
        end
        
        %% execute function, log its output, and return output of function
        function varargout = callfun(obj, fun, varargin)
            nargs             = eval(sprintf('nargout(@%s)', fun));
            varargout         = cell(1, nargs);
            cmd               = sprintf('%s(varargin{:})', fun);
            [T varargout{:}]  = evalc(cmd);
            T                 = sprintf('\t%s', regexprep(T, '\n', '\n\t'));
            obj.write(T(1:end-1));
        end
        
        %% error message
        function error(obj, msg, varargin)
            obj.write(['\nERROR: ' msg '\n'], varargin{:});
        end
        
        %% warn and log message
        function warn(obj, id, msg, varargin)
            fprintf('\n');
            warning(id, sprintf('%s\n', msg), varargin{:});
            fids = obj.fids(obj.fids ~= 1 & obj.fids ~= 2);
            obj.mfprintf(fids, ['\nWARNING: ' msg '\n'], varargin{:});
        end
        
        %% destructor
        function delete(obj)
            fids = obj.fids(obj.fids ~= 1 & obj.fids ~= 2);
            for i=1:length(fids), fclose(fids(i)); end
        end
    end
end