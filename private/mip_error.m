function mip_error(msg, varargin)
%     MIP_ERROR   Light wrapper around error
%         MIP_ERROR(MSG, VARARGIN)
% 
%     This function calls the 'error' function and then fclose('all').
%     
%     Created by Zarrar Shehzad on 2012-09-12.

if nargin == 0, msg = 'Program cannot proceed'; end

error(msg, varargin{:});
fclose('all');

end %  function