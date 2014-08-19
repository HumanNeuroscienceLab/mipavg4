function isWritable = can_write(varargin)

%CAN_WRITE Tests if a folder is writable.
%   CAN_WRITE(FOLDERNAME) returns 1 if the folder indicated by FOLDERNAME is
%   writable and 0 otherwise.
%
%   CAN_WRITE tests whether the directory, FOLDERNAME, is writable by
%   attempting to create a subdirectory within FOLDERNAME.  If the
%   subdirectory can be created then FOLDERNAME is writable and the
%   subdirectory is deleted.  CAN_WRITE differs from FILEATTRIB in that
%   CAN_WRITE captures network and share level permissions where FILEATTRIB
%   does not.
%
%   Example:
%
%   x = can_write('C:\WINDOWS');      returns whether or not the Windows
%                                     directory is writable
%
%   x = can_write('\\VT1\public');    returns whether or not the network  
%                                     share, 'public' is writable
%
%   x = can_write                     returns whether or not the present 
%                                     working directory is writable
%
%
%   See also FILEATTRIB

%   Written by Chris J Cannell
%   Contact ccannell@mindspring.com for questions or comments.
%   12/08/2005

if nargin > 1
    error('Too many input arguments');
elseif nargin == 1
    folderName = varargin{1};
else
    folderName = pwd;
end

% create a random folder name so no existing folders are affected
testDir = ['deleteMe_', num2str(floor(rand*1e12))];

% check if folderName exists
if exist(folderName, 'dir') == 7
    % test for write permissions by creating a folder and then deleting it
    [isWritable,message,messageid] = mkdir(folderName, testDir);
    % check if directory creation was successful
    if isWritable == 1
        % we have permission to write so delete the created test directory
        [status,message,messageid] = rmdir(fullfile(folderName, testDir));
        if status ~= 1
            disp(['Warning: Test folder "',...
                fullfile(folderName, testDir), '" could not be deleted']);
        end
    end
else
    error(['Folder "', folderName, '" does not exist!']);
end