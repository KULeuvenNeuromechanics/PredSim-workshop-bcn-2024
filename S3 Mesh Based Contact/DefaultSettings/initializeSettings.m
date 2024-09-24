function [S] = initializeSettings(varargin)
% --------------------------------------------------------------------------
%initializeSettings
%   This function creates the empty settings struct S up to the field 
%   above the field containing data. 
%
% OUTPUT:
%   - S -
%
% --------------------------------------------------------------------------

S = struct;

S.misc         = [];
S.post_process = [];
S.solver       = [];
S.weights      = [];

% save computername
S.misc.computername = getenv('COMPUTERNAME');

% save path to repo
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
S.misc.main_path = pathRepo;

S.post_process.rerun=0;

if ~isempty(varargin) %not used here

    [pathDefaultSettings,~,~] = fileparts(mfilename('fullpath'));
    [pathRepo,~,~] = fileparts(pathDefaultSettings);

    reference_path = fullfile(pathRepo,'Subjects',varargin{1},['settings_',varargin{1},'.m']);

    if isfile(reference_path)
        disp(['Initialising settings from "',reference_path,'".'])
        run(reference_path);
    else
        warning(['Could not initialise from "',reference_path,'". ',...
            'Ignoring input argument "',varargin{1},'".']);
    end


end







