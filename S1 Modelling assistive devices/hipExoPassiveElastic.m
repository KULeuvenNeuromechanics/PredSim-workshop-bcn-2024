function [exo] = hipExoPassiveElastic(init, settings_orthosis)
% --------------------------------------------------------------------------
% hipExoPassiveElastic
%   Defines a passive hip exoskeleton as a rotational spring 
% 
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'hipExoPassiveElastic'  i.e. name of this function           
%       - stiffness:        stiffness in Nm/rad
%       - left_right:       'l' for left or 'r' for right
%   Values are set via S.orthosis.settings{i} in main.m
%
%
% OUTPUT:
%   - exo -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 13/August/2024
% --------------------------------------------------------------------------

% create Orthosis object
exo = Orthosis('exo',init);

% read settings that were passed from main.m
hip_stiffness = settings_orthosis.stiffness; % stiffness in Nm/rad
side = settings_orthosis.left_right; % 'l' for left or 'r' for right

% ...

end