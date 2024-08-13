function [rubberBand] = rubberBandBetweenAnkles(init, settings_orthosis)
% --------------------------------------------------------------------------
% rubberBandBetweenAnkles
%   Defines a rubber band connecting both ankles.
% 
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'rubberBandBetweenAnkles'  i.e. name of this function           
%       - stiffness:        stiffness in N/m
%       - zero_length:      length @ 0 stiffness force, in m
%   Values are set via S.orthosis.settings{i} in main.m
%
%
% OUTPUT:
%   - rubberBand -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 13/August/2024
% --------------------------------------------------------------------------

% ...

end