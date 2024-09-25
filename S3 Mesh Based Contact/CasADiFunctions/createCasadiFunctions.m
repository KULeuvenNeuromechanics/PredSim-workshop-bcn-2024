function [f_casadi] = createCasadiFunctions(S)
% --------------------------------------------------------------------------
%createCasadiFunctions.m
%   Overview function from which al casadi functions are created
% 
% INPUT:
%   - S -
%   * setting structure S
% OUTPUT:
%   - f_casadi -
%   * Struct containing all casadi functions.

%% Create generic casadi functions
f_casadi = createCasadi_GenHelper(S);




