function [f_casadi] = createCasadi_GenHelper(S)
% --------------------------------------------------------------------------
% createCasadi_GenHelper
%   Function to create general Casadi functions.
%
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_casadi -
%   * Struct that contains all casadi functions.
%
% Original authors: Dhruv Gupta, Lars D'Hondt, Tom Buurke
% Original date: 01/12/2021
%
% Last edit by: Gil Serrancolí
% Last edit date: 29/08/2024
% --------------------------------------------------------------------------

import casadi.*

N_torq_act = 6;
N_pass_dof = 6;

%% Normalized sum of squared values
q=MX.sym('q',1,6);
expdata=MX.sym('expdata',1,6);
f_casadi.ftrackkinexp=Function('ftrackkinexp',{q,expdata},{sum((expdata([1:4 6])-q([1:4 6])).^2)});

force_exp=MX.sym('force_exp',1,2);
force_mod=MX.sym('force_mod',1,2);
f_casadi.ftrackforceexp=Function('ftrackforceexp',{force_mod,force_exp},{sum(((force_exp-force_mod)/1000).^2)});

uT=MX.sym('uT',1,6);
f_casadi.uTsqsum=Function('uTsqsum',{uT},{sum(uT.^2)});

ua=MX.sym('ua',1,6);
f_casadi.uasqsum=Function('uasqsum',{ua},{sum(ua.^2)});


