function [] = PostProcessing(S,model_info,f_casadi)
% --------------------------------------------------------------------------
% PostProcessing
%   This function calls subfunctions that post-process the simulation
%   results.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
% 
% Original author: Lars D'Hondt
% Original date: May/2022
%
% Last edit by: Gil Serrancolí
% Last edit date: September 2024
% --------------------------------------------------------------------------

% load results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'R');

%% Call post-processing subfunctions

% Saves a .mot file containing a full gait cycle, for visualisation in the
% OpenSim GUI.
[R] = PostProcess_write_motion_file(model_info,R);

%%
save(Outname,'R','-append');