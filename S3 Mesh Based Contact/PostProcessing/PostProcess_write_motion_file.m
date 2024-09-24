function [R] = PostProcess_write_motion_file(model_info,R)
% --------------------------------------------------------------------------
% PostProcess_write_motion_file
%   This function creates a motionfile with 2 steps of predicted gait.
% 
% INPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
%   - R -
%   * struct with simulation results
%
% OUTPUT:
%   - R -
%   * struct with simulation results
% 
% Original author: Lars D'Hondt
% Original date: 12/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% Two gait cycles
t_mesh = R.time.coll;
% Joint angles
q_opt_GUI = R.kinematics.q_opt_unsc_all.rad;
q_opt_GUI(:,1:3)=q_opt_GUI(:,1:3)*180/pi;
q.data=[t_mesh' q_opt_GUI];
q.labels=['time' R.colheaders.coordinates];

write_motionFile(q,[R.S.subject.save_folder '\' R.S.subject.name '.mot']);