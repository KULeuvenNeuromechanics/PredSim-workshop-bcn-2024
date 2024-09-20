% --------------------------------------------------------------------------
% Modelling assistive devices
% hands-on sessions
% Example 1. Changing the parameters of an existing orthosis
%
%   This script requires PredSim and its dependencies to be installed.
%   See https://github.com/KULeuvenNeuromechanics/PredSim
%
% Original author: Lars D'Hondt
% Original date: 12/August/2024
% --------------------------------------------------------------------------

clear
close all
clc

%% Initialize S

% Full path to the folder containing the PredSim code
pathPredSim = 'C:\GBW_MyPrograms\PredSim-dev';

% Add the PredSim folder to the matlab search path
addpath(pathPredSim)

% Add the folder with default settings to the matlab search path. Other
% folders will be added automatically.
addpath(fullfile(pathPredSim,'DefaultSettings'))

% Initialize the settings. The argument 'gait1018' tells the function to 
% load the settings in PredSim/Subjects/gait1018/settings_gait1018.m
%   Note: This can be different from S.subject.name and osim_path.
[S] = initializeSettings('gait1018');

% Add the path to the PredSim repo to the settings
S.misc.main_path = pathPredSim;

%% Required inputs
% Minimally required inputs.

% Name of the subject. This is used to refer to a folder in
% PredSim/Subjects where intermediate files are stored.
S.subject.name = 'gait1018';

% Path to folder where simulation results should be saved
S.subject.save_folder  = fullfile(pwd,'Results',S.subject.name); 

% Provide an initial guess for the kinematics.
%   Option 1: warm-start
% Use a motion file with inverse kinematics of an average gait cycle of an
% example subject.
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
% Indicate that the motion file containt exactly 1 gait cycle.
%   Note: Motion files with a PredSim result contain 2 gait cycles (200%)
S.subject.IG_selection_gaitCyclePercent = 100;
%   Option 2: cold-start
% 'quasi-random' moves the model forward at the target velocity with all
% coordinates in their default position.
% S.subject.IG_selection = 'quasi-random';

% Path to the opensim model of the subject
%   Note: The name of the file and folder do not have to match S.subject.name
osim_path = fullfile(S.misc.main_path,'Subjects','gait1018','gait1018.osim');

% Run the simulation as a batch job (parallel computing toolbox). This can
% be useful when running many simulations.
S.solver.run_as_batch_job = false;

%% Optional inputs
% See PredSim/README.md for the full list of settings and their defaults.

% Path to CasADi libraries (top folder)
%   Note: If casadi is already in the matlab search path, you can use 
%   S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath()
S.solver.CasADi_path = 'C:\GBW_MyPrograms\casadi\v3.6.5';

% PredSim tries to automatically detect the compiler command that needs to 
% be passed to OpenSimAD. If this fails, manually set the command.
% S.OpenSimADOptions.compiler = 'Visual Studio 17 2022';

% Diplay cmake outputs (part of OpenSimAD workflow) in command window. This
% is only useful for debugging.
S.OpenSimADOptions.verbose_mode = false;


% stuff related to specific example ...



%% Run predictive simulations

% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    [savename] = run_pred_sim(S,osim_path);
end



