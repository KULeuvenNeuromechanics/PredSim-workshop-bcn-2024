% --------------------------------------------------------------------------
% Modelling assistive devices
% hands-on sessions
% Example 2. Defining a custom orthosis
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


%% Required inputs
% Minimally required inputs.

% Name of the subject. This is used to refer to a folder in
% PredSim/Subjects where intermediate files are stored.
S.subject.name = 'gait1018';

% Path to folder where simulation results should be saved
S.misc.save_folder  = fullfile(pwd,'Results','example_2',S.subject.name); 

% Provide an initial guess for the kinematics.
%   Option 1: warm-start
% Use a motion file with inverse kinematics of an average gait cycle of an
% example subject.
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
% Indicate that the motion file containt exactly 1 gait cycle.
%   Note: Motion files with a PredSim result contain 2 gait cycles (200%)
S.solver.IG_selection_gaitCyclePercent = 100;
%   Option 2: cold-start
% 'quasi-random' moves the model forward at the target velocity with all
% coordinates in their default position.
% S.solver.IG_selection = 'quasi-random';

% Path to the opensim model of the subject
%   Note: The name of the file and folder do not have to match S.subject.name
osim_path = fullfile(S.misc.main_path,'Subjects','gait1018','gait1018.osim');

% Run the simulation as a batch job (parallel computing toolbox). This can
% be useful when running many simulations.
S.solver.run_as_batch_job = false;

%% Optional inputs
% See PredSim/README.md for the full list of settings and their defaults.

% PredSim tries to automatically detect the compiler command that needs to 
% be passed to OpenSimAD. If this fails, manually set the command.
% S.OpenSimADOptions.compiler = 'Visual Studio 17 2022';

% Display cmake outputs (part of OpenSimAD workflow) in command window. This
% is only useful for debugging.
% S.OpenSimADOptions.verbose_mode = true;

% When running the simulation as a batch job, pass the path to the current
% folder to the job such that it can find the function hipExoPassiveElastic
% or rubberBandBetweenAnkles.
% S.solver.batch_job_paths = {pwd};

%% 2.1 Add a hip exoskeleton to the model

% Select 'hipExoPassiveElastic'. Before running a simulation, this function
% needs to be completed
S.orthosis.settings{1}.function_name = 'hipExoPassiveElastic';

% Add the exoskeleton to the right side
S.orthosis.settings{1}.left_right = 'r';

% Set the stiffness
S.orthosis.settings{1}.stiffness = 10; % Nm/rad

% Add the same exoskeleton on the left side
S.orthosis.settings{2}.function_name = 'hipExoPassiveElastic';
S.orthosis.settings{2}.left_right = 'l';
S.orthosis.settings{2}.stiffness = 10; % Nm/rad


%% 2.2 Add a rubber band connecting both ankles of the model

% % Set velocity
% S.misc.forward_velocity = 2.67;
% 
% % Remove previously added orthoses. 
% S.orthosis.settings = {};
% 
% % Select 'rubberBandBetweenAnkles'. Before running a simulation, this function
% % needs to be completed
% S.orthosis.settings{end+1}.function_name = 'rubberBandBetweenAnkles';
% 
% % Set the stiffness
% S.orthosis.settings{end}.stiffness = 120; % N/m
% 
% % Set the length at which the rubber band is not stretched
% S.orthosis.settings{end}.zero_length = 0.25; % m

%% Run predictive simulations

% Start simulation
runPredSim(S,osim_path);




