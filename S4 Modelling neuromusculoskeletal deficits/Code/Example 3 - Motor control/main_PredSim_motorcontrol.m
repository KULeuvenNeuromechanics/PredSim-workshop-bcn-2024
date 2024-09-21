% --------------------------------------------------------------------------
% Modelling neuromusculoskeletal defictis
% hands-on sessions
% Example 3. Modelling motor control thorugh muscle synergies
%
%   This script requires PredSim and its dependencies to be installed.
%   See https://github.com/KULeuvenNeuromechanics/PredSim
%
% Original author: Míriam Febrer
% Original date: 05/September/2024
% --------------------------------------------------------------------------

clear
close all
clc

%% Choose the simulation that you want to run

% 1. Select the CP subject
% Subject CP1 ('BCN_CP1')
% Subject CP2 ('BCN_CP2')
subject_name = 'BCN_CP1';

% 2. Select the type of simulation
% No synergies ('NoSyn')
% Imposing the number of synergies ('SynN')
% Tracking synergy weights ('SynW')
motor_control = 'SynN';

%% Add folders to the matlab search path

% Full path to the folder containing the PredSim code -> update this path
% to yours
pathPredSim = 'C:\Users\febre\Documents\GitHub\PredSim_repo_results\PredSim';

% Add the PredSim folder to the matlab search path
addpath(pathPredSim)

% Add the folder with default settings to the matlab search path. Other
% folders will be added automatically.
addpath(fullfile(pathPredSim,'DefaultSettings'))

%% Initialize S

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
S.subject.name = subject_name;

% Path to folder where simulation results should be saved
S.misc.save_folder  = fullfile(pwd,'Results',S.subject.name);

% Since the model is not symmetric, we should do a full gait cycle
% simulation
S.misc.gaitmotion_type = 'FullGaitCycle';

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
% S.subject.IG_selection = 'quasi-random';

% Path to the opensim model of the subject
%   Note: The name of the file and folder do not have to match S.subject.name
osim_path = fullfile(pwd,'..\..\Models\',S.subject.name,[S.subject.name '_PredSim.osim']);

% Run the simulation as a batch job (parallel computing toolbox). This can
% be useful when running many simulations.
S.solver.run_as_batch_job = false;

%% Optional inputs
% See PredSim/README.md for the full list of settings and their defaults.

% Path to CasADi libraries (top folder)
%   Note: If casadi is already in the matlab search path, you can use
S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath();
% S.solver.CasADi_path = 'C:\GBW_MyPrograms\casadi\v3.6.5';

% PredSim tries to automatically detect the compiler command that needs to
% be passed to OpenSimAD. If this fails, manually set the command.
% S.OpenSimADOptions.compiler = 'Visual Studio 17 2022';

% Diplay cmake outputs (part of OpenSimAD workflow) in command window. This
% is only useful for debugging.
S.OpenSimADOptions.verbose_mode = false;


%% Synergy-related inputs
% Load the synergy analysis results.
load([S.subject.name,'_Syn.mat'])
switch motor_control
    case 'NoSyn'
        S.subject.synergies = 0; % 1 = implement muscle synergies
        S.subject.TrackSynW = 0; % 1 = track synergy weights
    case 'SynN'
        S.subject.synergies = 1; % 1 = implement muscle synergies
        S.subject.TrackSynW = 0; % 1 = track synergy weights
        S.subject.NSyn_r = SynN.R;
        S.subject.NSyn_l = SynN.L;
    case 'SynW'
        S.subject.synergies = 1; % 1 = implement muscle synergies
        S.subject.TrackSynW = 1; % 1 = track synergy weights
        S.subject.NSyn_r = SynN.R;
        S.subject.NSyn_l = SynN.L;
        S.subject.TrackSynW_side = 'RightLeft';
        % RIGHT
        % number of tracked synergies (may be different from the number of synergies)
        S.subject.TrackSynW_NSyn_r = SynN.R; 
        % Synergy weights to be tracked. These muscles correspond to
        % the muscles that have been measured experimentally
        S.subject.knownSynW_r = {'rect_fem_r', SynW.R(:,1)',...
            'vasti_r', SynW.R(:,2)',...
            'bifemsh_r', SynW.R(:,3)',...
            'hamstrings_r', SynW.R(:,4)',...
            'tib_ant_r', SynW.R(:,5)',...
            'gastroc_r', SynW.R(:,6)',...
            'soleus_r', SynW.R(:,7)',...
            'glut_max_r', SynW.R(:,8)'};
        % LEFT
        S.subject.TrackSynW_NSyn_l = SynN.L;
        S.subject.knownSynW_l = {'rect_fem_l', SynW.L(:,1)',...
            'vasti_l', SynW.L(:,2)',...
            'bifemsh_l', SynW.L(:,3)',...
            'hamstrings_l', SynW.L(:,4)',...
            'tib_ant_l', SynW.L(:,5)',...
            'gastroc_l', SynW.L(:,6)',...
            'soleus_l', SynW.L(:,7)',...
            'glut_max_l', SynW.L(:,8)'};
end



%% Run predictive simulations

% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    [savename] = run_pred_sim(S,osim_path);
end




