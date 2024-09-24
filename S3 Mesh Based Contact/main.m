%% Tracking Simulations of Knee Pose Using a Mesh-Based Contact model

% This script starts the tracking simulation. The required inputs are 
% necessary to start the simulations. Optional inputs, if left empty, will 
% be taken from getDefaultSettings.m.

% This assumes that you have the libraries of CasADi and the MATLAB path
% points to them. We also assume that you work on Windows 10 or newer.

% Authors: Gil Serrancolí, Mohanad Harba and Joan Badia
% Last edit: 24/09/2024

clear
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For the workshop, we will mainly change those parameters:
geom_mesh_type='100x188'; %either 49x188, 100x188
radiusAccFaces='1'; % '05' is for 0.5 cm, or 1 is for 1 cm 
N_meshes=10;
preoptimizekneedofs=2; %0 is considering all dofs as experimental 
% data or 0 values, 1 is preoptimizing knee ty, 2 is preoptimizing knee
% adduction and knee ty
tol_ipopt = 3;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));

% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% Initialize S
pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)

[S] = initializeSettings();

S.post_process.result_filename='PoseEstimation';

%% Required inputs
% name of the subject
S.subject.name='KneeProsthesis';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepo,'TrackSimResults',S.subject.name); 

% subject folder to save intermediate data
S.misc.subject_path = fullfile(S.misc.main_path,'Subjects',S.subject.name);

% number of mesh intervals
S.solver.N_meshes=N_meshes;

% Initial and final time to track data
S.initial_time=2.2292;
S.final_time  =3.3740;

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepo,'TrackSimResults',S.subject.name); 

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

% Either use 'quasi-random' or 'data-informed' initial guess
S.subject.IG_selection='data-informed';

%% Get Experimental Data
pathExperimentalData = fullfile(pathRepo,'ExperimentalData');
addpath(pathExperimentalData);

[S] = getExperimentalData(S);

%Define cost function weights
S.weights.q_dotdot = 0.1; %0.1
S.weights.trackq = 100; %100
S.weights.trackKCF = 10; %10
S.weights.minresidualforce = 10; %10
S.weights.minsecq = 1; %1

% smoothing parameter into the contact model (preventing high negative
% pressures)
S.options.kpress='k5e5'; 
% Mesh
S.options.geom_mesh=geom_mesh_type; %either 49x188, 100x188
S.options.radiusAccFaces=radiusAccFaces; % '05' is for 0.5 cm, or 1 is for 1 cm 
if strcmp(S.options.geom_mesh,'49x188')&&strcmp(S.options.radiusAccFaces,'05')
    name_npairs='pairs1407';
elseif strcmp(S.options.geom_mesh,'49x188')&&strcmp(S.options.radiusAccFaces,'1')
    name_npairs='pairs2139';
elseif strcmp(S.options.geom_mesh,'100x188')&&strcmp(S.options.radiusAccFaces,'05')
    name_npairs='pairs2882';
elseif strcmp(S.options.geom_mesh,'100x188')&&strcmp(S.options.radiusAccFaces,'1')
    name_npairs='pairs4191';
end

%wheter to scale knee ty coordinate from the expected mean and appropriate 
% range value (not containing 0 value)
S.options.scale_kneety_av=1;
% minimize force tracking error, or impose it as constraint
S.options.forcetrackingAs='costfun'; %'costfun' or 'constraint'
% use knee residual forces (med and lat)
S.options.useResidualforces=0;
S.options.preoptimizekneedofs=preoptimizekneedofs; %0 is considering all dofs as experimental 
% data or 0 values, 1 is preoptimizing knee ty, 2 is preoptimizing knee
% adduction and knee ty

% name of the .dll containing the contact model
% Inputs: 6 coordinates (3 translations and 3 rotations of the knee)
% Outputs: 2 forces (medial and lateral compartments)
S.misc.external_function=['f_bothcompartments_' S.options.kpress '_' S.options.geom_mesh '_at'  S.options.radiusAccFaces 'cm_' name_npairs  '_OpenSim.dll'];

% type of parallel computing
S.solver.parallel_mode = 'thread';
% number of threads in parallel mode
S.solver.N_threads = 4;

% maximal amount of itereations after wich the solver will stop
S.solver.max_iter = 10000;
S.solver.linear_solver='mumps';
% the power (10^-x) the error has to reach before the OCP can 
% be regarded as solved; a higher number gives a more precise answer, but
% requires more time
S.solver.tol_ipopt = tol_ipopt;
S.solver.warm_start_init_point='yes';

% Load specific configuration for ipopt
S.post_process.load_prev_opti_vars=0;

% Start simulation
if S.solver.run_as_batch_job
    add_track_sim_to_batch(S,osim_path)
else
    [savename] = run_track_sim(S,osim_path);
end

%% Plot results
if ~S.solver.run_as_batch_job
    
    % call plotting script
    plot_figures(S);
end
