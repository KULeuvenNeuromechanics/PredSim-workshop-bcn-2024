%% filepaths
Misc.DataPath = fullfile(pathParamEst,['Data_' Misc.subjectName]);
Misc.model_path = fullfile(Misc.DataPath,'Model',['BCN_' Misc.subjectName '_PredSim_noSphere.osim']);
Misc.OutPath = fullfile(pathParamEst,['ResultsTest_' Misc.subjectName]); % folder to store results
addpath(genpath(fullfile(pathParamEst,'MuscleRedundancySolver')));

Misc.timePointsFile = ['timepoints_all_' Misc.subjectName '.xlsx'];

%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
trialTypes = {'s2s','sit2stand','stand2sit','gait_normal','gait_fast','squat','cmj','IPSA_KE_slow','IPSA_KE_fast','IPSA_DF_slow','IPSA_DF_fast','pendulum_sit_r_iso','pendulum_sit_l_iso','pendulum_sit_r_pm','pendulum_sit_l_pm'};
trialTypes_pendulum = {'pendulum_sit_r_iso','pendulum_sit_l_iso','pendulum_sit_r_pm','pendulum_sit_l_pm'};
trialTypes_functional = {'s2s','sit2stand','stand2sit','gait_normal','gait_fast','squat','cmj'};
trialInfo = readtable(fullfile(Misc.DataPath,Misc.timePointsFile));
time = [trialInfo.IC1 trialInfo.IC2];
Misc.trialName = trialInfo.Trial';
Misc.IKfile = strcat(pathParamEst,'\',trialInfo.IKfile)';
Misc.IDfile = strcat(pathParamEst,'\',trialInfo.IDfile)';
Misc.EMGfile = strcat(pathParamEst,'\',trialInfo.EMGfile)';
Misc.side = lower(trialInfo.side)';
Misc.trialType = trialInfo.movement';
Misc.trialFP = trialInfo.FP;
Misc.trialsIPSA = ismember(Misc.trialType,'IPSA');
Misc.ipsaDOF = trialInfo.ipsaDOF';
Misc.ipsaMovement = trialInfo.ipsaMovement';
Misc.ipsaSpeed = trialInfo.ipsaSpeed';
% Misc.ipsaRangeStart = trialInfo.ipsaRangeStart;
% Misc.ipsaRangeEnd = trialInfo.ipsaRangeEnd;
Misc.ipsaInfo = strcat(trialInfo.ipsaDOF,'_',trialInfo.ipsaMovement,'_',trialInfo.ipsaSpeed)';
Misc.trials_left = find(ismember(Misc.side,'l'));
Misc.trials_right = find(ismember(Misc.side,'r'));

% select the DOFs you want to include in the optimization
Misc.DOF_functional = {'hip_flexion_','knee_angle_','ankle_angle_'};
for t=1:length(Misc.trialType)
    if ismember(Misc.trialType{t},trialTypes_functional)
        Misc.DofNames_Input{t,1}=strcat(Misc.DOF_functional,Misc.side{t});
    elseif ismember(Misc.trialType{t},trialTypes_pendulum)
        Misc.DofNames_Input{t,1}={['knee_angle_',Misc.side{t}]};
    else
        Misc.DofNames_Input{t,1}={Misc.ipsaDOF{t}};
    end
end

% Name of the results file
Misc.OutName = strcat(Misc.trialName,{'_'},Misc.trialType,{'_'},Misc.side,{'_'})';

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%% Settings

% Set the tendon stifness of all muscles
Misc.kT = [];      % default way to set tendon stiffenss (default values is 35)

% Provide the correct headers in case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
% Misc.EMGheaders = {'Time','bifemlh_r','tib_ant_r','per_long_r','lat_gas_r',...
%     'bifemsh_r','soleus_r','vas_lat_r','vas_med_r','per_brev_l','tib_ant_l','per_long_l',...
%     'lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l',...
%     'rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l','glut_med2_r','rect_fem_r'};
% User can use either Misc.EMGheaders or Misc.EMGFileHeaderCorrespondence.
% Misc.EMGFileHeaderCorrespondence is a n x 2 cell array where first column
% is the header name in the EMG file and second column is the name of the
% muscle in the .osim model corresponding to that header name. The order of
% rows of this matrix does not have to match the order in the EMG file and
% not all rows in this matrix need to be in all the EMG files of all trials
Misc.EMGFileHeaderCorrespondence = {'Time', 'Time';...
    'time', 'Time';...
    'RREF', 'rect_fem_r';...
    'RVAL', 'vasti_r';...
    'RBIF', 'hamstrings_r';...
    'RTIA', 'tib_ant_r';...
    'RGAS', 'gastroc_r';...
    'RSOL', 'soleus_r';...
    'RGLU', 'glut_max_r';...
    'LREF', 'rect_fem_l';...
    'LVAL', 'vasti_l';...
    'LBIF', 'hamstrings_l';...
    'LTIA', 'tib_ant_l';...
    'LGAS', 'gastroc_l';...
    'LSOL', 'soleus_l';...
    'LGLU', 'glut_max_l'};

% channels you want to use for EMG constraints
Misc.EMGSelection = {'rect_fem_r','vasti_r','hamstrings_r','tib_ant_r','gastroc_r','soleus_r','glut_max_r',...
    'rect_fem_l','vasti_l','hamstrings_l','tib_ant_l','gastroc_l','soleus_l','glut_max_l'};

Misc.EMG_MuscleCopies = {};

% If the input EMG needs to be normalized
Misc.normalizeEMG = 1;
Misc.normalizeToMRS = 0;

% information for the EMG constraint
Misc.EMGconstr  = 1;     		% Boolean to select EMG constrained option

Misc = getMuscleProperties(Misc.model_path,Misc); % NOTE: getShift was being used before adjusting kT, it has been corrected here

% % parameter estimation
% Misc.Estimate_TendonStiffness = Misc.allMuscleList(1:Misc.nAllMuscList); % Names of muscles of which tendon stifness is estimated
Misc.Estimate_TendonStiffness = []; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.5; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 1.5; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_TendonStiffness = {}; % Couple muscles that should have equal tendons stiffness

Misc.Estimate_OptimalFiberLength= Misc.allMuscleList(1:Misc.nAllMuscList); % Names of muscles of which tendon stifness is estimated
Misc.lb_lMo_scaling= 0.4; % Lower bound for scaling generic tendon stiffness
Misc.ub_lMo_scaling = 2.5; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_fiber_length = {}; % Couple muscles that should have equal tendons stiffness

Misc.Estimate_TendonSlackLength= Misc.allMuscleList(1:Misc.nAllMuscList); % Names of muscles of which tendon stifness is estimated
Misc.lb_lTs_scaling= 0.4; % Lower bound for scaling generic tendon stiffness
Misc.ub_lTs_scaling = 2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_slack_length = {}; % Couple muscles that should have equal tendons stiffness

% select muscles
Misc.MuscleNames_Input = []; % select muscles % [] means all muscles included

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option

% % Set weights
Misc.wVm    = 0.002; % weight on slack term

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 0;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0; 	% TO DO: we should report results of EMG driven simulation as well

Misc.opt_sides = {'left' 'right'};
Misc.Mesh_Frequency = 25;

%% Run muscle tendon estimator
[Results,DatStore,Misc] = solveMuscleRedundancy_BCN_workshop(time,Misc);
save(fullfile(Misc.OutPath,[Misc.subjectName '_' Misc.AnalysisID '_Results.mat']),'Results','DatStore','Misc');

%% Prep for PredSim
% Add spheres
modelWithSphereAtCorrectLocation = fullfile(Misc.DataPath,'Model','withSphere',['BCN_' Misc.subjectName '_PredSim.osim']);
contact_spheres = get_contact_spheres(modelWithSphereAtCorrectLocation);
modelWithUpdated_lMo_lTs = fullfile(Misc.OutPath,Misc.newModelFile);
modelWithUpdated_lMo_lTs_PredSim = fullfile(Misc.OutPath,['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst.osim']);
add_contact_spheres(modelWithUpdated_lMo_lTs,contact_spheres,modelWithUpdated_lMo_lTs_PredSim);

% write tendon stiffness estimates in the form that PredSim can read
tendon_stiff_scale_PredSim = [Misc.allMuscleList',num2cell(Results.Param.kT_scaling_paramopt)]';
tendon_stiff_scale_PredSim = tendon_stiff_scale_PredSim(:)';
save(fullfile(Misc.OutPath,['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst_tendon_stiff_scale.mat']),'tendon_stiff_scale_PredSim');

%% Moving towards PredSim
cd(pathPredSim)
predSimSubjectName = ['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst'];
predSimSubjectDir =  fullfile(pathPredSim,'PredSim','Subjects',predSimSubjectName);
mkdir(predSimSubjectDir)
copyfile(fullfile(Misc.OutPath,['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst_tendon_stiff_scale.mat']),...
    fullfile(predSimSubjectDir,[predSimSubjectName '_tendon_stiff_scale.mat']));
copyfile(modelWithUpdated_lMo_lTs_PredSim,...
    fullfile(predSimSubjectDir,[predSimSubjectName '.osim']));
Misc.pathPredSim = pathPredSim;
Misc.predSimSubjectName = predSimSubjectName;

%% Run PredSim
%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clearvars -except Misc
close all
clc
% path to the repository folder
pathRepo = fullfile(Misc.pathPredSim,'PredSim');
% path to the folder that contains the repository folder
pathRepoFolder = Misc.pathPredSim;

%% Initialize S
pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)

[S] = initializeSettings('gait1018');
S.misc.main_path = pathRepo;

addpath(fullfile(S.misc.main_path,'VariousFunctions'))

%% Required inputs
% name of the subject
S.subject.name = Misc.predSimSubjectName;

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.subject.IG_selection = 'quasi-random';
if strcmp(Misc.subjectName,'CP1')
    S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','CP4_T0_10_IK_adjusted.mot');
    S.subject.IG_selection_gaitCyclePercent = 100;
elseif strcmp(Misc.subjectName,'CP2')
    S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','CP16_T0_11_IK_adjusted.mot');
    S.subject.IG_selection_gaitCyclePercent = 100;
end

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = Misc.useCluster;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
% S.bounds.a.lower            = ;
% S.bounds.SLL.upper          = ;
% S.bounds.SLR.upper          = ;
% S.bounds.dist_trav.lower    = ;
% S.bounds.t_final.upper      = ;
% S.bounds.t_final.lower      = ;
% S.bounds.Qs                 = {'pelvis_tilt',-30,30,'pelvis_list',-30,30};


% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = 100;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;
% S.misc.msk_geom_bounds      = {{'knee_angle_r'},0,90,{'mtp_angle_'},-50,20};
% S.misc.default_msk_geom_bound = ;
% S.misc.msk_geom_bounds      = {{'knee_angle_r','knee_angle_l'},-120,10,'lumbar_extension',nan,30};
S.misc.gaitmotion_type = 'FullGaitCycle';

% % S.post_process
S.post_process.make_plot = 0;
% S.post_process.savename  = 'datetime';
% S.post_process.load_prev_opti_vars = 1;
% S.post_process.rerun   = 1;
% S.post_process.result_filename = '';

% % S.solver
% S.solver.linear_solver  = '';
% S.solver.tol_ipopt      = ;
% S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
% S.solver.N_threads      = 10;
% S.solver.N_meshes       = 50;
if Misc.useCluster
    S.solver.par_cluster_name = Misc.clusterName;
end
S.solver.CasADi_path    = Misc.casadiPath;

% % S.subject
% S.subject.mass              = ;
% S.subject.IG_pelvis_y       = 1;
S.subject.adapt_IG_pelvis_y = 1;
if strcmp(Misc.subjectName,'CP1')
    S.subject.v_pelvis_x_trgt   = 1.1240;
elseif strcmp(Misc.subjectName,'CP2')
    S.subject.v_pelvis_x_trgt   = 1.0314;
end
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = {{'soleus','_gas','per_','tib_','_dig_','_hal_','FDB'},0.9}; %,'FDB'
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.tendon_stiff_scale = {{'soleus','gastroc'},0.5};
load(fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '_tendon_stiff_scale.mat']));
S.subject.tendon_stiff_scale = tendon_stiff_scale_PredSim;
% S.subject.scale_MT_params = {{'soleus_l'},'FMo',0.9,{'soleus_l'},'alphao',1.1};
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
% S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
% S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};
% S.subject.set_limit_torque_coefficients_selected_dofs = ...
%     {{'knee_angle_r','knee_angle_l'},-[11.03 -11.33 -6.09 33.94]',-[0.13 -2.4]',...
%     {'mtp_angle_r','mtp_angle_l'},-[0.18 -70.08 -0.9 14.87]',-[65/180*pi 0]'};
% S.subject.base_joints_legs = 'hip';
% S.subject.base_joints_arms = [];
% S.subject.mtp_type          = '2022paper';

% % S.weights
% S.weights.E         = 0;
% % S.weights.E_exp     = ;
% S.weights.q_dotdot  = 1;
% S.weights.e_arm     = 500;
% S.weights.pass_torq = 1;
% S.weights.a         = 10*18;
% S.weights.slack_ctrl = ;
% S.weights.pass_torq_includes_damping = ;

% %S.OpenSimADOptions: required inputs to convert .osim to .dll
S.OpenSimADOptions.compiler = Misc.visualStudioVersion;
S.OpenSimADOptions.verbose_mode = 0; % 0 for no outputs from cmake

        
% S.bounds.distanceConstraints = [];
% S.bounds.distanceConstraints(end+1).point1 = 'toes_l';
% S.bounds.distanceConstraints(end).point2 = 'ground';
% S.bounds.distanceConstraints(end).direction = 'y';
% S.bounds.distanceConstraints(end).lower_bound = 0.05;
% S.bounds.distanceConstraints(end).upper_bound = [];

%% Run predictive simulations

% warning wrt pelvis heigt for IG
if S.subject.adapt_IG_pelvis_y == 0 && S.subject.IG_selection ~= "quasi-random"
    uiwait(msgbox(["Pelvis height of the IG will not be changed.";"Set S.subject.adapt_IG_pelvis_y to 1 if you want to use the model's pelvis height."],"Warning","warn"));
end

% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    [savename] = run_pred_sim(S,osim_path);
end

%% Plot results
if S.post_process.make_plot && ~S.solver.run_as_batch_job
    % set path to saved result
    result_paths{2} = fullfile(S.subject.save_folder,[savename '.mat']);
    % add path to subfolder with plotting functions
    addpath(fullfile(S.misc.main_path,'PlotFigures'))
    % call plotting script
    run_this_file_to_plot_figures
end



