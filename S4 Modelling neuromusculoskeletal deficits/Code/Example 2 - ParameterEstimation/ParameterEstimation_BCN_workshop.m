%% Example EMG constrained simulation for ankle, knee and hip in the sagital plane

% In this example we estimate parameters of multiple lower-limb muscles
% using an EMG driven simulation of the ankle-knee and hip.
% Additionally, this example handles different trials with different number
% of degrees of freedom together.
% This file is based on https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver/blob/master/Examples/Example_EMGWalking/EMGconstrained_variableDOF.m

%% clear variables and command window
clear all; clc; close all;

%% Required Inputs

% 1. Subject name
Misc.subjectName = 'CP1';

% 2. specify path of Seminar folder
pathSeminar = 'C:\GBW_MyPrograms\PredSim-workshop-bcn-2024\S4 Modelling neuromusculoskeletal deficits';

% (automated) path of ParameterEstimation folder
pathParamEst = fullfile(pathSeminar,'\Code\Example 2 - ParameterEstimation');

% 3. specify path of MuscleRedundancySolver folder
pathMRS = 'C:\GBW_MyPrograms\MuscleRedundancySolver';

% 4. specify path to CasADi folder
Misc.casadiPath = 'C:\GBW_MyPrograms\casadi_3_6_6';

% 5. Name given to the estimation you run
Misc.AnalysisID = 'v1'; 

% (automated) path to the data folder
Misc.DataPath = fullfile(pathSeminar,'Data',['Data_' Misc.subjectName]);

%% Setup 
Misc.model_path = fullfile(Misc.DataPath,'Model',['BCN_' Misc.subjectName '_PredSim_noSphere.osim']);
Misc.OutPath = fullfile(pathParamEst,['Results_' Misc.subjectName]); % folder to store results
addpath(genpath(pathMRS));

Misc.timePointsFile = ['timepoints_all_' Misc.subjectName '.xlsx'];

trialTypes = {'s2s','sit2stand','stand2sit','gait_normal','gait_fast','squat','cmj','IPSA_KE_slow','IPSA_KE_fast','IPSA_DF_slow','IPSA_DF_fast','pendulum_sit_r_iso','pendulum_sit_l_iso','pendulum_sit_r_pm','pendulum_sit_l_pm'};
trialTypes_pendulum = {'pendulum_sit_r_iso','pendulum_sit_l_iso','pendulum_sit_r_pm','pendulum_sit_l_pm'};
trialTypes_functional = {'s2s','sit2stand','stand2sit','gait_normal','gait_fast','squat','cmj'};
trialInfo = readtable(fullfile(Misc.DataPath,Misc.timePointsFile));
time = [trialInfo.IC1 trialInfo.IC2];
Misc.trialName = trialInfo.Trial';
Misc.IKfile = strcat(Misc.DataPath,'\',trialInfo.IKfile)';
Misc.IDfile = strcat(Misc.DataPath,'\',trialInfo.IDfile)';
Misc.EMGfile = strcat(Misc.DataPath,'\',trialInfo.EMGfile)';
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

Misc.opt_sides = {'left' 'right'};



%% Optional Inputs
% upper and lower bound for difference between simulated and measured
% muscle activity (EMG)
Misc.EMGbounds  = [-0.1 0.1]; 

% upper and lower bound for EMG scaling
Misc.BoundsScaleEMG = [0.2 5];

Misc.Estimate_TendonStiffness = Misc.allMuscleList(1:Misc.nAllMuscList); % Names of muscles of which tendon stifness is estimated
% Misc.Estimate_TendonStiffness = []; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.5; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 1.5; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_TendonStiffness = {}; % Couple muscles that should have equal tendons stiffness

Misc.Estimate_OptimalFiberLength= Misc.allMuscleList(1:Misc.nAllMuscList); % Names of muscles of which tendon stifness is estimated
Misc.lb_lMo_scaling= 0.4; % Lower bound for scaling generic tendon stiffness
Misc.ub_lMo_scaling = 2.5; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_fiber_length = {}; % Couple muscles that should have equal tendons stiffness

Misc.Estimate_TendonSlackLength= Misc.allMuscleList(1:Misc.nAllMuscList); % Names of muscles of which tendon stifness is estimated
Misc.lb_lTs_scaling= 0.5; % Lower bound for scaling generic tendon stiffness
Misc.ub_lTs_scaling = 2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_slack_length = {}; % Couple muscles that should have equal tendons stiffness

% select muscles
Misc.MuscleNames_Input = []; % select muscles % [] means all muscles included

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option

% Weights
Misc.wEMG   = 2.5; % weight on tracking EMG
Misc.wAct   = 0.2; % weight on reducing activations
Misc.wTres  = 14;  % weight on residual torques
Misc.wVm    = 0.002; % weight on slack term

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 0;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 0;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0; 	

Misc.Mesh_Frequency = 25;


%% Run muscle tendon parameter estimator
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
opt_fiber_length_scale_PredSim = [Misc.allMuscleList',num2cell(Results.Param.lMo_scaling_paramopt)]';
opt_fiber_length_scale_PredSim = opt_fiber_length_scale_PredSim(:)';
save(fullfile(Misc.OutPath,['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst_optimal_fiber_length_scale.mat']),'opt_fiber_length_scale_PredSim');

% write tendon stiffness estimates in the form that PredSim can read
tendon_slack_length_scale_PredSim = [Misc.allMuscleList',num2cell(Results.Param.lTs_scaling_paramopt)]';
tendon_slack_length_scale_PredSim = tendon_slack_length_scale_PredSim(:)';
save(fullfile(Misc.OutPath,['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst_tendon_slack_length_scale.mat']),'tendon_slack_length_scale_PredSim');

% write tendon stiffness estimates in the form that PredSim can read
tendon_stiff_scale_PredSim = [Misc.allMuscleList',num2cell(Results.Param.kT_scaling_paramopt)]';
tendon_stiff_scale_PredSim = tendon_stiff_scale_PredSim(:)';
save(fullfile(Misc.OutPath,['BCN_' Misc.subjectName '_' Misc.AnalysisID '_paramEst_tendon_stiffness_scale.mat']),'tendon_stiff_scale_PredSim');

