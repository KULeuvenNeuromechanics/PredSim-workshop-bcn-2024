%% Modelling neuromusculoskeletal defictis
% hands-on sessions
% Example 1. Modelling musculoskeletal impairments from the clinical exam
%
%   This script requires PredSim and its dependencies to be installed.
%   See https://github.com/KULeuvenNeuromechanics/PredSim
% 

% Original authors: Bram Van Den Bosch, Ellis Van Can
% Original date: September 20,2024

clc; clear; close all
%% 1. Add paths 

% ------- start edit -------
PredSim_path = 'C:\GBW_MyPrograms\PredSim'; % path to PredSim
Seminar_path = 'C:\GBW_MyPrograms\PredSim-workshop-bcn-2024\S4 Modelling neuromusculoskeletal deficits'; % path to Seminar code
casADiPath = 'C:\GBW_MyPrograms\casadi_3_6_5';
% ------- stop edit -------

addpath(genpath(PredSim_path));
addpath(genpath(Seminar_path));
addpath(genpath(casADiPath));
%% 2. Intialize settings
% Subject CP1 ('BCN_CP1')
% Subject CP2 ('BCN_CP2')

% ------- start edit -------
subject_name = 'BCN_CP1'; 
% ------- stop edit -------

osim_path = fullfile(Seminar_path,'Models',subject_name,[subject_name,'_PredSim.osim']);

% load model info
[f_lMT_vMT_dM, model_info,coordinates] = generatePolynomials(subject_name, osim_path, PredSim_path);

%% In step 3 and 4 you will work with the data from the Clinical Exam (CE)
% Passive range of motion and muscle strength scores are provided in the 'Clinical Exam'/BCN_CP#:

% NOTE: For convenience we put these two steps in this specific order but they can be performed independently
%% 3. Define muscle strength scaling factors for PredSim
% Scale the maximal (active) muscle force based on the strength scores in the CE

% Clinical Exam to strength scaling factor reference
%   CE        scaling factor
%   1           0.1
%   2           0.3
%   3           0.5
%   4           0.7
%   5           1

% EXAMPLE: 
% CE score 'strength_Hpext_R' = 3 
% Code:
% S.settings.muscle_strength = {{'glut_max_r'},0.5}


% ------- start edit -------

S.settings.muscle_strength = {...       
    {'glut_max_r'},1,...                % R_hip_ext
    {'iliopsoas_r'},1,...               % R_hip_flex    
    {'rect_fem_r' 'vasti_r'},1,...      % R_knee_ext
    {'hamstrings_r' 'bifemsh_r'},1,...  % R_knee_bend
    {'gastroc_r' 'soleus_r'},1,...      % R_ankle_pf
    {'tib_ant_r'},1,...                 % R_ankle_df
    {'glut_max_l'},1,...                % L_hip_ext
    {'iliopsoas_l'},1,...               % L_hip_flex    
    {'rect_fem_l' 'vasti_l'},1,...      % L_knee_ext
    {'hamstrings_l' 'bifemsh_l'},1,...  % L_knee_bend
    {'gastroc_l' 'soleus_l'},1,...      % L_ankle_pf
    {'tib_ant_l'},1};                   % L_ankle_df

% ------- stop edit -------


%% 4. Find muscle fiber scaling factors
% We scale optimal fiber length (lMo) as an estimate of muscle contractures (shortening of a muscle)
% You have to make code edits for part 4.1 and 4.2

%%% 4.2 Get the passive range of motion (pROM) scores from the clinical exam

% Do this for each muscle that deviates from typical values (see S4
% Modelling neuromusculoskeletal deficits/README.md)
% If CE value pROM = typical value, sf_lMo will remain 1

% Start with the most distal muscles.

% Do this for each side seperately

% ------- start edit -------

muscle_toScale = 'soleus'; % Options: 'soleus', 'gastrocnemii', % hamstrings
side = 'r'; % Options: 'l', 'r' 
CE_angle = 10 ; % Passive ROM angle Clinical Exam

% ------- stop edit -------

% For scaling of the hamstrings, you also need to provide gastrocs
% scaling factor as these also crosses the knee joint
if contains(muscle_toScale,'hamstrings')
    prompt = {'Enter scaling factors for gastrocnemii:'};
    dlgtitle = 'sf_lMo_gastrocnemii';
    dims = [1 60];
    definput = {'1'};
    sf_lMo_gastrocnemii = inputdlg(prompt, dlgtitle, dims, definput);
end


% Scale stiffness and damping of the joints (St) based on subject dimensions
% For privacy reasons these St's are already provided
% (length_subject * massa_subject)/(length_model * massa_model)
if contains(subject_name,'BCN_CP1')
    S.subject.St = 0.994369501;
elseif contains(subject_name,'BCN_CP2')
    S.subject.St = 0.909325513;
end

% Put the model in Clinical Exam position if this deviates from supine
if strcmp(muscle_toScale,'soleus')
            coord_name = 'ankle_angle';
            model_toDiffPos = 1;
            model_changeAngle = ['knee_angle_',side];
            model_changeDeg = -90;
elseif strcmp(muscle_toScale,'gastrocnemii')
            coord_name = ['ankle_angle'];
            model_toDiffPos = 0;
elseif strcmp(muscle_toScale,'hamstrings')
            coord_name = ['knee_angle'];
            model_toDiffPos = 1;
            model_changeAngle = ['hip_flexion_',side];
            model_changeDeg = 90;        
end

coord_name_side = [coord_name,'_',side];

% Find joint index in coordinates
idx_joint = find(strcmp(coordinates,coord_name_side));

% Create Qs (joint angles) and Qdot (angular velocity)
ROM = linspace((CE_angle-20)*pi/180,(CE_angle+20)*pi/180,25); % Range of motion
n = length(ROM);

Qs = zeros(n,length(coordinates));
Qdots = zeros(n,length(coordinates));
Qs(:,idx_joint) = ROM; 

if model_toDiffPos
Qs(:,strcmp(coordinates,model_changeAngle)) = ones(n,1)*(model_changeDeg*pi/180); %knee_angle_r 90°
end

%%% 4.2 Evaluates scaling factor

% Assume the eximator applies a torque of 15 Nm. 
% Then the passive torque around the joint should be 15 Nm at the end of
% the ROM observed during the clinical exam 

% Find the scaling factor (sf) that reflects this 15 Nm. 
% Then run the code and pick the right scaling factor from the figure 
% (find the line that goes through the crosssection of -15 % Nm and the
% clinical exam angle)

% To find the sf you have to play a little with sf_lMo (line ...)
% It can be useful to begin with a broader range, such as [0.7:0.1:1], and
% then narrow it down.

% NOTE: 15 NM is a rough estimate of the examinators torque and
% this can vary between examinators and clinical examinations

% ------- start edit -------

% Define scaling factor range
sf_lMo = flip([0.8:0.05:0.95]); 

% ------- stop edit -------


% Loop that evaluates the scaling factors
close all % close previous figs
for j = 1:length(sf_lMo)
    if strcmp(muscle_toScale,'soleus')
        scale.subject.scale_MT_params = {{['soleus_',side]},'lMo',sf_lMo(j)};
    elseif strcmp(muscle_toScale,'gastrocnemii')
        scale.subject.scale_MT_params = {{['gastroc_',side]},'lMo',sf_lMo(j)};
    elseif strcmp(muscle_toScale,'hamstrings')
        scale.subject.scale_MT_params = {{['gastroc_',side]},'lMo',sf_lMo_gastrocnemii,...
            {['hamstrings_',side],['bifemsh_',side]},'lMo',sf_lMo(j)};
     %%%% TO DO: decide if we want to include bifemsh in this
    end


    % muscle strength
    if ~isfield(scale.subject,'muscle_strength')
        scale.subject.muscle_strength = [];
    end
    
    % muscle stiffness
    if ~isfield(scale.subject,'muscle_pass_stiff_shift')
        scale.subject.muscle_pass_stiff_shift = [];
    end
    if ~isfield(scale.subject,'muscle_pass_stiff_scale')
        scale.subject.muscle_pass_stiff_scale = [];
    end
    
    % tendon stiffness
    if ~isfield(scale.subject,'tendon_stiff_scale')
        scale.subject.tendon_stiff_scale = [];
    end
   
    % % get strength scaling factors
    if ~isempty(scale)
        try
            % scale
            [model_info.muscle_info] = scale_MTparameters(scale,model_info.muscle_info);
        catch errmsg
            error(['Unable to extract strength scale factor because: ', errmsg.message]);
        end
    end
    
    
    % Solve the muscle-tendon force equilibrium for the given length and
    % activation to find the total force along the tendon.
    
    % n = 1; % number of timepoints
    options = optimset('Display','off');
    
    load Fvparam
    load Fpparam
    load Faparam
    
    FMo_in = struct_array_to_double_array(model_info.muscle_info.parameters,'FMo');
    lMo_in = struct_array_to_double_array(model_info.muscle_info.parameters,'lMo');
    lTs_in = struct_array_to_double_array(model_info.muscle_info.parameters,'lTs');
    alphao_in = struct_array_to_double_array(model_info.muscle_info.parameters,'alphao');
    vMmax_in = struct_array_to_double_array(model_info.muscle_info.parameters,'vMmax');
    tension = struct_array_to_double_array(model_info.muscle_info.parameters,'specific_tension');
    aTendon = struct_array_to_double_array(model_info.muscle_info.parameters,'tendon_stiff');
    shift = struct_array_to_double_array(model_info.muscle_info.parameters,'tendon_stiff_shift');
    stiffness_shift = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_pass_stiff_shift');
    stiffness_scale = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_pass_stiff_scale');
    strength = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_strength');
    
    a = ones(model_info.muscle_info.NMuscle,1)*0.01; % activation
    fse  = zeros(model_info.muscle_info.NMuscle,1); % tendon force-length characteristic 
    dfse = 0;
    
    vMT  = 0; % MT velocity
    
    MuscMoAsmp = 0; % constant pennation angle
    d = 0.01; % muscle damping
    
    lMT = zeros(model_info.muscle_info.NMuscle,n);
    FT = zeros(n,model_info.muscle_info.NMuscle); 
    M_muscle = zeros(model_info.muscle_info.NMuscle,n);
    Tau_pass = zeros(1,n);
    
    for i=1:n
        % evaluate casadi function to get MT lengths and moment arms
        [lMTj,vMTj,MAj] =  f_lMT_vMT_dM(Qs(i,:),Qdots(i,:));
        
        lMT(:,i) = full(lMTj);
        vMTj = full(vMTj);
        MA = full(MAj);
    
        l_MT = lMT(:,i);
        FT_tmp = fsolve('getHilldiffFun_usingCasadiFunction',fse,options,a,dfse,l_MT,vMT,FMo_in,lMo_in,...
                    lTs_in,alphao_in,vMmax_in,Fvparam,Fpparam,Faparam,tension,aTendon,shift,...
                    MuscMoAsmp,d,stiffness_shift,stiffness_scale,strength);
    
        FT(i,:) = FT_tmp.*FMo_in;
        M_muscle(:,i) = FT(i,:)'.*MA(:,idx_joint); % muscle moments
        [Tau_pass(:,i)] = getLimitTorque(coord_name, Qs(i,idx_joint)*180/pi, S.subject.St); % calculate limit torques
    end
    
    
    % calculate moments
   
    M_tot = sum(M_muscle)+Tau_pass;
    
    % plot
    
    f1 = gcf;
    figure(f1)
    plot(Qs(:,idx_joint)*180/pi,M_tot,'DisplayName',['sf lMo at ' num2str(sf_lMo(j)*100) '%']); 
    hold on;
end

title(['passive torque-angle relationship ',strrep(muscle_toScale, '_', ' '),' ', side]);
xline(CE_angle, 'HandleVisibility','off');
yline(-15, 'HandleVisibility','off');
ylabel('Torque (Nm)')
xlabel([strrep(coord_name_side, '_', ' '),' (°)'])
ylim([-20 0])

legend

%% 5. Define your estimated muscle fiber scaling factors for PredSim
% ------- start edit -------
S.subject.scale_MT_params = {{'soleus_r'},'lMo',1,...
                             {'gastroc_r'},'lMo',1,...
                             {'hamstrings_r'},'lMo',1,...
                              {'soleus_l'},'lMo',1,...
                             {'gastroc_l'},'lMo',1,...
                             {'hamstrings_l'},'lMo',1,};
% ------- stop edit -------

%% 6. Save S
save(fullfile(PredSim_path,'Subjects',subject_name,[subject_name,'.mat']),"S");
%% Now you are able to run predictive simulations with personalized muscle parameters!
