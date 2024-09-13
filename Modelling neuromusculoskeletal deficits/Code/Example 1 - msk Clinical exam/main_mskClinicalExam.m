%% Modelling neuromusculoskeletal defictis
% hands-on sessions
% Example 1. Modelling musculoskeletal impairments from the clinical exam
%
%   This script requires PredSim and its dependencies to be installed.
%   See https://github.com/KULeuvenNeuromechanics/PredSim
% 

% Original authors: Bram Van Den Bosch, Ellis Van Can
% Original date: September 13, 2024
clc; clear; close all
%% 1. Add paths 

% ------- start edit -------
PredSim_path = 'C:\GBW_MyPrograms\PredSim'; 
Seminar_path = 'C:\Users\u0167909\OneDrive - KU Leuven\Workshop BCN - Seminar model personalization';
% ------- stop edit -------

addpath(genpath(PredSim_path));
addpath(genpath(Seminar_path));
%% 2. Intialize settings
% Subject CP1 ('BCN_CP1')
% Subject CP2 ('BCN_CP2')
% ------- start edit -------
subject_name = 'BCN_CP1'; 
% ------- stop edit -------

osim_path = fullfile(Seminar_path,'Models',subject_name,[subject_name,'_PredSim.osim']);

% load model info
[f_lMT_vMT_dM, model_info,coordinates] = generatePolynomials(subject_name, osim_path, PredSim_path);
%% 3. Define muscle strength scaling factors for PredSim

% Scale the maximal (active) muscle force based on the clinical exam:
%   CE        scaling factor
%   1           0.1
%   2           0.3
%   3           0.5
%   4           0.7
%   5           1

% EXAMPLE: Clinical exam score 'strength_Hpext_R' = 3 
% will be converted to S.settings.muscle_strength = {{'glut_max_r'},0.7}

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


%% 4. find muscle fiber scaling factors
% We scale optimal fiber length (lMo) as an estimate of muscle contractures (shortening of a
% muscle)

% Assume the eximator applies a torque of 15 Nm
% Then the passive torque around the joint should be 15 Nm at the end of
% the ROM observed during the clinical exam 

% The goal of  is to find the scaling factor (sf) that reflects this. 
% In other words: find the line that goes through the crosssection of -15
% Nm and the clinical exam angle

% To find the sf you have to play a little with sf_lMo (line ...)
% It can be useful to begin with a broader range, such as [0.7:0.1:1], and
% then narrow it down.

% Be aware that 15 NM is a rough estimate of the examinators torque and
% this can vary between examinators and clinical examinations
% Therefore, after a few tries, when you have a line that is very close to the cross section, this is sufficient
% ------- start edit -------

% NOTE 1: First do distal muscles, for 'hamstrings' you also need to specify your
% found sf_lMo for the gastrocs as these also have an effect on knee ROM.
muscle_toScale = 'soleus'; % Options: 'soleus', 'gastrocnemii', % hamstrings
if contains(muscle_toScale,'hamstrings')
    sf_lMo_gastrocnemii = 0.8;
end

% NOTE 2: you have to do this for each side seperately. When clinical exam angle = typical
% value - sf can remain 1 and you dont have to do this estimation
side = 'r'; % Options: 'l', 'r' 
CE_angle = 10 ; % Passive ROM angle Clinical Exam

% Define scaling factor range
sf_lMo = flip([0.8:0.05:0.95]);  

% ------- stop edit -------

% get scaling factor
%(length_subject * massa_subject)/(length_model * massa_model)
if contains(subject_name,'BCN_CP1')
    S.subject.St = 0.994369501;
elseif contains(subject_name,'BCN_CP2')
    S.subject.St = 0.909325513;
end

%%% Scaling settings
% Angle + side
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

% joint index
idx_joint = find(strcmp(coordinates,coord_name_side));

% Create Qs (joint angles) and Qdot (angular velocity)
ROM = linspace((CE_angle-20)*pi/180,(CE_angle+20)*pi/180,25); % Range of motion
n = length(ROM);


Qs = zeros(n,length(coordinates));
Qdots = zeros(n,length(coordinates));
Qs(:,idx_joint) = ROM; 

% Put model to position in CLinical exam if position is other than supine
if model_toDiffPos
Qs(:,strcmp(coordinates,model_changeAngle)) = ones(n,1)*(model_changeDeg*pi/180); %knee_angle_r 90°
end

%%% Scaling factors
% set muscle parameters
close all
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
    
    a = ones(18,1)*0.01; % activation
    fse  = zeros(18,1); % tendon force-length characteristic 
    dfse = 0;
    
    vMT  = 0; % MT velocity
    
    MuscMoAsmp = 0; % constant pennation angle
    d = 0.01; % muscle damping
    
    lMT = zeros(18,n);
    FT = zeros(n,18); 
    M_muscle = zeros(18,n);
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
yline(-5, 'HandleVisibility','off');
yline(-10, 'HandleVisibility','off');
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
%% Lets run predictive simulations with personalized muscle parameters!

