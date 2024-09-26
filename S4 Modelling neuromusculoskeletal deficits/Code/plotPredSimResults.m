clear
close all
clc

%% Required settings
subject = 'CP1'; % subject name
paramEstModelName = {['BCN_' subject,'_v1_paramEst'], ['BCN_' subject,'_v1_paramEst']}; % names of the models
paramEstSuffix =  {'job39' 'job42'};
modelLegend = {'all MT parameters estimated analysis 1'; 'all MT parameters estimated analysis 2'}; % what the user wants to call each result.
pathSeminar = 'C:\GBW_MyPrograms\PredSim-workshop-bcn-2024\S4 Modelling neuromusculoskeletal deficits';
predsimResultsPath = {'C:\GBW_MyPrograms\PredSimResults'};

%% Setup
pathData = fullfile(pathSeminar,'Data');
modelLegends = [{'Generic - no personalization'},modelLegend{:}];
originalModelPredSimResultPath = {[pathData '\Data_' subject '\Model\PredSimResult']};
IKdataPath = [pathData '\Data_' subject '\IK'];
if strcmp(subject,'CP1')
    originalSuffix = 'job30';
    mass = 58.8;
    subNameCorrespondence = 'CP4';
elseif strcmp(subject,'CP2')
    originalSuffix = 'job23';
    mass = 54.4;
    subNameCorrespondence = 'CP16';
end

colors = [[0 0 0];winter(length(paramEstModelName))];
joints = {'pelvis_tilt', 'pelvis_tx', 'pelvis_ty', ...
    'hip_flexion_l', 'hip_flexion_r', ...
    'knee_angle_l', 'knee_angle_r', ...
    'ankle_angle_l', 'ankle_angle_r', ...
    'lumbar_extension'};
jointLabels = joints;
subName = ['BCN_' subject];
models = [subName paramEstModelName];
jobs = [originalSuffix paramEstSuffix];
resultPath = originalModelPredSimResultPath;
for i=1:length(paramEstModelName)
    resultPath = [resultPath predsimResultsPath];
end
addpath(genpath(predsimResultsPath{1}))

exp_IK_r = load(fullfile(IKdataPath,[subName '_exp_IK_right_2D/IK_gaitcycle_' subNameCorrespondence '_right.mat']));
exp_IK_l = load(fullfile(IKdataPath,[subName '_exp_IK_left_2D/IK_gaitcycle_' subNameCorrespondence '_left.mat']));
exp_2_r = exp_IK_r;
exp_2_l = exp_IK_l;
exp_2_r.IKlabel_right = exp_2_r.IKlabel;
exp_2_l.IKlabel_left = exp_2_l.IKlabel;
exp_2_r.IKgs_exp_mean_right = exp_2_r.IKgs_exp_mean;
exp_2_l.IKgs_exp_mean_left = exp_2_l.IKgs_exp_mean;
exp_2_r.IKgs_exp_std_right = exp_2_r.IKgs_exp_std;
exp_2_l.IKgs_exp_std_left = exp_2_l.IKgs_exp_std;

for m=1:length(models)
    clear R
    IK(m) = read_motionFile_v40([resultPath{m} '\' models{m} '\' models{m} '_' jobs{m} '.mot']);
    load([resultPath{m} '\' models{m} '\' models{m} '_' jobs{m} '.mat']);
    Rall{m} = R;
    ID(m).data = R.kinetics.T_ID;
    ID(m).labels = R.colheaders.coordinates;
    MIall{m} = model_info;
end

%% plot result
figure
p=0;
for j=1:length(joints)
    p=p+1;
    subplot(2,5,p)
    for c=1:length(exp_2_r.IKlabel_right)
        if strcmp(exp_2_r.IKlabel_right{c},joints{j})
            if strcmp(joints{j}(end-1:end),'_l')
                    meanCurve = exp_2_l.IKgs_exp_mean_left(:,c);
                    curve1 = [exp_2_l.IKgs_exp_mean_left(:,c) + 2*exp_2_l.IKgs_exp_std_left(:,c)]';
                    curve2 = [exp_2_l.IKgs_exp_mean_left(:,c) - 2*exp_2_l.IKgs_exp_std_left(:,c)]';
                else
                    meanCurve = exp_2_r.IKgs_exp_mean_right(:,c);
                    curve1 = [exp_2_r.IKgs_exp_mean_right(:,c) + 2*exp_2_r.IKgs_exp_std_right(:,c)]';
                    curve2 = [exp_2_r.IKgs_exp_mean_right(:,c) - 2*exp_2_r.IKgs_exp_std_right(:,c)]';
            end
            x=linspace(0,100,100);
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            filling = fill(x2, inBetween, [0.25 0.25 0.25]);
            filling.EdgeColor = [0.25 0.25 0.25];
            filling.FaceAlpha = 0.5;
            hold on
        end
    end
    hold on
    for m=1:length(models)
        if strcmp(joints{j}(end-1:end),'_l')
            idx_HS_l = find(diff(Rall{m}.ground_reaction.GRF_l(:,2)>(mass/3))==1);
            idx_l = [idx_HS_l:100 1:(idx_HS_l-1)];
            jointCurve = IK(m).data(idx_l,find(ismember(IK(m).labels,joints{j})));
        else
            jointCurve = IK(m).data(1:end/2,find(ismember(IK(m).labels,joints{j})));
        end
        plot(linspace(0,100,100),jointCurve,'Color',colors(m,:),'LineWidth',1.5);
        RMSE(j,m) = sqrt(sumsqr(meanCurve - jointCurve)/size(meanCurve,1));
        hold on
    end
    xlim([0 100])
    if p>5
        xlabel('% gait cycle')
    end
    if p==1 || p==6
        ylabel('Angle (deg)')
    end
    box off
    title(jointLabels{j},'interpreter','none')
    if j==4
        legend([{[subject ' experimental data']},modelLegends],'interpreter','none','FontSize',10)
    end
end
