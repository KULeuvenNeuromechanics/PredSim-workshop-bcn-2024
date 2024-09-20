clear
close all
clc

subject = 'CP2';
subName = ['BCN_' subject];
subNameCorrespondence = 'CP16';
pathParamEst = 'C:\Users\u0145647\OneDrive - KU Leuven\KU Leuven\BCN_workshop_full\PredSim-workshop-bcn-2024\Modelling neuromusculoskeletal deficits\ParameterEstimation';
predsimResultsPath = {'C:\GBW_MyPrograms\PredSim_2D_BCN\PredSimResults'};
addpath(genpath(predsimResultsPath{1}))
originalModelPredSimResultPath = {[pathParamEst '\Data_' subject '\Model\PredSimResult']};
IKdataPath = [pathParamEst '\Data_' subject '\IK'];
if strcmp(subject,'CP1')
    originalSuffix = 'job30';
    mass = 58.8;
elseif strcmp(subject,'CP2')
    originalSuffix = 'job23';
    mass = 54.4;
end
paramEstModelName = {[subName '_paramEst']};
paramEstSuffix = {'job29'};
modelLegends = {'linearly scaled MT parameters 1stepDown';...
    'all MT parameters estimated'};
colors = {'k','b'};
joints = {'pelvis_tilt', 'pelvis_tx', 'pelvis_ty', ...
    'hip_flexion_l', 'hip_flexion_r', ...
    'knee_angle_l', 'knee_angle_r', ...
    'ankle_angle_l', 'ankle_angle_r', ...
    'lumbar_extension'};
jointLabels = joints;
models = [subName paramEstModelName];
jobs = [originalSuffix paramEstSuffix];
for i=1:length(paramEstModelName)
    resultPath = [originalModelPredSimResultPath predsimResultsPath];
end

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
    if exist([resultPath{m} '\' models{m} '\' models{m} '_' jobs{m} '.mot'])
        IK(m) = read_motionFile_v40([resultPath{m} '\' models{m} '\' models{m} '_' jobs{m} '.mot']);
        load([resultPath{m} '\' models{m} '\' models{m} '_' jobs{m} '.mat'])
        Rall{m} = R;
        ID(m).data = R.kinetics.T_ID;
        ID(m).labels = R.colheaders.coordinates;
        MIall{m} = model_info;
    end
end
cycle = 'FGC';

%%
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
%             plot(linspace(0,100,100),IKgs_exp_mean(:,c),'k','LineWidth',1.5)
        end
    end
    hold on
    for m=1:length(models)
        if exist([resultPath{m} '\' models{m} '\' models{m} '_' jobs{m} '.mot'])
            if strcmp(joints{j}(end-1:end),'_l')
                idx_HS_l = find(diff(Rall{m}.ground_reaction.GRF_l(:,2)>(mass/3))==1);
                idx_l = [idx_HS_l:100 1:(idx_HS_l-1)];
                jointCurve = IK(m).data(idx_l,find(ismember(IK(m).labels,joints{j})));
            else
                jointCurve = IK(m).data(1:end/2,find(ismember(IK(m).labels,joints{j})));
            end
            if strcmp(joints{j},'knee_angle_l')
                jobs{m}
                mean(jointCurve)
            end
            plot(linspace(0,100,100),jointCurve,colors{m},'LineWidth',1.5);
            RMSE(j,m) = sqrt(sumsqr(meanCurve - jointCurve)/size(meanCurve,1));
            hold on
    %         text(40,40,['RMSE = ' num2str(RMSE(j,1))],'Color',colors{m})
        end
    end
    xlim([0 100])
%     xlabel('% gait cycle')
    box off
    title(jointLabels{j})
    if j==4
        legend(['experimental data',modelLegends'],'interpreter','none','FontSize',10)
    end
end

%%
clear x
clear x2
clear curve1
clear curve2
clear inBetween

load(fullfile(IKdataPath,[subName '_exp_IK_right_2D/EMG_gaitcycle_' subNameCorrespondence '.mat']));
figure,
for m=1:9
    subplot(2,5,m)
    for c=1:length(EMGlabel)
        if strcmp(EMGlabel{c},MIall{1}.muscle_info.muscle_names{m})
            curve1 = [EMGgs_exp_mean(:,c) + 2*EMGgs_exp_std(:,c)]';
            curve2 = [EMGgs_exp_mean(:,c) - 2*EMGgs_exp_std(:,c)]';
            x=linspace(0,100,101);
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            filling = fill(x2, inBetween, [0.75 0.75 0.75]);
            filling.EdgeColor = [0.75 0.75 0.75];
            filling.FaceAlpha = 0.5;
            hold on
        end
    end
    hold on
    for mm=1:length(models)
        if exist([resultPath{mm} '\' models{mm} '\' models{mm} '_' jobs{mm} '.mot'])
            idx = find(ismember(MIall{mm}.muscle_info.muscle_names,MIall{1}.muscle_info.muscle_names{m}));
            if ~isempty(idx)
                plot(linspace(0,100,100),Rall{mm}.muscles.a(:,idx),colors{mm})
            end
            hold on
        end
    end
    title(MIall{1}.muscle_info.muscle_names{m},'interpreter','none')
    if m==4
        legend(['experimental data',modelLegends'],'interpreter','none','FontSize',10)
    end
end
if strcmp(cycle,'FGC')
    figure,
    for ml=1:9
        subplot(2,5,ml)
        m=ml+9;
        for c=1:length(EMGlabel)
            if strcmp(EMGlabel{c},MIall{1}.muscle_info.muscle_names{m})
                curve1 = [EMGgs_exp_mean(:,c) + 2*EMGgs_exp_std(:,c)]';
                curve2 = [EMGgs_exp_mean(:,c) - 2*EMGgs_exp_std(:,c)]';
                x=linspace(0,100,101);
                x2 = [x, fliplr(x)];
                inBetween = [curve1, fliplr(curve2)];
                filling = fill(x2, inBetween, [0.75 0.75 0.75]);
                filling.EdgeColor = [0.75 0.75 0.75];
                filling.FaceAlpha = 0.5;
                hold on
            end
        end
        hold on
        for mm=1:length(models)
            if exist([resultPath{mm} '\' models{mm} '\' models{mm} '_' jobs{mm} '.mot'])
                idx = find(ismember(MIall{mm}.muscle_info.muscle_names,MIall{1}.muscle_info.muscle_names{m}));
                if ~isempty(idx)
                    plot(linspace(0,100,100),Rall{mm}.muscles.a(:,idx),colors{mm})
                end
                hold on
            end
        end
        title(MIall{1}.muscle_info.muscle_names{m},'interpreter','none')
        if ml==4
            legend(['experimental data',modelLegends'],'interpreter','none','FontSize',10)
        end
    end
end
