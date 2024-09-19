clear
close all
clc

subName = 'CP16';
trials = {'08' '11'};

timeLengths = [1.17 5.79;...
    1.01 4.37];

IC{1} = [1.32 2.42 4.59 5.74]; % Need to skip the second cycle
IC{2} = [1.06 2.17 3.26 4.32];

direction = {'p' 'p'};

modelName = 'BCN_CP2_PredSim_noSphere.osim';
DataPath = pwd;

%%
for f = 1:length(trials)
    clear trialID
    trialID = [subName '_T0_' trials{f}];
    time(1,1) = timeLengths(f,1);
    time(1,2) = timeLengths(f,2);
    clear IKsetup
    IKsetup = xmlread('generic_IK_Setup.xml');
    IKsetup.getElementsByTagName('InverseKinematicsTool').item(0).getAttributes.item(0).setValue(trialID);
    IKsetup.getElementsByTagName('time_range').item(0).getFirstChild.setNodeValue([num2str(time(1,1)) ' ' num2str(time(1,2))]);
    IKsetup.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue([DataPath '\' modelName]);
    IKsetup.getElementsByTagName('marker_file').item(0).getFirstChild.setNodeValue([DataPath '\' trialID '_gapFilled.trc']);
%     IKsetup.getElementsByTagName('output_motion_file').item(0).getFirstChild.setNodeValue([DataPath '\' trialID '_IK_long.mot']);
    IKsetup.getElementsByTagName('output_motion_file').item(0).getFirstChild.setNodeValue([DataPath '\' trialID '_IK.mot']);
    xmlwrite([DataPath '\Setup_' trialID '_IK.xml'], IKsetup);
    clear comIK
    comIK=['opensim-cmd run-tool Setup_' trialID '_IK.xml'];
    system(comIK)
end
% %%
% for f = 1:length(trials)
%     clear trialID
%     trialID = [subName '_T0_' trials{f}];
%     time(1,1) = timeLengths(f,1);
%     time(1,2) = timeLengths(f,2);
%     
%     clear ELsetup
%     ELsetup = xmlread('EL_generic_setup.xml');
%     ELsetup.getElementsByTagName('ExternalLoads').item(0).getAttributes.item(0).setValue(trialID);
%     ELsetup.getElementsByTagName('datafile').item(0).getFirstChild.setNodeValue([DataPath '\' trialID '.mot']);
%     ELsetup.getElementsByTagName('external_loads_model_kinematics_file').item(0).getFirstChild.setNodeValue([DataPath '\' trialID '_IK.mot']);
%     xmlwrite([DataPath '\ExternalLoads_' trialID '.xml'], ELsetup);
%     
%     clear IDsetup
%     IDsetup = xmlread('ID_generic_setup.xml');
%     IDsetup.getElementsByTagName('InverseDynamicsTool').item(0).getAttributes.item(0).setValue(trialID);
%     IDsetup.getElementsByTagName('results_directory').item(0).getFirstChild.setNodeValue([DataPath '\']);
%     IDsetup.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue([DataPath '\' modelName]);
%     IDsetup.getElementsByTagName('time_range').item(0).getFirstChild.setNodeValue([num2str(time(1,1)) ' ' num2str(time(1,2))]);
%     IDsetup.getElementsByTagName('external_loads_file').item(0).getFirstChild.setNodeValue([DataPath '\ExternalLoads_' trialID '.xml']);
%     IDsetup.getElementsByTagName('coordinates_file').item(0).getFirstChild.setNodeValue([DataPath '\' trialID '_IK.mot']);
%     IDsetup.getElementsByTagName('output_gen_force_file').item(0).getFirstChild.setNodeValue([trialID '_IK_ID.sto']);
%     xmlwrite([DataPath '\Setup_' trialID '_IK_ID.xml'], IDsetup);
%     clear comID
%     comID=['opensim-cmd run-tool Setup_' trialID '_IK_ID.xml'];
%     system(comID)
% end

%%
gs = 0;
for f=1:length(trials)
    clear trialID
    trialID = [subName '_T0_' trials{f}];
    clear IK
    IK = read_motionFile_v40([DataPath '\' trialID '_IK.mot']);
    for g = 1:length(IC{f})-1
        if g==2 && f==1
        else
            gs = gs + 1;
            clear timSpace
            timSpace = linspace(IC{f}(g),IC{f}(g+1),100);
            IKgs{gs} = interp1(IK.data(:,1),IK.data(:,2:IK.nc),timSpace);
%             IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1) = IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1) - IKgs{gs}(1,find(ismember(IK.labels,'pelvis_tz'))-1);
            IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1) = IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1) - IKgs{gs}(1,find(ismember(IK.labels,'pelvis_tx'))-1);
            if strcmp(direction{f},'n')
                IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1) = -1*IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1);
                IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1) = -1*IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1);
                IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tilt'))-1) = -1*IKgs{gs}(:,find(ismember(IK.labels,'pelvis_tilt'))-1);
                IKgs{gs}(:,find(ismember(IK.labels,'pelvis_list'))-1) = -1*IKgs{gs}(:,find(ismember(IK.labels,'pelvis_list'))-1);
                IKgs{gs}(:,find(ismember(IK.labels,'pelvis_rotation'))-1) = 180+IKgs{gs}(:,find(ismember(IK.labels,'pelvis_rotation'))-1);
            end
        end
    end
end
% gs = 0;
% for f=1:length(trials)
%     clear trialID
%     trialID = [subName '_T0_' trials{f}];
%     clear ID
%     ID = read_motionFile_v40([DataPath '\' trialID '_IK_ID.sto']);
%     for g = 1:length(IC{f})-1
%         gs = gs + 1;
%         clear timSpace
%         timSpace = linspace(IC{f}(g),IC{f}(g+1),100);
%         IDgs{gs} = interp1(ID.data(:,1),ID.data(:,2:ID.nc),timSpace);
% %         IDgs{gs}(:,find(ismember(ID.labels,'pelvis_tz'))-1) = IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1) - IDgs{gs}(1,find(ismember(IK.labels,'pelvis_tz'))-1);
% %         IDgs{gs}(:,find(ismember(ID.labels,'pelvis_tx'))-1) = IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1) - IDgs{gs}(1,find(ismember(IK.labels,'pelvis_tx'))-1);
% %         if strcmp(direction{f},'n')
% %             IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1) = -1*IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tz'))-1);
% %             IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1) = -1*IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tx'))-1);
% %             IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tilt'))-1) = -1*IDgs{gs}(:,find(ismember(IK.labels,'pelvis_tilt'))-1);
% %             IDgs{gs}(:,find(ismember(IK.labels,'pelvis_list'))-1) = -1*IDgs{gs}(:,find(ismember(IK.labels,'pelvis_list'))-1);
% %             IDgs{gs}(:,find(ismember(IK.labels,'pelvis_rotation'))-1) = 180+IDgs{gs}(:,find(ismember(IK.labels,'pelvis_rotation'))-1);
% %         end        
%     end
% end
%%
for c = 1:IK.nc-1
    figure
    for g=1:length(IKgs)
        IKc{c}(:,g) = IKgs{g}(:,c);
        plot(IKgs{g}(:,c))
        hold on
    end
    IKgs_exp_mean(:,c) = mean(IKc{c}(:,[1 3:end]),2);
    IKgs_exp_std(:,c) = std(IKc{c}(:,[1 3:end]),[],2);
    curve1 = [IKgs_exp_mean(:,c) + 2*IKgs_exp_std(:,c)]';
    curve2 = [IKgs_exp_mean(:,c) - 2*IKgs_exp_std(:,c)]';
    x=linspace(0,100,100);
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    filling = fill(x2, inBetween, [0.75 0.75 0.75]);
    filling.EdgeColor = [0.75 0.75 0.75];
    filling.FaceAlpha = 0.2;
    hold on
    plot(IKgs_exp_mean(:,c))
    hold on
    IKlabel{c} = IK.labels{c+1};
    title(IKlabel{c})
end
% for c = 1:ID.nc-1
%     figure
%     for g=1:length(IDgs)
%         IKc{c}(:,g) = IDgs{g}(:,c);
%         plot(IDgs{g}(:,c))
%         hold on
%     end
%     IDgs_exp_mean(:,c) = mean(IKc{c},2);
%     IDgs_exp_std(:,c) = std(IKc{c},[],2);
%     curve1 = [IDgs_exp_mean(:,c) + 2*IDgs_exp_std(:,c)]';
%     curve2 = [IDgs_exp_mean(:,c) - 2*IDgs_exp_std(:,c)]';
%     x=linspace(0,100,100);
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     filling = fill(x2, inBetween, [0.75 0.75 0.75]);
%     filling.EdgeColor = [0.75 0.75 0.75];
%     filling.FaceAlpha = 0.2;
%     hold on
%     plot(IDgs_exp_mean(:,c))
%     hold on
%     IDlabel{c} = ID.labels{c+1};
%     title(IDlabel{c})
% end
save(['IK_gaitcycle_' subName '_right.mat'],'IKgs_exp_mean','IKgs_exp_std','IKlabel')
    
    