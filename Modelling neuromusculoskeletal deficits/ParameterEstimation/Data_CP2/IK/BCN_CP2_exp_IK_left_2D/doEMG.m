clear
close all
clc

subName = 'CP16';
trialsGS = {'08' '11'};

% trialsnonGS = {'13';...
%     '16';...
%     'squat_T3';...
%     's2s_T3';...
%     'cmj_T2';...
%     'Pendulum_Sit_ISO_r_T1';...
%     'Pendulum_Sit_ISO_l_T5'};
trialsnonGS = {};
trials = [trialsGS trialsnonGS'];

IC{1} = [1.90 3.00 4.09 5.23];
IC{2} = [1.67 2.80 3.83 5.01];
% 
% IC{1} = [1.70 2.85 4.03];
% IC{2} = [1.70 2.86 4.02];
% IC{3} = [1.32 2.46 3.56];
% IC{4} = [1.12 2.26 3.42];
% IC{5} = [1.98 2.94];
% IC{6} = [1.72 2.67];
% IC{7} = [0.75 2.68];
% IC{8} = [1 1.64];
% IC{9} = [0.67 2.94];
% IC{10} = [2.3 7.5];
% IC{11} = [5.03 13.95];

nGC = 6;
DataPath = pwd;

%%
gs = 0;
for f=1:length(trials)
    clear trialID
    trialID = [subName '_T0_' trials{f} '_filtEMG'];
%     trialID = ['EMG_filt_' subName '_T0_' trials{f}];
    clear EMG
    EMG = read_motionFile_v40([DataPath '\' trialID '.mot']);
    if f==1
        EMGlabel_temp = EMG.labels(2:EMG.nc);
    end
    for g = 1:length(IC{f})-1
        gs = gs + 1;
        clear timSpace
        timSpace = linspace(IC{f}(g),IC{f}(g+1),101);
        clear EMGgs_temp
        EMGgs_temp = interp1(EMG.data(:,1),EMG.data(:,2:EMG.nc),timSpace);
        for c=1:EMG.nc-1
            for c2 = 2:EMG.nc
                if strcmp(EMGlabel_temp{c},EMG.labels{c2})
                    EMGgs{gs}(:,c) = EMGgs_temp(:,c2-1);
                end
            end
        end
    end
end

EMGcorrespondence = {'Time','Time';...
    'RREF','rect_fem_r';...
    'RVAL','vas_lat_r';...
    'RBIF','bifemlh_r';...
    'RMEH','semiten_r';...
    'RTIA','tib_ant_r';...
    'RGAS','med_gas_r';...
    'RSOL','soleus_r';...
    'RGLU','glut_med1_r';...
    'LREF','rect_fem_l';...
    'LVAL','vas_lat_l';...
    'LBIF','bifemlh_l';...
    'LMEH','semiten_l';...
    'LTIA','tib_ant_l';...
    'LGAS','med_gas_l';...
    'LSOL','soleus_l';...
    'LGLU','glut_med1_l'};
for i=1:length(EMGlabel_temp)
    EMGlabel{i} = EMGcorrespondence{find(ismember(EMGcorrespondence(:,1),EMGlabel_temp{i})),2};
end

figure
for c = 1:EMG.nc-1
    subplot(4,4,c)
    for g=1:length(EMGgs)
        EMGc{c}(:,g) = EMGgs{g}(:,c);
    end
    EMGc{c} = EMGc{c}/max(max(EMGc{c}));
    EMGgs_exp_mean(:,c) = mean(EMGc{c}(:,1:nGC),2);
    EMGgs_exp_std(:,c) = std(EMGc{c}(:,1:nGC),[],2);
    curve1 = [EMGgs_exp_mean(:,c) + 2*EMGgs_exp_std(:,c)]';
    curve2 = [EMGgs_exp_mean(:,c) - 2*EMGgs_exp_std(:,c)]';
    x=linspace(0,100,101);
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    filling = fill(x2, inBetween, [0.75 0.75 0.75]);
    filling.EdgeColor = [0.75 0.75 0.75];
    filling.FaceAlpha = 0.2;
    hold on
    plot(EMGgs_exp_mean(:,c),'k')
    hold on
%     for g=1:length(EMGgs)
    for g=1:nGC
        plot(EMGc{c}(:,g))
        hold on
    end
    title(EMGlabel{c})
end

save(['EMG_gaitcycle_' subName '.mat'],'EMGgs_exp_mean','EMGgs_exp_std','EMGlabel')
