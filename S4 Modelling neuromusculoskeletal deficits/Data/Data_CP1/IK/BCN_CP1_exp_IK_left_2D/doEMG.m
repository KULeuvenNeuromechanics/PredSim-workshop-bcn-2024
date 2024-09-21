clear
close all
clc

subName = 'CP15';
trialsGS = {'08' '09'};
trialsnonGS = {'16';...
    '19';...
    '37';...
    '33';...
    '38';...
    '121';...
    '112';...
    '100';...
    '58';...
    '67';...
    '91'};
trials = [trialsGS trialsnonGS'];

IC{1} = [2.58 3.84 5.13 6.5 7.79];
IC{2} = [1.15 2.51 3.8 5.04 6.26];
IC{3} = [3.79 4.94];
IC{4} = [2.77 3.85];
IC{5} = [1.24 5];
IC{6} = [1.93 6.16];
IC{7} = [2.12 3.07];
IC{8} = [2.98 4.10];
IC{9} = [1.85 4.90];
IC{10} = [0.86 4.7];
IC{11} = [0.06	5.92];
IC{12} = [1.4 7.09];
IC{13} = [0.06 5.2];

nGC = 8;

DataPath = pwd;

%%
gs = 0;
for f=1:length(trials)
    clear trialID
    trialID = [subName '_T0_' trials{f} '_filtEMG'];
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
