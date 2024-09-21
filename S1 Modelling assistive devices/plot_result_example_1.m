clear
close all
clc

results{1} = './Results/example_1/gait1018/gait1018_v4.mat';
legNames{1} = 'w/o exo';

results{2} = './Results/example_1/gait1018/gait1018_v3.mat';
legNames{2} = 'gain = 40';

results{3} = './Results/example_1/gait1018/gait1018_v5.mat';
legNames{3} = 'gain = 80';

figure
tiledlayout('flow')


for i=1:length(results)

    load(results{i},'R')
    
    % index of right ankle
    idx_ankle = find(strcmp(R.colheaders.coordinates,'ankle_angle_r'));
    
    nexttile(1)
    hold on
    plot(R.kinematics.Qs(:,idx_ankle), 'LineWidth',2)
    ylabel('Angle [°]')
    title('Ankle')
    
    nexttile(2)
    hold on
    plot(R.kinetics.T_ID(:,idx_ankle), 'LineWidth',2)
    ylabel('Moment [Nm]')
    title('Ankle - total')
    
    nexttile(3)
    hold on
    plot(R.orthosis.combined.T_bio(:,idx_ankle), 'LineWidth',2)
    ylabel('Moment [Nm]')
    title('Ankle - bio')
    
    nexttile(4)
    hold on
    plot(R.orthosis.combined.T_coord(:,idx_ankle), 'LineWidth',2)
    ylabel('Moment [Nm]')
    title('Ankle - exo')

    nexttile(5)
    hold on
    plot(R.muscles.a(:,strcmp(R.colheaders.muscles,'soleus_r')), 'LineWidth',2)
    ylabel('Activation[-]')
    title('Soleus')
    xlabel('gait cycle [%]')

    nexttile(6)
    hold on
    plot(R.muscles.a(:,strcmp(R.colheaders.muscles,'tib_ant_r')), 'LineWidth',2)
    ylabel('Activation[-]')
    title('Tib ant')
    xlabel('gait cycle [%]')


end

lg=legend(legNames,'FontSize',12, 'Orientation','horizontal');
lg.Layout.Tile = 'south';


