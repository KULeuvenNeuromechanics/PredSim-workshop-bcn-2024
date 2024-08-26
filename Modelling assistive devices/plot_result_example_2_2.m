


clear
close all
clc

% create empty arrays for results and legend names
results = {};
legendNames = {};

% add a result and its legend name to the end of the array
results{end+1} = '.\Results\example_2\gait1018\gait1018_v7.mat';
legendNames{end+1} = 'baseline (running)';

results{end+1} = '.\Results\example_2\gait1018\gait1018_v8.mat';
legendNames{end+1} = 'rubber band';


% create figure
fig1 = figure;
tiledlayout('flow')
colours = [0,0,0; parula(length(results))];

% loop over results
for i=1:length(results)

    load(results{i})
    
    figure(fig1);
    
    % kinematics
    nexttile(1)
    hold on
    plot(R.kinematics.Qs(:, strcmp(R.colheaders.coordinates, 'hip_flexion_r')),...
        'Color',colours(i,:), 'LineWidth',2, 'DisplayName',legendNames{i})
    
    if i==1
        ylabel('[°]')
        xlabel('[% gait cycle]')
        title('hip flexion')
        lg=legend('FontSize',12, 'Orientation','horizontal');
        lg.Layout.Tile = 'south';
    end
    
    nexttile(2)
    hold on
    plot(R.kinematics.Qs(:, strcmp(R.colheaders.coordinates, 'knee_angle_r')),...
        'Color',colours(i,:), 'LineWidth',2)
    
    if i==1
        ylabel('[°]')
        xlabel('[% gait cycle]')
        title('knee extension')
    end
    
    nexttile(3)
    hold on
    plot(R.kinematics.Qs(:, strcmp(R.colheaders.coordinates, 'ankle_angle_r')),...
        'Color',colours(i,:), 'LineWidth',2)
    
    if i==1
        ylabel('[°]')
        xlabel('[% gait cycle]')
        title('ankle dorsiflexion')
    end

    % kinetics
    nexttile(4)
    hold on
    plot(R.orthosis.combined.T_bio(:, strcmp(R.colheaders.coordinates, 'hip_flexion_r')),...
        'Color',colours(i,:), 'LineWidth',2)
    
    if i==1
        ylabel('[Nm]')
        xlabel('[% gait cycle]')
        title('hip flexion')
    end
    
    nexttile(5)
    hold on
    plot(-R.orthosis.combined.T_bio(:, strcmp(R.colheaders.coordinates, 'knee_angle_r')),...
        'Color',colours(i,:), 'LineWidth',2)
    
    if i==1
        ylabel('[Nm]')
        xlabel('[% gait cycle]')
        title('knee extension')
    end
    
    nexttile(6)
    hold on
    plot(R.orthosis.combined.T_bio(:, strcmp(R.colheaders.coordinates, 'ankle_angle_r')),...
        'Color',colours(i,:), 'LineWidth',2)
    
    if i==1
        ylabel('[Nm]')
        xlabel('[% gait cycle]')
        title('ankle dorsiflexion')
    end

    % GRF
    nexttile(7)
    hold on
    plot(R.ground_reaction.GRF_r(:,2), 'Color',colours(i,:), 'LineWidth',2)
    
    if i==1
        ylabel('[N]')
        xlabel('[% gait cycle]')
        title('vertical GRF')
    end

    % rubber band
    if isfield(R.orthosis,'separate')

        % Reconstruct full gait cycle for variables in R.orthosis.separate
        %   By default, PredSim will assume left-right symmetry and simulate only a
        %   half gait cycle (see S.misc.gaitmotion_type). During post-processing,
        %   the full gait cycle is reconstructed based on automatically detected
        %   symmetry of the model. Because left-right symmetry of orthoses is not
        %   automatically detected, R.orthosis.separate only contains the values
        %   for the half gait cycle.
        %   rubberBandBetweenAnkles is symmetric to itself, so the gait cycle can
        %   be reconstructed by repeating the half gait cycle.
        %   Note that the result was rearranged to start at right heel-strike.
        rubberBand_length = R.orthosis.separate{1, 1}.rubberBand_length;
        % indices of 1st half gait cycle
        idx1st = find(~isnan(rubberBand_length));
        % corresponding indices of 2nd half cycle
        idx2nd = rem(idx1st + R.S.solver.N_meshes - 1, 2*R.S.solver.N_meshes) + 1;
        % reconstruct full cycle
        rubberBand_length(idx2nd) = rubberBand_length(idx1st);
        
        % idem for force
        rubberBand_force = R.orthosis.separate{1, 1}.rubberBand_force;
        rubberBand_force(idx2nd) = rubberBand_force(idx1st);
        
        nexttile(8)
        hold on
        plot(rubberBand_length*100, 'Color',colours(i,:), 'LineWidth',2)
        ylabel('[cm]')
        xlabel('[% gait cycle]')
        title('rubber band length')
        
        nexttile(9)
        hold on
        plot(rubberBand_force, 'Color',colours(i,:), 'LineWidth',2)
        ylabel('[N]')
        xlabel('[% gait cycle]')
        title('rubber band force')

    end


end



