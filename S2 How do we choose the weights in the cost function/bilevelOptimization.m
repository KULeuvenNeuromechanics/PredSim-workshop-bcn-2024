clear; close all; clc;

% This script identifies the optimal weights for the two objective 
% functionals (J1,J2) used in 'solveMOCP' so as to best match the    
% "experimental" walking kinematic data of the biped model. The kinematic 
% data (segment angles) are stored in the provided file 'walkingData.mat'.
%
% This technique is an inverse optimal control method, designed to 
% determine the optimal control policy (i.e., the weights for the 
% selected cost functionals) that best replicates the experimental
% motion data.
%
% Author: Alessio Artoni

%% Use pre-computed ideal and nadir points and load kinematic data

% Use pre-computed ideal and nadir points (by 'exploreParetoFront.m') for
% the predictive simulation of walking of the planar five-link biped model:
idealPoint = [239.8502; 19.2975];
nadirPoint = [6321.1; 192.578];

% Load the "experimental" kinematic data (absolute angles of the model's 
% segments over time) obtained with a specific "motor contol policy" (i.e, 
% unknown, to-be-determined weights for the two cost functionals of the
% model):
load('walkingData.mat'); % Contains a structure called 'walkingData'

%% Bilevel optimization (inverse optimal control): weight identification

% We identify the weights for the two cost functionals (hence the Pareto 
% point) by minimizing the deviations between predicted model kinematics
% and experimental kinematic data (see function 'kinDeviation').

% The derivative-free 'patternsearch' algorithm is used. Other algorithms,
% ideally for global derivative-free optimization, could be tested.
% An initial guess of 0.5 is used for each weight.
% Function evaluations are parallelized (needs Parallel Computing Toolbox).
% 'MeshTolerance' indicates the (approximate) level of accuracy sought
% in the final solution.
% We have one bound constraint for each weight (between 0 and 1).
% We have one linear equality constraint: w1 + w2 = 1.
w_guess = [0.5,0.5];
options = optimoptions("patternsearch",'Display','iter','PlotFcn',...
          @psplotbestf,'UseParallel',true,MaxIterations=150,...
          MeshTolerance=0.1);
      
% Objective function
fun = @(w)kinDeviation(walkingData,w,idealPoint,nadirPoint);

% Compute solution
w_sol = patternsearch(fun,w_guess,[],[],[1,1],[1],[0,0],[1,1],options);

fprintf('\nThe identified weights are w1 = %.4f and w2 = %.4f.\n',...
        w_sol(1), w_sol(2));

% Note that this specific problem is actually one-dimensional and has just
% one bound constraint.
% w1_guess = 0.5;
% fun = @(w1)kinDeviation(walkingData,[w1,1-w1],idealPoint,nadirPoint);
% w1_sol = patternsearch(fun,w1_guess,[],[],[],[],0,1,options);

%% Animate and plot the results

% Animate predicted motion
sol = solveMOCP(w_sol,true,false,idealPoint,nadirPoint);

% Create 2D matrices with 5 columns containing the segment angles in deg
expMat = 180/pi*[walkingData.q1', walkingData.q2', walkingData.q3',...
                 walkingData.q4', walkingData.q5'];
predMat = 180/pi*[sol.q1_opt', sol.q2_opt', sol.q3_opt', sol.q4_opt',...
                  sol.q5_opt'];

% Plot predicted vs. experimental segment angles (absolute, not relative)
figure()
plotHandles = [];
for i = 1:5
    % Create a subplot in a 1-row, 5-column grid
    subplot(1, 5, i); % 1 row, 5 columns, ith subplot   
    % Plot predicted and experimental data
    hPred=plot(sol.time,predMat(:, i), 'b-', 'DisplayName', 'Predicted');
    hold on
    hExp=plot(sol.time,expMat(:, i), 'r--', 'DisplayName', 'Experimental');    
    % Store handles for legend (only need to store once)
    if i == 1
        plotHandles = [hPred, hExp];
    end    
    xlabel('Time [s]'); 
    % Add a title for each subplot
    title(['q' num2str(i) ' [deg]']);    
    hold off; % Release current plot for the next subplot
end
% Add common legend outside the subplots
legend(plotHandles, {'Predicted', 'Experimental'},...
       'Orientation', 'horizontal', 'Position', [0.5 0.02 0 0]);
% Desired figure size
width = 1000;
height = width*9/16;
% Get screen size (returns vector [x, y, screenWidth, screenHeight])
screenSize = get(0, 'ScreenSize');
% Calculate x and y coordinates to center the figure
x = (screenSize(3) - width) / 2;  % Center horizontally
y = (screenSize(4) - height) / 2; % Center vertically
% Set the figure position and size
set(gcf, 'Position', [x, y, width, height]);















