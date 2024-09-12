clear; close all; clc;

% This script starts by estimating the ideal and nadir objective vectors
% for the two objective functionals in 'solveMOCP':
%   - J1 = time integral of the squared 2-norm of joint torques
%   - J2 = time integral of the squared 2-norm of angular accelerations
% Ideal and nadir vectors are needed to normalize these two objectives (and
% properly solve an MOCP).
%
% Then, a user-specified number of points belonging to the Pareto front
% are computed and visualized. Parallelization is used to improve speed 
% with the 'parfor' loop (requires the Parallel Computing Toolbox).
%
% Author: Alessio Artoni

%% Specify the number of Pareto points to be computed

nParetoPoints = 50;


%% Estimate IDEAL and NADIR points (objective vecors)

disp('Estimating ideal and nadir points...');

idealPoint=[0;0];
nadirPoint=[0;0];

% Solve OCP with weights [1,0] (to obtain the minimum of J1)
sol = solveMOCP([1,0],false,false);
idealPoint(1)=sol.J1_original_opt;
nadirPoint(2)=sol.J2_original_opt;

% Solve OCP with weights [0,1] (to obtain the minimum of J2)
sol = solveMOCP([0,1],false,false);
idealPoint(2)=sol.J2_original_opt;
nadirPoint(1)=sol.J1_original_opt;

%% Explore the Pareto front
% Obtain a number of points belonging to the Pareto front using equally
% spread weights.

fprintf('Calculating the requested %d points on the Pareto front...\n', nParetoPoints);

normalized_ParetoPoints = zeros(nParetoPoints,2);
original_ParetoPoints = zeros(nParetoPoints,2);
weights = zeros(nParetoPoints,2);

% Here we use parallelization with 'parfor' (for speed)
parfor k=1:nParetoPoints
    w1 = (k-1)/(nParetoPoints-1);
    weights(k,:) = [w1,1-w1];
    sol = solveMOCP(weights(k,:),false,false,idealPoint,nadirPoint);
    normalized_ParetoPoints(k,:) = [sol.J1_normalized_opt, sol.J2_normalized_opt];
    original_ParetoPoints(k,:) = [sol.J1_original_opt, sol.J2_original_opt];
end

% Plot the normalized Pareto points
figure()
hold on
p1 = plot(normalized_ParetoPoints(2:end-1,1), normalized_ParetoPoints(2:end-1,2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
p2 = plot(normalized_ParetoPoints([1 end],1), normalized_ParetoPoints([1 end],2), 's', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
p3 = plot(0,0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
p4 = plot(1,1, 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
legend([p1, p2, p3, p4], {'Pareto-optimal objective vectors', 'Individual minima', ...
    'Ideal objective vector', 'Nadir objective vector'}, 'Location', 'northeast');
axis equal
axis tight
xlabel('J1 (normalized)');
ylabel('J2 (normalized)');
title('Points on the Pareto front (normalized objective space)');
grid on

% Plot the Pareto points in their original scaling (units of measurement)
figure()
hold on
p5 = plot(original_ParetoPoints(2:end-1,1), original_ParetoPoints(2:end-1,2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
p6 = plot(original_ParetoPoints([1 end],1), original_ParetoPoints([1 end],2), 's', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
p7 = plot(idealPoint(1), idealPoint(2), 'd', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
p8 = plot(nadirPoint(1), nadirPoint(2), 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
legend([p5, p6, p7, p8], {'Pareto-optimal objective vectors', 'Individual minima', ...
    'Ideal objective vector', 'Nadir objective vector'}, 'Location', 'northeast');
%axis equal
axis tight
xlabel('J1  [(Nm)^2s]  (original units)');
ylabel('J2  [s^{-3}]  (original units)');
title('Points on the Pareto front (non-normalized objective space)');
grid on

disp('Done.');
