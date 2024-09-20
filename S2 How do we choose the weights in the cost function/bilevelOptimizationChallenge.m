clear; close all; clc;

% This script identifies the optimal weights for the five objective 
% functionals used in 'solveMOCPchallenge' so as to obtain a control policy
% that minimizes the maximum (absolute) value of all joint torques of the
% biped model over time.
%
% Ideal and nadir points are not provided here since the objective
% functionals are intrinsically well scaled.
%
% Author: Alessio Artoni

%% Bilevel optimization (inverse optimal control): weight identification

% The derivative-free 'patternsearch' algorithm is used. Other algorithms,
% ideally for global derivative-free optimization, could be tested.
% An initial guess of 0.2 is used for each weight.
% Function evaluations are parallelized (needs Parallel Computing Toolbox).
% 'MeshTolerance' indicates the (approximate) level of accuracy sought 
% in the final solution.
% We have one bound constraint for each weight (between 0 and 1).
% We have one linear equality constraint: w1+w2+w3+w4+w5 = 1.
w_guess = 0.2*ones(1,5);
options = optimoptions("patternsearch",'Display','iter','PlotFcn',...
          @psplotbestf,'UseParallel',true,MaxIterations=150,...
          MeshTolerance=0.01)

% Objective function
fun = @(w)solveMOCPchallenge(w,false,false);

% Compute solution
w_sol = patternsearch(@(w) fun(w).max_torque,w_guess,[],[],ones(1,5),[1],...
                      zeros(1,5),1.0*ones(1,5),options);

sol = solveMOCPchallenge(w_sol,true,true);

fprintf(['\nThe identified weights are w1 = %.4f, w2 = %.4f, w3 = %.4f, '...
         'w4 = %.4f, and w5 = %.4f.\n'], w_sol(1), w_sol(2), w_sol(3),...
         w_sol(4), w_sol(5));






















