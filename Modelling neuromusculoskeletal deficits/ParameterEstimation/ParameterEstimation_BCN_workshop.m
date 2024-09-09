%% Example EMG constrained simulation for ankle, knee and hip in the sagital plane

% In this example we estimate parameters of multiple lower-limb muscles
% using an EMG driven simulation of the ankle-knee and hip.
% Additionally, this example handles different trials with different number
% of degrees of freedom together.
% This file is based on https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver/blob/master/Examples/Example_EMGWalking/EMGconstrained_variableDOF.m

%% clear variables and command window
clear all; clc; close all;

%% Parameters to change for BCN workshop

Misc.subjectName = 'CP1';

% Weights
Misc.wEMG   = 2.5; % weight on tracking EMG
Misc.wAct   = 0.2; % weight on reducing activations
Misc.wTres  = 14;  % weight on residual torques

% upper and lower bound for difference between simulated and measured
% muscle activity (EMG)
Misc.EMGbounds  = [-0.1 0.1]; 

% upper and lower bound for EMG scaling
Misc.BoundsScaleEMG = [0.2 5];

pathParamEst = 'C:\Users\u0145647\OneDrive - KU Leuven\KU Leuven\BCN_Workshop\ParameterEstimation';
pathPredSim = 'C:\GBW_MyPrograms\PredSim_2D_BCN';

Misc.AnalysisID = 'v1'; % Name given to the estimation you run

Misc.casadiPath = 'C:\GBW_MyPrograms\casadi_3_6_5';
Misc.useCluster = 0;
if Misc.useCluster
    Misc.clusterName = 'LocalFourByTwo';
end
Misc.visualStudioVersion = 'Visual Studio 17 2022';

%%
runParamEstAndPredSim