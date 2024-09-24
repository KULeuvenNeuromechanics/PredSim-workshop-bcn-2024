function [S] = getExperimentalData(S)
% --------------------------------------------------------------------------
% getExperimentalData
%
% This file loads the experimental data of one normal gait trial of the
% repository Grand Challenge Competition to Predict In Vivo Knee Loads 
% (fourth edition), full data available from 
% https://simtk.org/projects/kneeloads
% For further description of the repository and data read>
% BJ Fregly, TF Besier, DG Lloyd, SL Delp, SA Banks, MG Pandy, DD D'Lima. 
% Grand challenge competition to predict in vivo knee loads. J Orthop Res 
% 2012, 30(4):503-13. doi: 10.1002/jor.22023.
%
% Here we use the flexion angle and the inf-sup. forces to be tracked
%
% Author: Gil Serrancolí
% Last edit: September 2024

IKdata=importdata('KneeProsthesis_states2track.sto');

ContactForces=importdata('KCF_ref.mot');

S.expdata.IKdata=IKdata;
S.expdata.ContactForces=ContactForces;
S.expdata.ContactForces.data(:,2:3)=S.expdata.ContactForces.data(:,2:3)*4.4482216153; 
%data was in lbs, convert it to N

S.expdata.Qsall.inDeg='false';








end