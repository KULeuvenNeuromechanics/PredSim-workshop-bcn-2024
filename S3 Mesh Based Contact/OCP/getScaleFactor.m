function [scaling] = getScaleFactor(S,model_info,bounds_nsc)
% --------------------------------------------------------------------------
% getScaleFactor
%   This script provides scaling factors for the optimisation variables.
%   Scale factors are based on the bound with highest absolute value, so
%   optimisation variables remain within the interval [-1,1].
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - bounds_nsc -
%   * boundaries for all optimisation variables (not scaled)
%
% OUTPUT:
%   - scaling -
%   * scale factors for all optimisation variables
% 
% Original author: Lars D'Hondt
% Original date: 6/April/2023
%
% Last edit by: Gil Serrancolí
% Last edit date: 5/September/2024
% --------------------------------------------------------------------------


coordinate_names = model_info.coord_names.all;
NCoord = model_info.nq;

%% Default: based on bounds
% Qs
scaling.Qs = max(abs(bounds_nsc.Qs.lower),abs(bounds_nsc.Qs.upper));
scaling.knee_ty.a=0.042;
scaling.knee_ty.b=0.007;
% Qdots
scaling.Qdots = max(abs(bounds_nsc.Qdots.lower),abs(bounds_nsc.Qdots.upper));
% Qdotdots
scaling.Qdotdots = max(abs(bounds_nsc.Qdotdots.lower),abs(bounds_nsc.Qdotdots.upper));
% residual forces
scaling.residualforces=1000;