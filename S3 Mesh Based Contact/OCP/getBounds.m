function [bounds_nsc] = getBounds(S,model_info)
% --------------------------------------------------------------------------
% getBounds
%   This script provides bounds for the (not scaled) optimisation variables.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - bounds_nsc -
%   * boundaries for all optimisation variables (not scaled)
%
% 
%
% Last edit by: Gil Serrancolí
% Last edit date: September 2024
% --------------------------------------------------------------------------

% Get the names of the coordinates
coordinate_names = model_info.coord_names.all;
NCoord = model_info.nq;

%% Initialise
coords = ["Qs","Qdots","Qdotdots"];
for i=coords
    bounds_nsc.(i).lower = nan(1,NCoord);
    bounds_nsc.(i).upper = nan(1,NCoord);
end

%% Lower and upper bounds
if S.options.scale_kneety_av
    bounds_nsc.Qs.lower=[-1.3 -0.1 -0.2 -0.02 -1 -0.02];
    bounds_nsc.Qs.upper=[ 0.2  0.1  0.2  0.02  1  0.02];
else
    bounds_nsc.Qs.lower=[-1.3 -0.1 -0.2 -0.02 0.03 -0.02];
    bounds_nsc.Qs.upper=[ 0.1  0.1  0.2  0.02 0.06  0.02];
end

bounds_nsc.Qdots.lower=-10*ones(1,6);
bounds_nsc.Qdots.upper= 10*ones(1,6);

bounds_nsc.Qdotdots.lower=[-1000 -100 -100 -10 -10 -10];
bounds_nsc.Qdotdots.upper=[ 1000  100  100  10  10  10];

% Not needed for the workshop:
bounds_nsc.residualforces.lower=-1000;
bounds_nsc.residualforces.upper= 1000;

bounds_nsc.uT.lower=-200*ones(1,6);
bounds_nsc.uT.upper= 200*ones(1,6);

bounds_nsc.E = 1e3;
bounds_nsc.E = 1e9;

end % end of function