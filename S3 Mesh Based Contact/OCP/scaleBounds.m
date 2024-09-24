function [bounds] = scaleBounds(S,model_info,bounds_nsc,scaling)
% --------------------------------------------------------------------------
% scaleBounds
%   This script scales bounds for the scaled optimisation variables.
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
%   - scaling -
%   * scale factors for all optimisation variables
%
% OUTPUT:
%   - bounds -
%   * boundaries for all scaled optimisation variables
%
% 
% Original author: Lars D'Hondt
% Original date: 6/April/2023
%
% Last edit by: Gil Serrancolí
% Last edit date: September 2024
% --------------------------------------------------------------------------

bounds = bounds_nsc;

% Qs
if S.options.scale_kneety_av
    bounds.Qs.lower([1:4 6]) = (bounds_nsc.Qs.lower([1:4 6]))./scaling.Qs([1:4 6]);
    bounds.Qs.upper([1:4 6]) = (bounds_nsc.Qs.upper([1:4 6]))./scaling.Qs([1:4 6]);
    bounds.Qs.lower(5) = bounds_nsc.Qs.lower(5);
    bounds.Qs.upper(5) = bounds_nsc.Qs.upper(5);

else
    bounds.Qs.lower = (bounds_nsc.Qs.lower)./scaling.Qs;
    bounds.Qs.upper = (bounds_nsc.Qs.upper)./scaling.Qs;
end
% Qdots
bounds.Qdots.lower = (bounds_nsc.Qdots.lower)./scaling.Qdots;
bounds.Qdots.upper = (bounds_nsc.Qdots.upper)./scaling.Qdots;
% Qdotdots
bounds.Qdotdots.lower = (bounds_nsc.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds_nsc.Qdotdots.upper)./scaling.Qdotdots;

% residual forces 
bounds.residualforces.lower =(bounds_nsc.residualforces.lower)./scaling.residualforces;
bounds.residualforces.upper =(bounds_nsc.residualforces.upper)./scaling.residualforces;

end