function [f_lMT_vMT_dM, model_info,coordinates] = generatePolynomials(subject_name, osim_path, pathPredSim)
% --------------------------------------------------------------------------
% generatePolynomials
%   Generate polynomials to describe musculoskeletal geometry of a model.
%
%
% INPUT:
%   - subject_name -
%   * name of the
% 
%   - osim_path -
%   * Full path to opensim model
%
%   - pathPredSim -
%   * Full path to PredSim repository
%
% OUTPUT:
%   - f_lMT_vMT_dM -
%   * casadi function to evaluate the polynomial. 
%       [lMT, vMT, dM] = f_lMT_vMT_dM(Qs, Qdots)
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Lars D'Hondt
% Original date: 13 September 2024

% Last edit by: Ellis Van Can
% Last edit date: September 13, 2024
% --------------------------------------------------------------------------

addpath(pathPredSim)
addpath(fullfile(pathPredSim,'DefaultSettings'))

[S] = initializeSettings('gait1018');
S.subject.name = subject_name;

% pass something to required settings, doesn't really matter what
S.misc.save_folder = pwd;
S.solver.IG_selection = 'quasi-random';

% Settings that are not specified get their default value
S = getDefaultSettings(S,osim_path);

addpath([S.misc.main_path '\VariousFunctions'])
addpath([S.misc.main_path '\PreProcessing'])
addpath([S.misc.main_path '\CasadiFunctions'])

% run preprocessing
[S,model_info] = PreProcessing(S,osim_path);

% create casadi function for polynomials
f_lMT_vMT_dM = createCasadi_MSKGeometry(S,model_info);

% get model coordinates
import org.opensim.modeling.*
model = Model(osim_path);
coordinateSet = model.getCoordinateSet();
numCoordinates = coordinateSet.getSize();
coordinates = cell(1,numCoordinates);
for i = 0:numCoordinates-1
    coordinate = coordinateSet.get(i);
    coordinates{i+1} = char(coordinate.getName());
end
end
