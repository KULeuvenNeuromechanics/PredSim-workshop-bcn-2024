function [rubberBand] = rubberBandBetweenAnkles_sol(init, settings_orthosis)
% --------------------------------------------------------------------------
% rubberBandBetweenAnkles
%   Defines a rubber band connecting both ankles.
% 
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'rubberBandBetweenAnkles'  i.e. name of this function           
%       - stiffness:        stiffness in N/m
%       - zero_length:      length @ 0 stiffness force, in m
%   Values are set via S.orthosis.settings{i} in main.m
%
%
% OUTPUT:
%   - rubberBand -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 13/August/2024
% --------------------------------------------------------------------------


% Note:
%   In the calculation below, variables with suffix "_vec" are vectors
%   expressed in the global reference frame.
%   See Modelling assistive devices/Hands-on examples.md/2.2 for a sketch
%   of the vector definitions.

% create Orthosis object
rubberBand = Orthosis('rubberBand',init);

% rubber band wraps around ankle at 40 cm below the knee joint centre
pos_band_wrt_knee = [0,-0.4,0]; % expressed in tibia body frame

% get position of point where rubber band wraps around left ankle
pos_ankle_left_vec = rubberBand.var_point('ankle_left', 'tibia_l', pos_band_wrt_knee); 

% get position of point where rubber band wraps around right ankle
pos_ankle_right_vec = rubberBand.var_point('ankle_right', 'tibia_r', pos_band_wrt_knee); 

% length of the rubber band
band_length_vec = pos_ankle_right_vec - pos_ankle_left_vec;
band_length = norm(band_length_vec);

% elongation of the rubber band
band_elongation = band_length - settings_orthosis.zero_length;

% orientation of the rubber band (unit vector)
band_orientation_vec = band_length_vec / norm(band_length_vec);

% force in the rubber band (magnitude)
band_force = settings_orthosis.stiffness * band_elongation;

% The rubber band only provides tensile forces. Set the force to 0 if the
% elongation is negative.
% If zero_length is below 0.16
if settings_orthosis.zero_length > 0.16
    band_force = band_force.*smoothIf(band_elongation, 5e-3, 0);
end

% force in the rubber band on the left ankle
band_force_on_left_ankle_vec = band_force * band_orientation_vec;

% force in the rubber band on the right ankle
band_force_on_right_ankle_vec = -band_force_on_left_ankle_vec;

% apply force to left ankle
rubberBand.addBodyForce(band_force_on_left_ankle_vec, 'band_force_left',...
    'tibia_l', pos_band_wrt_knee, 'ground');

% apply force to right ankle
rubberBand.addBodyForce(band_force_on_right_ankle_vec, 'band_force_right',...
    'tibia_r', pos_band_wrt_knee, 'ground');

% select intermediate variables that we want to be evaluated during
% post-processing of the solution
rubberBand.addVarToPostProcessing(band_length, 'rubberBand_length')
rubberBand.addVarToPostProcessing(band_force, 'rubberBand_force')



end