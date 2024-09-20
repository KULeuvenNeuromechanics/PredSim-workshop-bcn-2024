function [Tau_pass] = getLimitTorque(coord_name, q, scale)
% --------------------------------------------------------------------------
% getLimitTorque
%   Returns the passive torque attributed to limit torques at q:
%       Tau_pass = K_pass(1)*exp(K_pass(2)*(q - theta_pass(2))) + ...
%                  K_pass(3)*exp(K_pass(4)*(q - theta_pass(1)))
%
%   Coefficient values are taken from: 
%       Anderson III, Frank Clayton. A dynamic optimization solution for a 
%       complete cycle of normal gait. The University of Texas at Austin, 1999.
%
% INPUT:
%   - coord_name -
%   * string of coordinate 
%
%   - q -
%   * angle at which to evaluate the limit torques
%
%   - scale -
%   * scaling factor for limit torque amplitude
%
% OUTPUT:
%   - Tau_pass -
%   * passive limit torque at angle q 
%
% Original author: Lars D'Hondt (PredSim)
% Original date: 12/April/2022
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 02/June/2023
% --------------------------------------------------------------------------

%% Predefined 
K.hip_flexion = [-2.44 5.05 1.51 -21.88];
theta.hip_flexion = [-0.6981 1.81];

K.hip_adduction = [-0.03 14.94 0.03 -14.94];
theta.hip_adduction = [-0.5 0.5];

K.hip_rotation = [-0.03 14.94 0.03 -14.94];
theta.hip_rotation = [-0.92 0.92];

K.knee_angle = [-6.09 33.94 11.03 -11.33];
theta.knee_angle = [-2.4 0.13];

K.ankle_angle = [-2.03 38.11 0.18 -12.12];
theta.ankle_angle = [-0.74 0.52];

K.subtalar_angle = [-60.21 16.32 60.21 -16.32];
theta.subtalar_angle = [-0.65 0.65];

K.mtp_angle = [-0.9 14.87 0.18 -70.08];
theta.mtp_angle = [0 65/180*pi];

K.lumbar_extension = [-0.35 30.72 0.25 -20.36];
theta.lumbar_extension = [-0.5235987755982988 0.17];

K.lumbar_bending = [-0.25 20.36 0.25 -20.36];
theta.lumbar_bending = [-0.3490658503988659 0.3490658503988659];

K.lumbar_rotation = [-0.25 20.36 0.25 -20.36];
theta.lumbar_rotation = [-0.3490658503988659 0.3490658503988659];


%% calculate passive torque
% names of coordinates with predefined bounds
coeff_names = fieldnames(K);

% field index corresponding to 
idx = find(strcmp(coeff_names,coord_name));

K_pass = K.(coeff_names{idx});
theta_pass = theta.(coeff_names{idx});

q = q*pi/180;

Tau_pass = K_pass(1)*scale*exp(K_pass(2)*(q - theta_pass(2))) + ...
           K_pass(3)*scale*exp(K_pass(4)*(q - theta_pass(1)));

end