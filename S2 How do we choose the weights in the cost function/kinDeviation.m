function dev = kinDeviation(walkingData,weights,ideal,nadir)

% This function calculates the sum of the squared deviations between
% the segment angles of the biped model read from the external MAT-file 
% 'walkingData' and the corresponding angles obtained by solving a 
% predictive MOCP, i.e., using the function 'solveMOCP' with given 
% 'weights' and 'ideal' & 'nadir' objective vectors.
%
% Author: Alessio Artoni

sol = solveMOCP(weights,false,false,ideal,nadir);

dev = sum((sol.q1_opt - walkingData.q1).^2 + ...
          (sol.q2_opt - walkingData.q2).^2 + ...
          (sol.q3_opt - walkingData.q3).^2 + ...
          (sol.q4_opt - walkingData.q4).^2 + ...
          (sol.q5_opt - walkingData.q5).^2);

    