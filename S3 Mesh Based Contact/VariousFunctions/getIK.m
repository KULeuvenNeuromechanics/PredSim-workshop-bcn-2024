function Qs_IK = getIK(IKdata,model_info)
% --------------------------------------------------------------------------
% getIK
%   This function returns the coordinate positions (radian) given the path 
%   to the inverse kinematics file and the joints of interest.
%   
% INPUT:
%   - path_IK -
%   * path to the inverse kinematics file
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - Qs_IK -
%   * struct with inverse kinematics results (time, coordinates and labels)
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
%   updated to work for generic set of coordinates
% Last edit by: Gil Serrancolí
% Last edit date: 05/09/2024
% --------------------------------------------------------------------------

% Get the names of the coordinates
coordinate_names = model_info.coord_names.all;

% Load the reference date
Qsall = IKdata;

%% Build the struct with inverse kinematics
% First column has the timestamps
Qs_IK.time = Qsall.data(:,strcmp(Qsall.colheaders,{'time'}));
Qs_IK.all(:,1) = Qs_IK.time;
Qs_IK.colheaders{1,1} = 'time';

% Add the IK of every joint, in the defined order
for i = 1:length(coordinate_names)
    coordinate = coordinate_names{i};
    idx_coord_Qsall = find(contains(Qsall.colheaders,[coordinate '/value']),1,'first');
    if isempty(idx_coord_Qsall)
        Qs_IK.(coordinate) = zeros(size(Qs_IK.time));
    else
        % if strcmp(Qsall.inDeg,'yes')
        %     if sum(model_info.ExtFunIO.jointi.translations(:) == model_info.ExtFunIO.coordi.(coordinate))
        %         Qs_IK.(coordinate) = Qsall.data(:,idx_coord_Qsall);
        %     else
        %         Qs_IK.(coordinate) = Qsall.data(:,idx_coord_Qsall).*(pi/180);
        %     end
        % else
        %data is assumed to be in radians
            Qs_IK.(coordinate) = Qsall.data(:,idx_coord_Qsall);
        % end

    end
    Qs_IK.all(:,i+1) = Qs_IK.(coordinate);
    Qs_IK.colheaders{1,i+1} = coordinate;
end

% Low-pass filter
order = 4;
cutoff_low = 6;
fs=1/mean(diff(Qs_IK.all(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qs_IK.allfilt = Qs_IK.all;
Qs_IK.allfilt(:,2:end) = filtfilt(af,bf,Qs_IK.allfilt(:,2:end));

end