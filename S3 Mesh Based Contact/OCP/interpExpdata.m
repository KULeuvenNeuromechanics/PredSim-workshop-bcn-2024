function expdata = interpExpdata(S,model_info,d)
% Interpolate experimental data according to collocation scheme
% Author: Gil Serrancolí
% Last edit: September 2024

N = S.solver.N_meshes; % number of mesh intervals
nq = model_info.nq;
coordinate_names = model_info.coord_names.all;

Qs_exp_raw = getIK(S.expdata.IKdata,model_info);
Qs_exp.allfilt = Qs_exp_raw.allfilt;
Qs_exp.time = Qs_exp_raw.time;
Qs_exp.colheaders = Qs_exp_raw.colheaders;
time_IC = [S.initial_time,S.final_time];

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs_exp.allfilt));
Qs_spline.data(:,1) = Qs_exp.allfilt(:,1);
Qdots_spline.data = zeros(size(Qs_exp.allfilt));
Qdots_spline.data(:,1) = Qs_exp.allfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs_exp.allfilt));
Qdotdots_spline.data(:,1) = Qs_exp.allfilt(:,1);
for i = 2:size(Qs_exp.allfilt,2)
    Qs_exp.datafiltspline(i) = spline(Qs_exp.allfilt(:,1),Qs_exp.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs_exp.datafiltspline(i),Qs_exp.allfilt(:,1),1);
end

%identify initial f0 and final ff frames
dt=median(diff(Qs_exp_raw.time));
f0=find((S.initial_time<Qs_exp_raw.time+dt/2)&(S.initial_time>=Qs_exp_raw.time-dt/2));
ff=find((S.final_time<Qs_exp_raw.time+dt/2)&(S.final_time>=Qs_exp_raw.time-dt/2));

for i=1:nq
    coordinate = coordinate_names{i};
    coord_idx = i;
    expdata.IKdata.Qs_all(:,coord_idx) = Qs_spline.data(:,strcmp(Qs_exp.colheaders(1,:),coordinate));
    expdata.IKdata.Qdots_all(:,coord_idx) = Qdots_spline.data(:,strcmp(Qs_exp.colheaders(1,:),coordinate));
    expdata.IKdata.Qdotdots_all(:,coord_idx) = Qdotdots_spline.data(:,strcmp(Qs_exp.colheaders(1,:),coordinate));
end
    
% Interpolation
Qs_time = Qs_spline.data(:,strcmp(Qs_exp.colheaders(1,:),'time'));
time_expi.Qs(1) = f0;
time_expi.Qs(2) = ff;
step = (Qs_time(time_expi.Qs(2))-Qs_time(time_expi.Qs(1)))/(N);
interval = Qs_time(time_expi.Qs(1)):step:Qs_time(time_expi.Qs(2));

expdata.IKdata.Qs = interp1(round(Qs_time,4),expdata.IKdata.Qs_all,S.tgrid);
expdata.IKdata.Qdots = interp1(round(Qs_time,4),expdata.IKdata.Qdots_all,S.tgrid);
expdata.IKdata.Qdotdots = interp1(round(Qs_time,4),expdata.IKdata.Qdotdots_all,S.tgrid);

f_col=1:length(S.tgrid_ext);
f_col(1:(d+1):end)=[];
expdata.IKdata.Qs_col = interp1(round(Qs_time,4),expdata.IKdata.Qs_all,S.tgrid_ext(f_col));
expdata.IKdata.Qdots_col = interp1(round(Qs_time,4),expdata.IKdata.Qdots_all,S.tgrid_ext(f_col));
expdata.IKdata.Qdotdots_col = interp1(round(Qs_time,4),expdata.IKdata.Qdotdots_all,S.tgrid_ext(f_col));

%% Filter and interpolate contact forces
dt=median(diff(S.expdata.ContactForces.data(:,1)));
fs=1/dt;
[B,A]=butter(3,6/(fs/2));
KCF_exp_filt=filtfilt(B,A,S.expdata.ContactForces.data(:,2:3));
expdata.ContactForces=interp1(S.expdata.ContactForces.data(:,1),KCF_exp_filt,S.tgrid);
expdata.ContactForces_col=interp1(S.expdata.ContactForces.data(:,1),KCF_exp_filt,S.tgrid_ext(f_col));


end