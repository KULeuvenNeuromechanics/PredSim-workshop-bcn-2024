function [] = plot_figures(S)
% --------------------------------------------------------------------------
% plot_figures
%   This functions plots the main figures of the results. The tracking of
%   kinematics, the tracking of knee contact forces (medial and lateral)
%   and make a video of the movement with the level of pressures.

import casadi.*

options.video = 0; % Create .mp4 video of the animations
options.vidname = 'testvid1';

% load results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'R');
time=R.time;

%% Kinematics
figure(1);
d=R.S.solver.d;
N=S.solver.N_meshes;
Qs_exp_all(1:(d+1):N*(d+1)+1,:)=R.S.expdata.IKdata.Qs;
for i=1:N
    for j=1:d
        Qs_exp_all((i-1)*(d+1)+1+j,:)=R.S.expdata.IKdata.Qs_col((i-1)*d+j,:);
    end
end
Qs_exp_all(:,1:3)=Qs_exp_all(:,1:3)*180/pi;
Qs_mod_all_rad= R.kinematics.q_opt_unsc_all.rad;
Qs_mod_all=Qs_mod_all_rad;
Qs_mod_all(:,1:3)=Qs_mod_all_rad(:,1:3)*180/pi;
for i=1:size(Qs_exp_all,2)
    subplot(2,3,i)
    plot(R.time.coll, Qs_mod_all(:,i),'r','LineWidth',2);
    hold all;
    if i==1
        plot(R.time.coll,Qs_exp_all(:,i),'--b','Linewidth',2);
    end
    title(strrep(R.colheaders.coordinates{i},'_',' '));
    if i<=3
        ylabel('angle [º]');
    else
        ylabel('translation [m]');
    end
    xlabel('time [s]')
    if i==1
        legend('model','reference','orientation','horizontal','box','off','Position',[0.137,0.885,0.207,0.043]);
    elseif i==2
        legend('model','orientation','horizontal','box','off','Position',[0.527,0.887,0.094,0.043]);
    end
end
set(gcf,'Position',[79.7,117.7,860,420]);

%% Contact forces
figure(2);
contactforces_exp_all(1:(d+1):N*(d+1)+1,:)=R.S.expdata.ContactForces;
for i=1:N
    for j=1:d
        contactforces_exp_all((i-1)*(d+1)+1+j,:)=R.S.expdata.ContactForces_col((i-1)*d+j,:);
    end
end

subplot(1,2,1);
plot(R.time.coll, contactforces_exp_all(:,1),'LineWidth',2);
hold all;
ylabel('Force [N]');
xlabel('time [s]');
tgrid_ext=R.time.coll;
tgrid_col=tgrid_ext;
tgrid_col(1:(d+1):end)=[];
plot(tgrid_col, R.ContactForces(:,1),'LineWidth',2);
ylim([0 1200]);
title('Medial force')

subplot(1,2,2);
plot(R.time.coll, contactforces_exp_all(:,2),'LineWidth',2);
hold all;
ylabel('Force [N]');
xlabel('time [s]');
tgrid_ext=R.time.coll;
tgrid_col=tgrid_ext;
tgrid_col(1:(d+1):end)=[];
plot(tgrid_col, R.ContactForces(:,2),'LineWidth',2);
ylim([0 1200]);
legend('experimental','model','Orientation','horizontal','box','off');
set(gcf,'Position',[370.3,212.3,806,306])
title('Lateral force');

%% Contact pressures visualization




%% Compute rotation of the tibia with respect to the femur

%Function for rotation 3x3
psi=MX.sym('psi',1);
theta=MX.sym('theta',1);
phi=MX.sym('phi',1);
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];

R=R1*R2*R3;
fR3x3=Function('fR',{psi,theta,phi},{R});

%function for translation 4x4
x=MX.sym('x',1);
y=MX.sym('y',1);
z=MX.sym('z',1);
Rtrans=[1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1];
ftrans=Function('ftrans',{x,y,z},{Rtrans});

%function for translation 0.042 4x4
Rtrans0042=[1 0 0 0; 0 1 0 0.042; 0 0 1 0; 0 0 0 1];

%function for translation and change of orientation 4x4
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R_aux=R1*R2*R3;
Rrottrans=MX(4,4);
Rrottrans(1:3,1:3)=R_aux;
Rrottrans(4,4)=1;
Rrottrans(2,4)=0;
frottrans=Function('ftrans',{psi,theta,phi},{Rrottrans});

%function for rotation 4x4
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R_aux=R1*R2*R3;
Rrot=MX(4,4);
Rrot(1:3,1:3)=R_aux;
Rrot(4,4)=1;
fR=Function('fR',{psi,theta,phi},{Rrot}); %4x4 matrix
f_Raux=Function('f_Raux',{x,y,z,psi,theta,phi},{R_aux}); %only the 3x3 matrix

%translation in tibia coord system
aux=R*[-x;-y;-z];
Rtranstib_4x4=MX(4,4);
Rtranstib_4x4(1,1)=1;
Rtranstib_4x4(2,2)=1;
Rtranstib_4x4(3,3)=1;
Rtranstib_4x4(4,4)=1;
Rtranstib_4x4(1:3,4)=aux;

Rtib=Rtranstib_4x4*Rtrans0042*Rrot;
fTransf=Function('fTransf',{x,y,z,psi,theta,phi},{Rtib});

if options.video

    F = struct('cdata', cell(1,size(pmat,2)), 'colormap', cell(1,size(pmat,2)));
    vidt = 0;

end

centerfem=[0 0.042 0];

fem=stlread('.\Subjects\KneeProsthesis\Femoral_Component_simpler.stl');
if strcmp(S.options.geom_mesh,'100x188')
    tib=stlread('.\Subjects\KneeProsthesis\Tibial Insert simpler_superior_100.stl');
elseif strcmp(S.options.geom_mesh,'49x188')
    tib=stlread('.\Subjects\KneeProsthesis\Tibial_Insert_simpler_superior.stl');
else
    %no other cases implemented
    keyboard;
end

femPoints=fem.Points+centerfem;
pointsfem=find(((femPoints(:,2))>0.00375)&((femPoints(:,2))<0.05));
pointstib=find((tib.Points(:,2)>0.00375)&(tib.Points(:,2)<0.04));
nfacesFem=size(fem.ConnectivityList,1);
nfacesTib=size(tib.ConnectivityList,1);

[m,Ifem]=sort(fem.ConnectivityList(:,1));
facesFem=fem.ConnectivityList(Ifem,:);
[m,Itib]=sort(tib.ConnectivityList(:,1));
facesTib=tib.ConnectivityList(Itib,:);

facesFem_1=[];
facesTib_1=[];
k=1;

for i=1:nfacesFem
    j=1;
    found=false;
    while ~found & (j<size(pointsfem,1))
        if any(facesFem(i,:)==pointsfem(j))
            facesFem_1(k,:)=facesFem(i,:);
            k=k+1;
            found=true;
        end
        j=j+1;
    end
end
k=1;
for i=1:nfacesTib
    j=1;
    found=false;
    while ~found & (j<size(pointstib,1))
        if any(facesTib(i,:)==pointstib(j))
            facesTib_1(k,:)=facesTib(i,:);
            k=k+1;
            found=true;
        end
        j=j+1;
    end
end

facesnot_fem = zeros(size(facesFem,1)-size(facesFem_1,1),3);
facesnot_tib = zeros(size(facesTib,1)-size(facesTib_1,1),3);
LIAfem = ismember(facesFem,facesFem_1,'rows');
LIAtib = ismember(facesTib,facesTib_1,'rows');
knot = 1;

for i=1:size(facesFem,1)
    if LIAfem(i) == 0
        facesnot_fem(knot,:) = facesFem(i,:);
        knot = knot+1;
    end
end
knot = 1;
for i=1:size(facesTib,1)
    if LIAtib(i) == 0
        facesnot_tib(knot,:) = facesTib(i,:);
        knot = knot+1;
    end
end

fem1 = triangulation(vertcat(facesFem_1,facesnot_fem),fem.Points+centerfem);
tib1 = triangulation(vertcat(facesTib_1,facesnot_tib),tib.Points);

for i=1:size(Qs_mod_all,1)

    pos_or=Qs_mod_all_rad   (i,[4:6, 1:3]);

    %kinematics tibia model
    Mtransformtib=full(fTransf(pos_or(1),pos_or(2),pos_or(3),-pos_or(4),-pos_or(5),pos_or(6)));
    for j=1:size(tib.Points,1)
        pointsTib(j,:,i)=Mtransformtib*[tib.Points(j,:) 1]';
    end
end


% Calculate pressures and colors to color using log scale
%Get values of the pressures from the solution
pmat = zeros(size(Qs_mod_all_rad,1), nfacesTib);
F1=external('F', ['Subjects\KneeProsthesis\' strrep(S.misc.external_function,'.dll','_forDebug.dll')]);

for i = 1:size(Qs_mod_all_rad,1)
    out1 = full(F1(Qs_mod_all_rad(i,:)));
    for j = 1:nfacesTib
        if out1(j+2) < 1
            pmat(i,j) = 1;
        else
            pmat(i,j) = out1(j+2);
        end
    end
    FMed(i)=full(out1(1));
    FLat(i)=full(out1(2));
end

colormat = zeros(size(pmat));
% colormat = log10(pmat)./log10(max(pmat, [], 'all'));
colormat = (pmat)./(max(pmat, [], 'all'));

%% save pmat
load(Outname,'R');
R.pmat=pmat;
save(Outname,'R','-append')

%% Draw the femoral and tibial components with the faces of the tibial 
% components in red if there was significant contact, otherwise in green
fcon = figure(3); 
fcon.Color = 'w';
fcon.Position = [0 0 840 1080];

for t=1:size(Qs_mod_all,1)
    clf
    hold on
    patch('Faces',fem1.ConnectivityList,'Vertices',fem1.Points,'Facecolor',[0.3 0.3 0.3],'Facealpha',0.2);
    for i = 1:nfacesTib
        face_i = tib1.ConnectivityList(i,:);
        points_i = pointsTib(face_i,1:3,t);
        patch('Faces',[1 2 3],'Vertices',points_i,...
            'Facecolor',[colormat(t,i) 1-colormat(t,i) 0],'Facealpha',1,'EdgeColor', [67 124 23]./255);
    end
    view([-195 -60]);
    axis equal
    axis off

    title(['time ', num2str(round(time.coll(t),2))]);

    pos_or=Qs_mod_all_rad   (t,1:6);
    Mtransformtib=full(fTransf(pos_or(1),pos_or(2),pos_or(3),-pos_or(4),-pos_or(5),pos_or(6)));

    drawnow;

    if options.video == 1
        vidt = vidt + 1;
        F(vidt) = getframe(gcf);
    end

end

if options.video == 1
    vid = VideoWriter([options.vidname, '.mp4'],'MPEG-4');
    vid.FrameRate = 5;
    open(vid)
    writeVideo(vid,F)
    close(vid)
end