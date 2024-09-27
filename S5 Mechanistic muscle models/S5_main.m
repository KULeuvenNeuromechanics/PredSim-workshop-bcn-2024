clear all; close all; clc
%This hands-on tutorial is an introduction into mechanistic muscle models.
%In the tutorial, you will compare the 2-state crossbridge model of Huxley
%(1957) to existing phenomenological Hill-type muscle models (used in
%PredSim & OpenSim). While mechanistic crossbridge models have not yet been
%incoorporated into the PredSim framework (this is work in progress), this
%tutorial provides custom code to run forward simulations with these
%models. The mechanistic muscle modeling framework is based on van der Zee
%et al., (2024), and the associated CaFaXC muscle model.
%
% Author: Tim J van der Zee
% Contact: tim.vanderzee@kuleuven.be
%
% Last edit: 9/27/2024
%
% make sure to clone both the PredSim-workshop-bcf-2024 repository and the
% CaFaXC repository, and specify their locations below

% [EDIT START]
CaFaXC_folder = uigetdir('','Select CaFaXC folder');
PredSim_folder = uigetdir('','Select PredSim folder');
% [EDIT STOP]

addpath(genpath(CaFaXC_folder))
addpath(genpath(PredSim_folder))

%% Part 0: Hill-type force-velocity
% Here, we will evaluate and plot a force-velocity curve, typically used in PredSim
%
% Task 0: Run this section

% define a force-velocity curve
parms.ce.es = [-0.3183 -8.1492 -0.3741 0.8856];
parms.ce.vmaxrel = 10;
parms.func.fv = @(a, F, Fisom, parms) parms.ce.vmaxrel/(2*parms.ce.es(2)) * (-exp(parms.ce.es(4)/parms.ce.es(1) - F/(parms.ce.es(1)*a*Fisom)) + exp(F/(parms.ce.es(1)*a*Fisom) - parms.ce.es(4)/parms.ce.es(1)) - 2*parms.ce.es(3));

% evaluate the curve
fv.FHill = linspace(0,2, 100);
fv.vHill = parms.func.fv(1,fv.FHill, 1, parms);

% plot
if ishandle(1), close(1); end; figure(1)
figure(1)
plot(fv.vHill, fv.FHill,'linewidth',2); hold on
xlabel('Shortening velocity (L_{opt})')
ylabel('Force (F_{max})')
box off
ylim([0 2])
title('Hill-type force-velocity relation')

%% Part 1: The effect of rate constants on crossbridge distributions
% Here, you will choose crossbridge cycling rate constants and evaluate
% their effect on the crossbridge distribution and its first-order moment
%
% Task 1: Adjust the rate constants in XBparams.CB.f and XBparams.CB.g and re-run this Section. 
%
% Question 1: How do the values of XBparams.CB.f and XBparams.CB.g
% influence what is displayed in Subplots A-D?
%
% Subplot A) Steady-state crossbridge distribution shown for 4 types of velocity conditions.
%
% Subplot B) Crossbridge rate functions as specified by the rate constants. 
% note: three different zones can be distinguished

% Subplot C) Cumulative first-order moment of the crossbridge distribution (A).
%
% Subplot D) First-order moment evaluated across the entire strain interval
% note: there is a dashed line displayed in the background, why would this line be displayed?

% define crossbridge functions and parameters
parms = cfxc.gen_funcs(parms);

% scaling parameters
parms.CB.h = 12e-9; % crossbridge reach [m]
parms.CB.s = 2.6e-6;  % sarcomere length [m]
parms.CB.xi = -5:.01:5; % normalized strain vector

% [EDIT START]
% Adjust the crossbridge rate constants, assuming that parms.CB.g(1) = parms.CB.f
parms.CB.f = 1000*rand(1);
parms.CB.g = [parms.CB.f 1000*rand(1) 1000*rand(1)]; 
% [EDIT STOP]

% code below evaluates and plots the results
if ~exist('k','var')
    k = 0;
end

K = mod(k,4);
ks = repmat([.05 .5; .5 .5; .5 .05; .05 .05],100,1);

if ishandle(K+2), close(K+2), end; figure(K+2)
set(gcf,'name', ['XBparams.CB.f = ',num2str(parms.CB.f), ' - XBparams.CB.g = [', num2str(parms.CB.g),']'])
set(gcf,'units','normalized','position',[ks(K+1,:) .45 .3])

subplot(222)
plot(parms.CB.xi, parms.CB.f_func(parms),'k-','linewidth',2); hold on
plot(parms.CB.xi, parms.CB.g_func(parms),'k--','linewidth',2); box off
xline(0,'k--')
xline(1,'k--')
title('B. Crossbridge rate functions')
ylabel('Rate (s^{-1})'); 
xlabel('Crossbridge strain')
xlim([-1 2])

vsel = (-1:.5:1) * parms.ce.vmaxrel;
colors = [linspace(1,0,length(vsel))' zeros(length(vsel),1) linspace(0,1,length(vsel))'];

for i = 1:length(vsel)
    % calculate CB force velocity
    [~,n] = cfxc.CB_force_velocity(vsel(i), parms);

    subplot(221)
    plot(parms.CB.xi, n,'color',colors(i,:)); hold on
    xlim([-1 2])
    box off
    title('A. Steady-state crossbridge distribution')
    xlabel('Crossbridge strain')
    ylabel('Fraction attached')
    
    subplot(223)
    plot(parms.CB.xi, cumtrapz(parms.CB.xi, n.*parms.CB.xi)*4,'color',colors(i,:)); hold on
    xlim([-1 2])
    box off
    title('C. Cumulative first-order moment')
    xlabel('Crossbridge strain')
    ylabel('Moment')
    ylim([-1 2])
    
    subplot(224)
    plot(vsel(i), trapz(parms.CB.xi, n.*parms.CB.xi)*4,'o','color',colors(i,:),'markerfacecolor',colors(i,:)); hold on
    box off
    title('D. First-order moment')
    xlabel('Shortening velocity (L_{opt}/s)')
    ylabel('Moment')
    ylim([-1 2])
end

subplot(221)
xline(0,'k--')
xline(1,'k--')

subplot(223)
xline(0,'k--')
xline(1,'k--')

subplot(224)
plot(fv.vHill, fv.FHill,'k:','linewidth',2); hold on

k = k+1;
for i = 1:(K+1)
    figure(i+1)
    
    if (K+1) == 1
        subplot(221)
        legend('Fast shortening','Slow shortening','Isometric','Slow lengthening', 'Fast lengthening', 'location','best')
        legend boxoff
        
        subplot(222)
        legend('Attachment','Detachment','location','best') 
        legend boxoff
    end
        
end

%% Part 2: The effect of rate constants on force-velocity 
% Here, you will choose crossbridge cycling rate constants and evaluate
% their effect on the force-velocity relation of the crossbridge model
%
% Task 2: Adjust the rate constants in XBparams.CB.f and XBparams.CB.g and re-run this Section. 
%
% Question 2: How do the values of XBparams.CB.f and XBparams.CB.g
% influence the force-velocity relation?

% [EDIT START]                       
parms.CB.f = 1000*rand(1);
parms.CB.g = [parms.CB.f 1000*rand(1) 1000*rand(1)]; 
% [EDIT STOP]

% calculate CB force velocity
FCB = cfxc.CB_force_velocity(fv.vHill, parms);
RMSE = sum((fv.FHill - FCB).^2);
disp(['RMSE = ', num2str(RMSE)])

% plot
if ishandle(10), close(10);end
figure(10)
subplot(2,3,1:3)
plot(fv.vHill, [fv.FHill; FCB],'linewidth',2); hold on
xlabel('Shortening velocity (L_{opt}/s)')
ylabel('Force (F_{max})')
title('Force-velocity relation')
box off
ylim([0 2])
legend('Hill-type', 'crossbridge', 'location', 'best')
legend boxoff

% sensitivity analysis around the guess
s = logspace(-2,2,10);
titles = {'f_1','g_2','g_3'};
for j = 1:3
    newparms = parms;
    
    for i = 1:length(s)
        if j == 1
            newparms.CB.f = parms.CB.f * s(i);
        elseif j == 2
            newparms.CB.g(2) = parms.CB.g(2) * s(i);
        elseif j == 3
            newparms.CB.g(3) = parms.CB.g(3) * s(i);
        end
        
    FCB = cfxc.CB_force_velocity(fv.vHill, newparms);
    RMSE(j,i) = sum((fv.FHill - FCB).^2);
    end

    subplot(2,3,j+3)
    loglog(s, RMSE(j,:),'linewidth',2)
    xline(1,'k-')
    box off
    title(titles{j})
    xlabel('Scale factor')
    ylabel('RMSE')
end

%% Part 3: The effect of rate constants on force-velocity for distribution-moment model
% Here, you will choose crossbridge cycling rate constants and evaluate
% their effect on the force-velocity relation of the distribution-moment crossbridge model (Zahalak, 1981). 
% Because of differences between the original model and the approximation, we need to adjust the rate constants slightly
% for this model to fit the Hill-type force-velocity relation
% Hint: if you yielded a low RMSE in Part 2, start with rates constants you
% ended up with. 
%
% Task 3: Adjust the rate constants in XBparams.CB.f and XBparams.CB.g and re-run this Section. 
%
% Question 3: What is the main difference between the original and the
% approximated model?

% [EDIT START]                       
parms.CB.f = 1000*rand(1);
parms.CB.g = [parms.CB.f 1000*rand(1) 1000*rand(1)]; 
% [EDIT STOP]

% calculate DM model force-velocity
parms.CB.analytical = 1;
parms.CB.mu = 1;

% evaluate
[~, fv] = cfxc.fit_CB_on_Hill(parms, fv,[]);

% plot
if ishandle(1), close(1); end
figure(1)
plot(fv.vHill(:), [fv.FHill(:) fv.FCB(:,1) fv.FCB(:,3)],'linewidth',2); hold on
xlabel('Shortening velocity (L_{opt})')
ylabel('Force (F_{max})')
box off
ylim([0 2])
title('Hill-type force-velocity relation')

legend('Hill-type','Original crossbridge','Approximated crossbridge','location','best')
legend boxoff

%% Part 4: evaluate force-velocity during imposed joint movement
% The distribution-moment model has been interfaced with a tendon and
% parallel elastic tissue, and connected to a skeletal model. 
% 
% Here, you will impose a knee movement and evaluate the quadriceps forces predicted by the crossbridge model.
% You specifiy a knee joint angular velocity, that results in 
% muscle-tendon complex length changes, which are divided between contractile and
% elastic elements according to their relative stiffnesses.

% Task 4: Adjust the imposed joint angular velocity (omega_knee)
%
% Question 4: What are the difference between the Hill-type model and the
% crossbridge model?

% load some parameters
load('quad_parms.mat','parms')
parms.ce.es = [-0.3183 -8.1492 -0.3741 0.8856];
parms.ce.vmaxrel = 10;
parms.func.fv = @(a, F, Fisom, parms) parms.ce.vmaxrel/(2*parms.ce.es(2)) * (-exp(parms.ce.es(4)/parms.ce.es(1) - F/(parms.ce.es(1)*a*Fisom)) + exp(F/(parms.ce.es(1)*a*Fisom) - parms.ce.es(4)/parms.ce.es(1)) - 2*parms.ce.es(3));

if ishandle(20), close(20); end; figure(20)
figure(20)
color = get(gca,'colororder');

% [EDIT START]
% choose a isokinetic velocity
omega_knee = 1000*rand(1); % imposed angular velocity [deg/s]
parms.CB.f = 1000*rand(1);
parms.CB.g = [parms.CB.f 1000*rand(1) 1000*rand(1)]; 
% [EDIT STOP]

% conditions
phi_knee = 70; % starting angle [deg]
dphi_knee = 30; % angular excursion [deg]
momarm = .05; % assumed muscle moment arm
vmtc_imposed = omega_knee * pi/180 * momarm; % [m/s]
parms.exp.phi = phi_knee;
parms = cfxc.calc_x0(parms);

% stim parameters
parms.exp.stim_type = 'constant';
parms.exp.a = 1;
parms.exp.A = 1;
parms.set.optimum = 1;
parms.set.sim_mtc = 1;
parms.set.odeopt = [];
      
for j = 1:2
    if j == 1
        parms.type = 'Hill-type';
        X0 = [parms.ce.amin parms.exp.l0 parms.exp.lmtc];
    else
        parms.type = 'crossbridge_v2';
        X0 = [parms.ce.amin parms.exp.l0 parms.exp.x0(1) parms.exp.x0(3) parms.exp.lmtc];
    end
    
    tic
    parms.exp.vmtc = 0;  
    [~,x0] = ode45(@cfxc.sim_muscle, [0 .1], X0, parms.set.odeopt, parms);
    
    parms.exp.vmtc = vmtc_imposed;
    tmax = min([abs(dphi_knee/omega_knee) .1]);
    
    [t,x] = ode45(@cfxc.sim_muscle, [0 tmax], x0(end,:), parms.set.odeopt, parms);
    T(j) = toc;

    dx = nan(length(t), size(x,2));
    for i = 1:length(t)
        dx(i,:) = cfxc.sim_muscle([], x(i,:)', parms)';
    end
    
    % calc some variables
    lce = x(:,2);
    vce = dx(:,2);
    lmtc = x(:,end);
    Fpe = parms.func.fpe(lce, parms);
    lse = lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;
    Fce = Fse - Fpe;
    
    figure(20)
    subplot(241)
    plot(t, lce); hold on; box off
    title('CE length')
    xlabel('Time (s)')
    ylabel('Length (m)')
    
    subplot(242)
    plot(t, vce,'color',color(j,:)); hold on; box off   
    title('CE velocity')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    
    subplot(245)
    plot(t, lmtc); hold on; box off
    title('MT length')
    xlabel('Time (s)')
    ylabel('Length (m)')
    
    subplot(246)
    plot(t, Fce); hold on; box off
    title('CE force')
    xlabel('Time (s)')
    ylabel('Force (N)')
    
    subplot(2,4,[3,4,7,8])
    plot(vce/parms.ce.lceopt, Fce/parms.ce.Fmax, '.-','color',color(j,:).^2); hold on; box off
    
    Fss(j) = Fce(end);
    vss(j) = vce(end);
    
end

subplot(2,4,[3,4,7,8])
color = get(gca,'colororder');

% recalc fv
fv.FHill = linspace(0,2, 100);
fv.vHill = parms.func.fv(1,fv.FHill, 1, parms);
[~, fv] = cfxc.fit_CB_on_Hill(parms, fv,[]);

% plot fv
plot(fv.vHill, fv.FHill,'color',color(1,:),'linewidth',2)
plot(fv.vHill, fv.FCB(:,3)', 'color',color(2,:),'linewidth',2);
xlim([min(fv.vHill) max(fv.vHill)])
xlabel('Shortening velocity (L_{opt}/s)')
title('Force-velocity')

plot(vss(1)/parms.ce.lceopt, Fss(1)/parms.ce.Fmax,'o', 'color',color(1,:),'markerfacecolor',color(1,:))
plot(vss(2)/parms.ce.lceopt, Fss(2)/parms.ce.Fmax,'o', 'color',color(2,:),'markerfacecolor',color(2,:))
ylim([0 3])

legend('Hill-type','Crossbridge','location','best')
legend boxoff

                       
                       