function [] = OCP_formulation(S,model_info,f_casadi)
% --------------------------------------------------------------------------
% OCP_formulation
%   This function formulates the OCP and calls the solver
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
% OUTPUT:
%   - This function returns no outputs -
% 
% --------------------------------------------------------------------------

disp('Start formulating OCP...')
t0 = tic;

%% User inputs (typical settings structure)
% settings for optimization
N = S.solver.N_meshes; % number of mesh intervals
W = S.weights; % weights optimization
nq = model_info.nq;

%% Load external functions
import casadi.*
% The external function performs inverse dynamics and computes the pressure 
% at each face of the contact geometry mesh through the OpenSim/Simbody C++ 
% API. This external function is compiled as a dll from which we create a 
% Function instance using CasADi in MATLAB. More details about the external 
% function can be found in the documentation.
pathmain = pwd;

% Loading external functions.
setup.derivatives =  'AD'; % Algorithmic differentiation
cd(S.misc.subject_path)
F  = external('F',S.misc.external_function);
cd(pathmain);

%% Collocation Scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
d = 3; % degree of interpolating polynomial
S.solver.d=d;
method = 'legendre'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Time grid
% Mesh points
tgrid = linspace(S.initial_time,S.final_time,N+1);
dtime = zeros(1,d+1);
for i=1:(d+1)
    dtime(i)=tau_root(i)*((S.final_time-S.initial_time)/N);
end
% Mesh points and collocation points
tgrid_ext = zeros(1,(d+1)*N+1);
for i=1:N
    tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
end
tgrid_ext(end)=S.final_time;
S.tgrid=tgrid;
S.tgrid_ext=tgrid_ext;

%% Get bounds and initial guess

bounds_nsc = getBounds(S,model_info);
scaling = getScaleFactor(S,model_info,bounds_nsc);
bounds = scaleBounds(S,model_info,bounds_nsc,scaling);
S.expdata = interpExpdata(S,model_info,d);

if strcmp(S.subject.IG_selection,'quasi-random')
    %not implemented here
    guess = getGuess_QR_opti(S,model_info,scaling,d);
else
    guess = getGuess_DI_opti(S,model_info,scaling,d,F);
end

%TO DO
% if (S.misc.visualize_bounds)
%     visualizebounds
% end

%% OCP create variables and bounds
% using opti
opti = casadi.Opti();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define states
% Qs at mesh points
Qs = opti.variable(nq,N+1);
opti.subject_to(bounds.Qs.lower'*ones(1,N+1) < Qs < bounds.Qs.upper'*ones(1,N+1));
opti.set_initial(Qs, guess.Qs');
% Qs at collocation points
Qs_col = opti.variable(nq,d*N);
opti.subject_to(bounds.Qs.lower'*ones(1,d*N) < Qs_col < ...
    bounds.Qs.upper'*ones(1,d*N));
opti.set_initial(Qs_col, guess.Qs_col');
% Qdots at mesh points
Qdots = opti.variable(nq,N+1);
opti.subject_to(bounds.Qdots.lower'*ones(1,N+1) < Qdots < ...
    bounds.Qdots.upper'*ones(1,N+1));
opti.set_initial(Qdots, guess.Qdots');
% Qdots at collocation points
Qdots_col = opti.variable(nq,d*N);
opti.subject_to(bounds.Qdots.lower'*ones(1,d*N) < Qdots_col < ...
    bounds.Qdots.upper'*ones(1,d*N));
opti.set_initial(Qdots_col, guess.Qdots_col');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define controls

% Define "slack" controls
% Time derivative of Qdots (states) at collocation points
A_col = opti.variable(nq, d*N);
opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
    bounds.Qdotdots.upper'*ones(1,d*N));
opti.set_initial(A_col, guess.Qdotdots_col');

% Define knee residual forces
if S.options.useResidualforces
    knee_res = opti.variable(2, d*N);
    opti.subject_to(bounds.residualforces.lower*ones(2,d*N) < knee_res < ...
        bounds.residualforces.upper*ones(2,d*N));
    opti.set_initial(knee_res, guess.residualforces');
end

%% OCP: collocation equations
% Define CasADi variables for states
Qsk         = MX.sym('Qsk',nq);
Qsj         = MX.sym('Qsj',nq,d);
Qskj        = [Qsk Qsj];
Qdotsk      = MX.sym('Qdotsk',nq);
Qdotsj      = MX.sym('Qdotsj',nq,d);
Qdotskj     = [Qdotsk Qdotsj];

Qsj_nsc_exp= MX.sym('Qsj_nsc',nq,d); 
KCF_exp      = MX.sym('KCF_exp',2,d);
knee_resj   = MX.sym('knee_resj',2,d);

% if nq.torqAct > 0
%     a_ak        = MX.sym('a_ak',nq.torqAct);
%     a_aj        = MX.sym('a_akmesh',nq.torqAct,d);
%     a_akj       = [a_ak a_aj];
% end
% % Define CasADi variables for controls
% if nq.torqAct > 0
%     e_ak    = MX.sym('e_ak',nq.torqAct);
% end

% Define CasADi variables for "slack" controls
Aj          = MX.sym('Aj',nq,d);

J           = 0; % Initialize cost function
eq_constr   = {}; % Initialize equality constraint vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step
h = (S.final_time-S.initial_time)/N;
% Loop over collocation points
for j=1:d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unscale variables
    if S.options.scale_kneety_av
        Qskj_nsc=MX.zeros(nq,d+1);
        Qskj_nsc([1:4 6],:) = Qskj([1:4 6],:).*(scaling.Qs([1:4 6])'*ones(1,size(Qskj,2)));
        Qskj_nsc(5,:) = Qskj(5,:)*scaling.knee_ty.b+scaling.knee_ty.a;
    else
        Qskj_nsc = Qskj.*(scaling.Qs'*ones(1,size(Qskj,2)));
    end
    Qdotskj_nsc = Qdotskj.*(scaling.Qdots'*ones(1,size(Qdotskj,2)));
    Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));
    knee_resj_nsc = knee_resj.*(scaling.residualforces);

    % Expression for the state derivatives at the collocation points
    Qsp_nsc      = Qskj_nsc*C(:,j+1);
    Qdotsp_nsc   = Qdotskj_nsc*C(:,j+1);
    % if nq.torqAct > 0
    %     a_ap         = a_akj*C(:,j+1);
    % end

    % Append collocation equations
    % Dynamic constraints are scaled using the same scale
    % factors as the ones used to scale the states
    
    % Skeleton dynamics (implicit formulation)
    qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
    eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.Qs';
    eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./scaling.Qdots';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Add contribution to the cost function
    J = J + ...
        W.q_dotdot  * B(j+1)*(f_casadi.uasqsum(Aj(:,j)))*h + ...
        W.trackq    * B(j+1)*(sum((Qskj_nsc(1,j+1)-Qsj_nsc_exp(1,j)).^2))*h+...
        W.minsecq   * B(j+1)*(sum((Qskj_nsc([2 3 4 6],j+1)).^2))*h;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign Qs
    F_ext_input = Qskj_nsc(:,j+1);
    
    % Evaluate external function
    out = F(F_ext_input);
    KCFMed=out(1);
    KCFLat=out(2);

    if strcmp(S.options.forcetrackingAs,'costfun')
        J = J+  W.trackKCF   * B(j+1) * (f_casadi.ftrackforceexp([KCFMed; KCFLat],KCF_exp(:,j)))*h;
    elseif strcmp(S.options.forcetrackingAs,'constraint')
        if S.options.useResidualforces
            eq_constr{end+1}=f_casadi.ftrackforceexp([KCFMed+knee_resj_nsc(1,j); KCFLat+knee_resj_nsc(2,j)],KCF_exp(:,j));
            J = J + W.minresidualforce * B(j+1) * (sum(knee_resj(:,j).^2))*h; 
        else
            eq_constr{end+1}=f_casadi.ftrackforceexp([KCFMed; KCFLat],KCF_exp(:,j));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
end % End loop over collocation points

eq_constr = vertcat(eq_constr{:});

% Casadi function to get constraints and objective
if S.options.useResidualforces
    coll_input_vars_def = {Qsk,Qsj,Qdotsk,Qdotsj,Aj,Qsj_nsc_exp,KCF_exp,knee_resj};
else
    coll_input_vars_def = {Qsk,Qsj,Qdotsk,Qdotsj,Aj,Qsj_nsc_exp,KCF_exp};
end
% if nq.torqAct > 0
%     coll_input_vars_def = [coll_input_vars_def,{a_ak,a_aj,e_ak}];
% end
f_coll = Function('f_coll',coll_input_vars_def,...
    {eq_constr,J});

% Repeat function for each mesh interval and assign evaluation to multiple threads
f_coll_map = f_coll.map(N,S.solver.parallel_mode,S.solver.N_threads);

% evaluate function with opti variables
if S.options.useResidualforces
    coll_input_vars_eval = {Qs(:,1:end-1), Qs_col, Qdots(:,1:end-1), Qdots_col,...
        A_col,S.expdata.IKdata.Qs_col',S.expdata.ContactForces_col', knee_res};
else
    coll_input_vars_eval = {Qs(:,1:end-1), Qs_col, Qdots(:,1:end-1), Qdots_col,...
            A_col,guess.Qs_col'.*scaling.Qs',S.expdata.ContactForces_col'};
end

[coll_eq_constr,Jall] = f_coll_map(coll_input_vars_eval{:});

% equality constraints
opti.subject_to(coll_eq_constr == 0);

% Loop over mesh points
if S.options.scale_kneety_av
    Qs_nsc=MX.zeros(size(Qs));
    Qs_nsc([1:4 6],:)=Qs([1:4 6],:).*scaling.Qs([1:4 6])';
    Qs_nsc(5,:)=Qs(5,:)*scaling.knee_ty.b+scaling.knee_ty.a;
    Qs_col_nsc=MX.zeros(size(Qs_col));
    Qs_col_nsc([1:4 6],:)=Qs_col([1:4 6],:).*scaling.Qs([1:4 6])';
    Qs_col_nsc(5,:)=Qs_col(5,:)*scaling.knee_ty.b+scaling.knee_ty.a;
end
for k=1:N
    % Variables within current mesh interval
    % States
    if S.options.scale_kneety_av
        Qskj = [Qs_nsc(:,k), Qs_col_nsc(:,(k-1)*d+1:k*d)];
    else
        Qskj = [Qs(:,k), Qs_col(:,(k-1)*d+1:k*d)];
    end
    Qdotskj = [Qdots(:,k), Qdots_col(:,(k-1)*d+1:k*d)];

    % Add equality constraints (next interval starts with end values of
    % states from previous interval)
    if S.options.scale_kneety_av
        opti.subject_to(Qs_nsc(:,k+1) == Qskj*D); % otherwise unscaled
    else
        opti.subject_to(Qs(:,k+1) == Qskj*D); % scaled 
    end
    opti.subject_to(Qdots(:,k+1) == Qdotskj*D); % scaled

end % End loop over mesh points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opti.minimize(sum(Jall));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['...OCP formulation done. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')

%%

if ~S.post_process.load_prev_opti_vars
    % Create NLP solver

    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.max_iter              = S.solver.max_iter;
    options.ipopt.linear_solver         = S.solver.linear_solver;
    options.ipopt.tol                   = 1*10^(-S.solver.tol_ipopt);
    options.ipopt.constr_viol_tol       = 1*10^(-S.solver.tol_ipopt);
    options.ipopt.warm_start_init_point = S.solver.warm_start_init_point;
    opti.solver('ipopt', options);
    % timer
    
    disp('Starting NLP solver...')
    disp(' ')
    t0s = tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve problem
    % Opti does not use bounds on variables but constraints. This function
    % adjusts for that.
    [w_opt,g_opt,stats,lam_x,lam_g,lam_p,llb,uub,new_g] = solve_NLPSOL(opti,options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ')
    disp(['...Exit NLP solver. Time elapsed ' num2str(toc(t0s),'%.2f') ' s'])
    disp(' ')
    disp(' ')
    % Create setup
    setup.tolerance.ipopt = S.solver.tol_ipopt;
    setup.bounds = bounds;
    setup.bounds_nsc = bounds_nsc;
    setup.scaling = scaling;
    setup.guess = guess;
    
    Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
    save(Outname,'w_opt','stats','setup','model_info','S');

else % S.post_process.load_prev_opti_vars = true
    
    % Advanced feature, for debugging only: load w_opt and reconstruct R before rerunning the post-processing.
    Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
    disp(['Loading vector with optimization variables from previous solution: ' Outname])
    clear 'S'
    load(Outname,'w_opt','stats','setup','model_info','R','S');
    scaling = setup.scaling;
    if exist('R','var')
        S = R.S;
    end
    clear 'R'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Essential post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read from the vector with optimization results
disp('Retrieving solution...')
disp(' ')


starti = 1;
Qs_opt = reshape(w_opt(starti:starti+nq*(N+1)-1),nq,N+1)';
starti = starti + nq*(N+1);
Qs_col_opt = reshape(w_opt(starti:starti+nq*(d*N)-1),nq,d*N)';
starti = starti + nq*(d*N);
Qdots_opt = reshape(w_opt(starti:starti+nq*(N+1)-1),nq,N+1)';
starti = starti + nq*(N+1);
Qdots_col_opt = reshape(w_opt(starti:starti+nq*(d*N)-1),nq,d*N)';
starti = starti + nq*(d*N);
qdotdot_col_opt =reshape(w_opt(starti:starti+nq*(d*N)-1),nq,(d*N))';
starti = starti + nq*(d*N);
if S.options.useResidualforces
    knee_res_opt = reshape(w_opt(starti:starti+2*(d*N)-1),2,(d*N))';
    starti = starti + 2*(d*N);
end
if starti - 1 ~= length(w_opt)
    disp('error when extracting results')
end

% Combine results at mesh and collocation points
Qs_mesh_col_opt=zeros(N*(d+1)+1,nq);
Qs_mesh_col_opt(1:(d+1):end,:)= Qs_opt;
Qdots_mesh_col_opt=zeros(N*(d+1)+1,nq);
Qdots_mesh_col_opt(1:(d+1):end,:)= Qdots_opt;

for k=1:N
    rangei = k*(d+1)-(d-1):k*(d+1);
    rangebi = (k-1)*d+1:k*d;
    Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
    Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
end

%% Unscale results
% States at mesh points
% Qs (1:N-1)
if S.options.scale_kneety_av
    q_opt_unsc.rad(:,[1:4 6]) = Qs_opt(1:end-1,[1:4 6]).*repmat(...
        scaling.Qs([1:4 6]),size(Qs_opt(1:end-1,[1:4 6]),1),1);
    q_opt_unsc.rad(:,5) = Qs_opt(1:end-1,5).*scaling.knee_ty.b+...
        scaling.knee_ty.a;
else
    q_opt_unsc.rad = Qs_opt(1:end-1,:).*repmat(...
        scaling.Qs,size(Qs_opt(1:end-1,:),1),1);
end
% Convert in degrees
q_opt_unsc.deg = q_opt_unsc.rad;
q_opt_unsc.deg(:,model_info.jointi.rotations) ...
    = q_opt_unsc.deg(:,model_info.jointi.rotations).*180/pi;
% Qs (1:N)
if S.options.scale_kneety_av
    q_opt_unsc_all.rad(:,[1:4 6]) = Qs_opt(:,[1:4 6]).*repmat(scaling.Qs([1:4 6]),size(Qs_opt(:,[1:4 6]),1),1);
    q_opt_unsc_all.rad(:,5) = Qs_opt(:,5).*scaling.knee_ty.b+scaling.knee_ty.a;
else
    q_opt_unsc_all.rad = Qs_opt.*repmat(scaling.Qs,size(Qs_opt,1),1);
end
% Convert in degrees
q_opt_unsc_all.deg = q_opt_unsc_all.rad;
q_opt_unsc_all.deg(:,model_info.jointi.rotations) ...
    = q_opt_unsc_all.deg(:,model_info.jointi.rotations).*180/pi;
% Qdots (1:N-1)
qdot_opt_unsc.rad = Qdots_opt(1:end-1,:).*repmat(...
    scaling.Qdots,size(Qdots_opt(1:end-1,:),1),1);
% Convert in degrees
qdot_opt_unsc.deg = qdot_opt_unsc.rad;
qdot_opt_unsc.deg(:,model_info.jointi.rotations) ...
    = qdot_opt_unsc.deg(:,model_info.jointi.rotations).*180/pi;
% Qdots (1:N)
qdot_opt_unsc_all.rad =Qdots_opt.*repmat(scaling.Qdots,size(Qdots_opt,1),1);

% States at collocation points
% Qs
if S.options.scale_kneety_av
    q_col_opt_unsc.rad(:,[1:4 6]) = Qs_col_opt(:,[1:4 6]).*repmat(scaling.Qs([1:4 6]),size(Qs_col_opt(:,[1:4 6]),1),1);
    q_col_opt_unsc.rad(:,5) = Qs_col_opt(:,5)*scaling.knee_ty.b+scaling.knee_ty.a;
else
    q_col_opt_unsc.rad = Qs_col_opt.*repmat(scaling.Qs,size(Qs_col_opt,1),1);
end
% Convert in degrees
q_col_opt_unsc.deg = q_col_opt_unsc.rad;
q_col_opt_unsc.deg(:,model_info.jointi.rotations) ...
    = q_col_opt_unsc.deg(:,model_info.jointi.rotations).*180/pi;
% Qdots
qdot_col_opt_unsc.rad = Qdots_col_opt.*repmat(...
    scaling.Qdots,size(Qdots_col_opt,1),1);
% Convert in degrees
qdot_col_opt_unsc.deg = qdot_col_opt_unsc.rad;
qdot_col_opt_unsc.deg(:,model_info.jointi.rotations) ...
    = qdot_col_opt_unsc.deg(:,model_info.jointi.rotations).*180/pi;
% "Slack" controls at collocation points
% Time derivative of Qdots
qdotdot_col_opt_unsc.rad = ...
    qdotdot_col_opt.*repmat(scaling.Qdotdots,size(qdotdot_col_opt,1),1);
% Convert in degrees
qdotdot_col_opt_unsc.deg = qdotdot_col_opt_unsc.rad;
qdotdot_col_opt_unsc.deg(:,model_info.jointi.rotations) ...
    = qdotdot_col_opt_unsc.deg(:,model_info.jointi.rotations).*180/pi;
if S.options.useResidualforces
    % residual forces
    knee_res_opt_unsc = knee_res_opt * scaling.residualforces;
end


%% Decompose optimal cost
J_opt           = 0;
Qdotdot_cost    = 0;
trackq_cost     = 0;
minsecq_cost    = 0;
trackKCF_cost    = 0;
trackresidualforce_cost = 0;
count           = 1;
h_opt           = (S.final_time-S.initial_time)/N;
for k=1:N
    for j=1:d
        % objective function
        J_opt = J_opt + ...
            W.q_dotdot*B(j+1)   *(f_casadi.uasqsum(qdotdot_col_opt(count,:)))*h_opt + ...
            W.trackq  *B(j+1)   *(sum((q_col_opt_unsc.rad(count,1)-S.expdata.IKdata.Qs_col(count,1)).^2))*h+...
            W.minsecq *B(j+1)   *(sum((q_col_opt_unsc.rad(count,[2 3 4 6])).^2))*h;

        %  W.trackq    * B(j+1)*(sum((Qskj_nsc(1,j+1)-Qsj_nsc_exp(1,j)).^2))*h+...
        % W.minsecq   * B(j+1)*(sum((Qskj_nsc([2 3 4 6],j+1)).^2))*h;
            % W.trackq  * B(j+1)  * ...
            % (f_casadi.ftrackkinexp(q_col_opt_unsc.rad(count,:), guess.Qs_col(count,:)'.*scaling.Qs'))*h_opt;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign Qs
        F_ext_input_opt = q_col_opt_unsc.rad(count,:);
           
        % Evaluate external function
        outi = F(F_ext_input_opt);
        KCFMed_opti=full(outi(1));
        KCFLat_opti=full(outi(2));
        KCFMed_opt(count)=KCFMed_opti;
        KCFLat_opt(count)=full(outi(2));
        %F is:
        % f_bothcompartments=Function('F',{x,y,z,psi,theta,phi},{KCFMed,KCFLat});
        % while the order of Qs is psi,theta,phi,x,y,z
        

        Qdotdot_cost = Qdotdot_cost + W.q_dotdot*B(j+1)*...
            (f_casadi.uasqsum(qdotdot_col_opt(count,:)))*h_opt;
        trackq_cost=trackq_cost+ W.trackq  *B(j+1)   *...
           (sum((q_col_opt_unsc.rad(count,1)-S.expdata.IKdata.Qs_col(count,1)).^2))*h;
        minsecq_cost= minsecq_cost + W.minsecq *B(j+1) *...
           (sum((q_col_opt_unsc.rad(count,[2 3 4 6])).^2))*h;
        % trackq_cost = trackq_cost + W.trackq  * B(j+1)  * ...
        %     (f_casadi.ftrackkinexp(q_col_opt_unsc.rad(count,:), ...
        %     guess.Qs_col(count,:)'.*scaling.Qs'))*h_opt;
        if strcmp(S.options.forcetrackingAs,'constraint')
            trackKCF_costf=0;
            if S.options.useResidualforces
                J_opt = J_opt+ W.minresidualforce * B(j+1) * (sum(knee_res_opt(count,:)).^2)*h;
                trackresidualforce_cost=trackresidualforce_cost+W.minresidualforce * B(j+1) * (sum(knee_res_opt(count,:)).^2)*h;
            else
                trackresidualforce_cost=0;
            end
            trackKCF_cost=0;
        else
            J_opt = J_opt+  W.trackKCF   * B(j+1) * (f_casadi.ftrackforceexp([KCFMed_opti; KCFLat_opti],S.expdata.ContactForces_col(count,:)))*h_opt;
            trackKCF_cost = trackKCF_cost + W.trackKCF   * B(j+1) * ...
                (f_casadi.ftrackforceexp([KCFMed_opti; KCFLat_opti],...
                S.expdata.ContactForces_col(count,:)))*h_opt;
            trackresidualforce_cost=0;
        end
        count = count + 1;
    end
end
J_optf = full(J_opt);
Qdotdot_costf = full(Qdotdot_cost);
trackq_costf = full(trackq_cost);
minsecq_costf= full(minsecq_cost);
trackKCF_costf = full(trackKCF_cost);
trackresidualforcef = full(trackresidualforce_cost);

contributionCost.absoluteValues = [Qdotdot_costf,trackq_costf,minsecq_costf,trackKCF_costf,trackresidualforcef];
contributionCost.relativeValues = [Qdotdot_costf,trackq_costf,minsecq_costf,trackKCF_costf,trackresidualforcef]...
    ./J_optf*100;
contributionCost.relativeValuesRound2 = ...
    round(contributionCost.relativeValues,2);
contributionCost.labels = {'joint accelerations','track qs','track KCF','trackresidual'};

% assertCost should be 0
assertCost = abs(J_optf - (Qdotdot_costf + trackq_costf + minsecq_costf + trackKCF_costf + trackresidualforcef));

assertCost2 = abs(stats.iterations.obj(end) - J_optf);

if assertCost > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing optimal cost wrt sum of terms')
    disp(['   Difference = ' num2str(assertCost)])
end
if assertCost2 > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing optimal cost wrt stats')
    disp(['   Difference = ' num2str(assertCost2)])
end

%% Save the results
% Structure Results_all
R.S = S;
R.objective = contributionCost;
R.time.mesh = tgrid;
R.time.coll = tgrid_ext;
R.colheaders.coordinates = model_info.coord_names.all;
R.colheaders.objective = contributionCost.labels;
R.kinematics.Qs = Qs_opt;
R.kinematics.Qdots = Qdots_opt;
R.kinematics.q_opt_unsc= q_opt_unsc.rad;
R.kinematics.qdot_opt_unsc = qdot_opt_unsc.rad;
R.kinematics.q_col_opt_unsc = q_col_opt_unsc.rad;
R.kinematics.qdot_col_opt_unsc = qdot_col_opt_unsc.rad;
R.kinematics.qdotdot_col_opt_unsc = qdotdot_col_opt_unsc.rad;
R.kinematics.q_opt_unsc_all.rad=[];
for i=1:N
    R.kinematics.q_opt_unsc_all.rad=[R.kinematics.q_opt_unsc_all.rad; q_opt_unsc.rad(i,:); q_col_opt_unsc.rad((i-1)*d+1:i*d,:)];
end
R.kinematics.q_opt_unsc_all.rad=[R.kinematics.q_opt_unsc_all.rad; q_opt_unsc_all.rad(N+1,:)];
R.ContactForces=[KCFMed_opt' KCFLat_opt'];

% save results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
disp(['Saving results as: ' Outname])
save(Outname,'w_opt','stats','setup','R','model_info');

end