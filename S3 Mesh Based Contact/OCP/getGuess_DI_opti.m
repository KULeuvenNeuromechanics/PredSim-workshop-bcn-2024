function guess = getGuess_DI_opti(S,model_info,scaling,d,F)
% --------------------------------------------------------------------------
% getGuess_DI_opti
%   This script provides an inital guess for the design variables.
%   The guess is data-informed (DI). We use experimental data to provide an
%   initial guess of the joint variables but set constant values to the 
%   muscle variable and the arm variables. We use a pre-defined final time 
%   that is function of the imposed speed.
%   
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - scaling -
%   * scale factors for all optimisation variables
% 
%   - d -
%   * degree of the interpolating polynomial of the collocation scheme
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: Gil Serrancolí
% Last edit date: 24/09/2024
% --------------------------------------------------------------------------

N = S.solver.N_meshes; % number of mesh intervals
nq = model_info.nq;
coordinate_names = model_info.coord_names.all;

%%
guess.Qs=S.expdata.IKdata.Qs;
if S.options.preoptimizekneedofs==1
    fprintf('Pre-optimizing knee sup-inf translation... \n');
    guess.Qs(:,5)=OptimizeForkneety(guess.Qs,S.expdata.ContactForces,N,d,F);
elseif S.options.preoptimizekneedofs==2
    fprintf('Pre-optimizing knee sup-inf translation and knee abd-add angle... \n');
    guess.Qs(:,[2 5])=OptimizeForkneeaddty(guess.Qs,S.expdata.ContactForces,N,d,F);
end
% guess.Qs(:,5)=0.043;
guess.Qdots=S.expdata.IKdata.Qdots;
guess.Qdotdots=S.expdata.IKdata.Qdotdots;

%% Collocation points
guess.Qs_col = zeros(d*N,nq);
guess.Qdots_col = zeros(d*N,nq);
guess.Qdotdots_col = zeros(d*N,nq);
guess.Qs_col=S.expdata.IKdata.Qs_col;
if S.options.preoptimizekneedofs==1
    guess.Qs_col(:,5)=OptimizeForkneety(guess.Qs_col,S.expdata.ContactForces_col,N,d,F);
    fprintf('Knee sup-inf translation pre-optimized... \n');
elseif S.options.preoptimizekneedofs==2
    guess.Qs_col(:,[2 5])=OptimizeForkneeaddty(guess.Qs_col,S.expdata.ContactForces_col,N,d,F);
    fprintf('Knee sup-inf translation and knee abd-add angle pre-optimized... \n');
end
guess.Qdots_col = S.expdata.IKdata.Qdots_col;
guess.Qdotdots_col = S.expdata.IKdata.Qdotdots_col;

guess.residualforces=zeros(N*d,2);


%% Scaling
guess.Qs = guess.Qs./repmat(scaling.Qs,N+1,1);
if S.options.scale_kneety_av
    guess.Qs(:,5)=(guess.Qs(:,5) - scaling.knee_ty.a)/scaling.knee_ty.b;
end
guess.Qdots = guess.Qdots./repmat(scaling.Qdots,N+1,1);
guess.Qdotdots = guess.Qdotdots(1:N,:)./repmat(scaling.Qdotdots,N,1);

guess.Qs_col = guess.Qs_col./repmat(scaling.Qs,N*d,1);
if S.options.scale_kneety_av
    guess.Qs_col(:,5)=(guess.Qs_col(:,5)-scaling.knee_ty.a)/scaling.knee_ty.b;
end
guess.Qdots_col = guess.Qdots_col./repmat(scaling.Qdots,N*d,1);
guess.Qdotdots_col = guess.Qdotdots_col(1:N*d,:)./repmat(scaling.Qdotdots,N*d,1);

end

function knee_ty=OptimizeForkneety(guess_Qs,KCF_exp,N,d,F);

nq=size(guess_Qs,2);
for i=1:size(guess_Qs,1)
    % X=[];
    % for j=1:nq
    %     X=[X guess_Qs(i,j)' 0];
    % end
    outi=F(guess_Qs(i,:));
    KCF_Mi=outi(1);
    KCF_Li=outi(2);

    
    q0=0.042;
    options=[];
    options.Display='none';
    qopt = lsqnonlin(@(x)fun(x,guess_Qs(i,:),F,KCF_exp(i,:),nq),q0,0.03,0.055,options);
    newqs(i,:)=guess_Qs(i,:);
    newqs(i,5)=qopt;
    outi= F(newqs(i,:));
    KCF_Mi=outi(1);
    KCF_Li=outi(2);
    KCF_M(i,:)=full(KCF_Mi);
    KCF_L(i,:)=full(KCF_Li);
    f=fun(qopt,newqs(i,:),F,KCF_exp(i,:),nq);
knee_ty(i)=qopt;    
    
end

end

function out=OptimizeForkneeaddty(guess_Qs,KCF_exp,N,d,F);

nq=size(guess_Qs,2);
for i=1:size(guess_Qs,1)
    % X=[];
    % for j=1:nq
    %     X=[X guess_Qs(i,j)' 0];
    % end
    outi=F(guess_Qs(i,:));
    KCF_Mi=full(outi(1));
    KCF_Li=full(outi(2));
    
    % if i==1
        q0_unsc=[0 0.042];
    % else
    %     q0_unsc=newqs(i-1,[2 5]);
    %     q0_unsc(2)=q0_unsc(2)-1e-3;
    % end
    q0(1)=q0_unsc(1)/0.1;
    q0(2)=(q0_unsc(2)-0.042)/0.007;
    options=[];
    options.Display='none';
    qopt = lsqnonlin(@(x)fun2(x,guess_Qs(i,:),F,KCF_exp(i,:),nq),q0,[-1 -1],[1 1],options);
    newqs(i,:)=guess_Qs(i,:);
    newqs(i,2)=qopt(1)*0.1;
    newqs(i,5)=qopt(2)*0.007+0.042;
    out= F(newqs(i,:));
    KCF_Mi=out(1);
    KCF_Li=out(2);
    KCF_M(i,:)=full(KCF_Mi);
    KCF_L(i,:)=full(KCF_Li);
    f=fun2(qopt,newqs(i,:),F,KCF_exp(i,:),nq);
    knee_add(i)=newqs(i,2);
    knee_ty(i)=newqs(i,5);    
    
end
knee_add=knee_add';
knee_ty=knee_ty';
out=[knee_add knee_ty];
end

function f=fun(x,Qs,F,KCFexp,nq)
    Qs(5)=x;
    out= F(Qs);
    KCF_M=full(out(1));
    KCF_L=full(out(2));
    
    f=sum(((KCFexp-[KCF_M KCF_L])/100).^2);
end

function f=fun2(x,Qs,F,KCFexp,nq)
    Qs(2)=x(1)*0.1;
    Qs(5)=x(2)*0.007+0.042;
    out= F(Qs);
    KCF_M=full(out(1));
    KCF_L=full(out(2));
    
    f=(((KCFexp-[KCF_M KCF_L])/1000));
end


