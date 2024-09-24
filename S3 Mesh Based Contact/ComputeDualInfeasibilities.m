% Function to call the dual infeasibility
% this assumes that you are within the OCP_formulation, and you have called
% opti and solved the problem, ending or not to an optimal solution
%
% Author: Gil Serrancolí
% Last edit: September 2024 


lamg=lam_g;
lamx=lam_x;
lamp=lam_p;

x_val=w_opt;
gradf=Function('gradf',{opti.x},{jacobian(opti.f,opti.x)});
gradfval=full(gradf(x_val));

eqcons_pos=find(llb==uub);
% ineqcons_pos=find(full(opti.debug.value(opti.lbg))~=full(opti.debug.value(opti.ubg)));
ineqcons_pos=find(llb~=uub);
gg=new_g;
gg_eq=gg(eqcons_pos);
jacob_eqf=Function('jacob_eqf',{opti.x},{jacobian(gg_eq,opti.x)});
jacob_eq_val=full(jacob_eqf(x_val));


dual_infeasibilities=gradfval'+jacob_eq_val'*lamg(eqcons_pos)+lamx; %Eq 4a in Wachter and Biegler 2006 has negative sign for last term. However, this expression (with positive sign) is consistent with CasADi/Ipopt results.


% dual_infeasibilities=gradfval'+jacob_eq_val'*lamg(eqcons_pos)+lamg(ineqcons_pos); %Eq 4a in Wachter and Biegler 2006 has negative sign for last term. However, this expression (with positive sign) is consistent with CasADi/Ipopt results.
