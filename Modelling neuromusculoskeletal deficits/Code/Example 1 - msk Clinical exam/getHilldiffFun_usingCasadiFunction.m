function Hilldiff_fun = getHilldiffFun(fse,a,dfse,lMT,vMT,FMo_in,lMo_in,...
    lTs_in,alphao_in,vMmax_in,Fvparam,Fpparam,Faparam,tension,aTendon,shift,...
    MuscMoAsmp,d,stiffness_shift,stiffness_scale,strength)

N_fun = length(lMT);
Hilldiff_fun = zeros(size(fse));

for m = 1:N_fun
    [Hilldiff_tmp,~,~,~,~,~,~] = ...
        ForceEquilibrium_FtildeState_all_tendon(a(m),fse(m),dfse,lMT(m),vMT,...
        FMo_in(m),lMo_in(m),lTs_in(m),alphao_in(m),vMmax_in(m),...
        Fvparam,Fpparam,Faparam,...
        tension(m),aTendon(m),shift(m),MuscMoAsmp,d,stiffness_shift(m),stiffness_scale(m),strength(m));

    Hilldiff_fun(m) = Hilldiff_tmp;
end

end