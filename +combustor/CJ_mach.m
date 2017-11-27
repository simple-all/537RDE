function Mcj = CJ_mach(h,f,T1,gamma1,gamma2,cp1,cp2)
F=@fun;
xo=2;
Mcj=fsolve(F,xo,optimset('LargeScale','off','Display','off','TolFun',1e-20,'MaxIter',1e5,'MaxFunEvals',1e5));
function F = fun(Mcj)
    T01=T1*(1+0.5*(gamma1-1)*(Mcj^2));
    F=sqrt((gamma1/gamma2)*(((cp1/cp2)*T01+f*h/cp2)/(T01)))*N_fun(gamma1,Mcj)-1/(sqrt(2*(gamma2+1)));
end

function N = N_fun(gamma,M)
N=M*((1+0.5*(gamma-1)*(M^2))^0.5)/(1+gamma*(M^2));
end

end