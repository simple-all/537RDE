function [vcj,pr,T2,gamma2]=CJ_detonation(f,T1)
error=1;

h=120e6; % Heating value of the hydrogen 120e6 J/kg
R=287;   % gas constant for air [J/(kg K)]

cp1=combustor.cp_air(T1);
gamma1=cp1/(cp1-R);

cp2=cp1;
gamma2=gamma1;

while error>1e-10
    Mcj = combustor.CJ_mach(h,f,T1,gamma1,gamma2,cp1,cp2);
    Mcj=real(Mcj);
    T01=T1*(1+0.5*(gamma1-1)*(Mcj^2));
    T02=(1/cp2)*(f*h+T01*cp1);
    T2=T02/(0.5*(gamma2+1));
    
    cp2new=combustor.cp_air(T2);
    error=(cp2-cp2new)/cp2;
    
    cp2=(1+0.95*error)*cp2new;
    gamma2=cp2/(cp2-R);
    
    error=abs(error);
end

[D1,~]=combustor.D_dD(Mcj,gamma1);
[D2,~]=combustor.D_dD(1,gamma2);

p02_p01=sqrt((T02*gamma1)/(T01*gamma2))*(D1/D2);
pr=p02_p01*((1+0.5*(gamma1-1)*(Mcj^2))^(gamma1/(gamma1-1)))/((0.5*(gamma2+1))^(gamma2/(gamma2-1)));

vcj=Mcj*sqrt(gamma1*R*T1);
end