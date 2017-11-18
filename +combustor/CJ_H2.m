clear all
close all
clc

error=1;

h=120e6; % Heating value of the hydrogen 120e6 J/kg
f=0.029; % fuel to air ratio [-]
T1=750;  % burner inlet temperature [K]
R=287;   % gas constant for air [J/(kg K)]

cp1=cp_air(T1);
gamma1=cp1/(cp1-R);

cp2=cp1;
gamma2=gamma1;

while error>1e-10
    Mcj = CJ_mach(h,f,T1,gamma1,gamma2,cp1,cp2);
    Mcj=real(Mcj);
    T01=T1*(1+0.5*(gamma1-1)*(Mcj^2));
    T02=(1/cp2)*(f*h+T01*cp1);
    T2=T02/(0.5*(gamma2+1));
    
    cp2new=cp_air(T2);
    error=(cp2-cp2new)/cp2;
    
    cp2=(1+0.95*error)*cp2new;
    gamma2=cp2/(cp2-R);
    
    fprintf('T2 = %.0f         Mcj = %.3f\n',T2,Mcj)
    error=abs(error);
end

[D1,~]=D_dD(Mcj,gamma1);
[D2,~]=D_dD(1,gamma2);

p02_p01=sqrt((T02*gamma1)/(T01*gamma2))*(D1/D2);
pr=((1+0.5*(gamma1-1)*(Mcj^2))^(gamma1/(gamma1-1)))/((0.5*(gamma2+1))^(gamma2/(gamma2-1)));
fprintf('pr = %.3f',pr)