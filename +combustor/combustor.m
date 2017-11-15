clear all
close all
clc

n=1e3;

m3=8; % inlet
gamma=1.4; % inlet
M3=2.5; % inlet
p3=2*101325; % inputs CEA
T3=1000; % inputs CEA
pt3=p3*(1+0.5*(gamma-1)*(M3^2))^(gamma/(gamma-1)); % inlet

Tmax=3000; % CEA
pr=6; % CEA
gamma=1.3; % CEA
R=287; % CEA

pmin=p3;
pmax=pr*p3;

theta=linspace(0,2*pi,n);
pc=zeros(1,n);

for i=1:n
    pc(i) = p_det(theta(i),pr,pmin);
end


Tc=zeros(1,n);
for i=1:n
    Tc(i)=T_det (Tmax,pc(i),pmax,gamma);
end

figure
plot(theta,pc)

figure
plot(theta,Tc)

Rmax=280; % Geometry
Rmin=200; % Geometry
f =0.029; % fuel to air ratio
m4=m3*(1+f);
vx=m4/(0.5*(((Rmax)^2)-((Rmin)^2))*K(gamma,R,pmax,pr,Tmax));