clear;
clc;
% flow_propertycalc.m;
% Taylor-Maccoll conical flow solution using equations given by
% Anderson, 1990, "Modern Compressible Flow (with Historical % Perspective)"
% 2nd edition, McGraw-Hill, p301-303. %
% Script co-ordinates the overall solution. Outputs are the
% isentropic property ratios and some static/stagnation values for the
% flow field at the cone angle %
%
% David BT Sercombe, 3 January 2004 % David R Buttsworth,
% Preset values for free stream Mach no. (M) Cone Angle (thetac) and
% Ratio of Specific Heats (g) are specific to the experimental cone$ in use.
M=7;
thetac=10;
g=1.4;
% calculates the conical shock angle 
thetas=conical.find_cone_shock_angle(M,thetac,g);
% returns the flow field solution for the given cone angle 
[v,mn1]=conical.taymacsol2(M,thetas,g);
% output and answer handling script
% resolves the radial and angular velocity components into a % single ray velocity value 
vdash=sqrt((v(:,1).^2)+(v(:,2).^2));
% converts the velocity values into Mach numbers 
mach=sqrt(2./(((vdash.^(-2))-1).*(g-1)));
% find the Mach number of the ray nearest the cone 
mfin=mach(length(mach));
% calculates the Temperature ratio for the ray nearest the cone 
To_T=1+(((g-1)/2).*mfin.^2);
% calculates the Pressure ratio for the ray nearest the cone 
Po_P=(1+(((g-1)/2).*mfin.^2)).^(g/(g-1));
% calculates the Density ratio for the ray nearest the cone 
rho_rh=(1+(((g-1)/2).*mfin.^2)).^(1/(g-1));
% Calculated the Stagnation pressure ratio (after/before conical shock) 
Po2_Po1=((((g+1)/2.*mn1.^2)./(1+(((g-1)/2).*mn1.^2))).^(g/(g-1)))./...
((((2*g/(g+1)).*mn1.^2)-((g-1)/(g+1))).^(1/(g-1)));
% asks user for the free stream stagnation pressure
%Po1=input('What is the Free Stream Staganation Pressure (P_o1)in MPa?:'); % calculates the post shock stagnation pressure
%Po2=Po2_Po1*Po1;
% answer printing scripts
fprintf('\nFor a Cone Angle of %.0f degrees and a Free Stream Mach no. ... of %.0f \nThe Shock Wave Angle is %.2f degrees \nThe Temperature Ratio ... (To_2/T) is %.4f at the Cone Angle \nThe Pressure Ratio (Po_2/P) is %.4f... \nThe Density Ratio (po_2/p) is %.4f\n',thetac,M,thetas,To_T,Po_P,rho_rh);
%fprintf('\nThe Stagnation Pressure Ratio (After/Before Shock) (Po_2/Po_1) is %.4f \nHence the Stagnation Pressure behind the shockwave MPa',Po2_Po1,Po2)
% EOF
