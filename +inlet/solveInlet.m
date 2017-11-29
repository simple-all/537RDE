function [M2, P2, T2, mdot, coneLength, q, ramDrag] = solveInlet(inletDiameter, coneAngle, M0, Pr, altitude, mDrop)
%SOLVEINLET Summary of this function goes here
%   Detailed explanation goes here
gamma = 1.4; % Ratio of specific heats
R_air = 287.058; % [J/kg*K] Air gas constant

[T0, P0, rho0] = atmosphere.atmosphere_metric(altitude, 1);
Pt0 = aeroBox.isoBox.calcStagPressure('mach', M0, 'gamma', gamma, 'Ps', P0);
Tt = aeroBox.isoBox.calcStagTemp('mach', M0, 'gamma', gamma, 'Ts', T0);

a0 = sqrt(gamma * R_air * T0); % [m/s] Free stream sound speed
u0 = a0 * M0; % [m/s] Free stream velocity

q = 0.5 * gamma * P0 * M0^2;

inletArea = pi * (inletDiameter / 2)^2; % [m^2]
mdot = u0 * inletArea * rho0; % [kg/s] Mass flow through inlet
% After shock properties

% calculates the conical shock angle
shockAngle=conical.find_cone_shock_angle(M0,coneAngle,gamma);
% returns the flow field solution for the given cone angle
% [v,mn1]=conical.taymacsol2(M0,shockAngle,gamma);
% % output and answer handling script
% % resolves the radial and angular velocity components into a % single ray velocity value
% vdash=sqrt((v(:,1).^2)+(v(:,2).^2));
% % converts the velocity values into Mach numbers
% mach=sqrt(2./(((vdash.^(-2))-1).*(gamma-1)));
% % find the Mach number of the ray nearest the cone
% M1=mach(length(mach));  
% % calculates the Temperature ratio for the ray nearest the cone
% To_T=1+(((gamma-1)/2).*M1.^2);
% % calculates the Pressure ratio for the ray nearest the cone
% Po_P=(1+(((gamma-1)/2).*M1.^2)).^(gamma/(gamma-1));
% % calculates the Density ratio for the ray nearest the cone
% rho_rh=(1+(((gamma-1)/2).*M1.^2)).^(1/(gamma-1));
% % Calculated the Stagnation pressure ratio (after/before conical shock)
% Po2_Po1=((((gamma+1)/2.*mn1.^2)./(1+(((gamma-1)/2).*mn1.^2))).^(gamma/(gamma-1)))./...
%     ((((2*gamma/(gamma+1)).*mn1.^2)-((gamma-1)/(gamma+1))).^(1/(gamma-1)));
% 
% % Get density behind the shock
% turnAngle = shockAngle;
% rho1 = rho0 * (((gamma + 1) * M0^2 * sind(turnAngle)^2) / ((gamma - 1) * M0^2 * sind(turnAngle)^2 + 2));
% T1 = T0 * (((2 * gamma * M0^2 * sind(turnAngle)^2 - (gamma - 1)) * ((gamma - 1) * M0^2 * sind(turnAngle)^2 + 2)) / ((gamma + 1)^2 * M0^2 * sind(turnAngle)^2));
% Tt = aeroBox.isoBox.calcStagTemp('mach', M1, 'gamma', gamma, 'Ts', T1);
% a1 = sqrt(gamma * R_air * T1);
% u1 = M1 * a1;

coneLength = (inletDiameter / 2) / tand(shockAngle);
ramDrag = mdot * u0; % Ram drag
M2 = mDrop * M0;

% Isolator
Pt2 = Pr * Pt0;
P2 = aeroBox.isoBox.calcStaticPressure('mach', M2, 'Pt', Pt2, 'gamma', gamma);
T2 = aeroBox.isoBox.calcStaticTemp('mach', M2, 'Tt', Tt, 'gamma', gamma);
end

