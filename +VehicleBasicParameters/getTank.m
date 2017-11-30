%% Calculates the weight and volume of a pill style tank given...
% 1. Normal Operating Pressure (P_operating) [psi]
% 2. Diameter of the tank (d) [inches]
% 3. Overall length of tank (l)
% 4. Fluid Density (fluid_density) [lb/ft3]

function [Weight_tank, Volume_tank, Weight_fluid, Surface_Area] = getTank(P_operating, d, l, fluid_density)
r = d/2; %[in]
a = l - (2*d); %[in]
Surface_Area = 2*pi*r*(2*r + a); %[sq.in]
Volume_tank = pi*r^2*((4/3)*r+a); %[cu.in]
Weight_tank = (3*P_operating*Volume_tank)/(1*10^6); %[lbf] %Burst pressure is 3 times max normal operating
Weight_fluid = Volume_tank * (fluid_density/1728); %[lbf]
end