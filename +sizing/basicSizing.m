clear

M0 = 7; % Free stream Mach
q = 1500 * 47.8802589; % [Pa]
mdot_air = 3; % [kg/s] Air mass flow rate
M2 = 2; % Isolator exit mach
P2 = 101325 * 2; % [Pa] Isolator exit static pressure
coneAngle = 10; % [deg] Inlet cone half angle


% Inlet sizing
[inletDiameter, inletGap, inletSystemLength, T2, Pr, Pr_is, altitude, P0] = inlet.genInlet(M0, q, mdot_air, M2, P2, coneAngle);
if (Pr_is >= 1)
    warning('Invalid design! Isolator recovery pressure too high');
end
fprintf('Operating at %0.3f km and a pressure recovery of %0.3f\n', altitude / 1e3, Pr);


% [Pa] Atmpsheric pressure
altPa = 1e3 * [101.33 99.49 97.63 95.91 94.19 92.46 90.81 89.15 87.49 85.91 85.44 81.22 78.19 75.22 72.40 69.64 57.16 46.61 37.65 30.13 23.93 18.82 14.82 11.65 9.17 7.24 4.49 2.80 1.76 1.12 0.146 2.2e-2 1.09e-6 4.98e-7 4.8e-10];
% [m] Altitude
alth = 0.3048 * [0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 6000 7000 8000 9000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000 70000 80000 90000 100000 150000 200000 300000 500000 2000000];
% [K] Atmospheric temperature
altT = 273 + [15 14 13 12 11 10 9 8 7 6 5 3 1 -1 -3 -5 -14 -24 -34 -44 -54 -57 -57 -57 -57 -57 -55 -52 -59 -46 -46 -46 -46 -46 -46];

fprintf('Altitude: %0.3f km\n', altitude / 1e3);
T0 = interp1(alth, altT, altitude, 'linear');

areaInlet = pi * (inletDiameter / 2)^2;
aInlet = sqrt(1.4 * 287 * T0);
uInlet = aInlet * M0;
ramDrag = mdot_air * uInlet;

% Combustor
D_outer = 11 * 0.0254;
D_inner = 9 * 0.0254;

Pmin = 20e6; % [Pa] Initial guess at minimum pressure
FAR = 0.029; % Fuel to Air ratio

numDets = 1;

% From CEA
Pr = 6.909; % Pressure ratio
Tmax = 3510.36; % [K] Temperature of burned gas
v_cj = 2053.2; % [m/s] CJ Velocity
R = 8314/24.373; % [J/kg*K] Gas constant of burned gas
gamma = 1.2004; % Ratio of specific heats
tsteps = 1000;
P0 = 3e3; % Atmospheric air pressure


Pr = 5.19; % Pressure ratio
Tmax = 3542.82; % [K] Temperature of burned gas
v_cj = 2023.2; % [m/s] CJ Velocity
R = 8314/24.107; % [J/kg*K] Gas constant of burned gas
gamma = 1.1854; % Ratio of specific heats
tsteps = 1000;
P0 = 3e3; % Atmospheric air pressure
M = 2; % Combustion mach

[Isp, F, Pmin] = combustor.solveRDE(Pr, Pmin, Tmax, v_cj, R, D_outer, D_inner, mdot_air, FAR, gamma, tsteps, P0, numDets, M);

fprintf('Isp: %0.2f s\nThrust: %0.3f kN\nMinimum chamber pressure: %0.3f kPa\n', Isp, F / 1000, Pmin / 1000);
fprintf('Total Thrust: %0.3f kN\n', (F - ramDrag) / 1000);




