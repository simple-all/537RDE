D_outer = 11 * 0.0254;
D_inner = 9 * 0.0254;

Pmin = 20e6; % [Pa] Initial guess at minimum pressure
mdot_air = 8; % kg/s
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

[Isp, F, Pmin] = combustor.solveRDE(Pr, Pmin, Tmax, v_cj, R, D_outer, D_inner, mdot_air, FAR, gamma, tsteps, P0, numDets);
