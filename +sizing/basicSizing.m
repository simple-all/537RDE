M0 = 6; % Free stream Mach
q = 1500 * 47.8802589; % [Pa]
mdot = 2.85; % [kg/s] Air mass flow rate
M2 = 2.3; % Isolator exit mach
P2 = 101325 * 2; % [Pa] Isolator exit static pressure
coneAngle = 10; % [deg] Inlet cone half angle


% Inlet sizing
[inletDiameter, inletGap, inletSystemLength, T2, Pr, Pr_is, altitude, P0] = inlet.genInlet(M0, q, mdot, M2, P2, coneAngle);
if (Pr_is >= 1)
    warning('Invalid design! Isolator recovery pressure too high');
end
fprintf('Operating at %0.3f km and a pressure recovery of %0.3f\n', altitude / 1e3, Pr);







