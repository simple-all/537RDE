clear;
close all;
clc;

% Start at altitude of 25.963 km at 1,625.6 m/s

craft = sim.Vehicle(1425.6, 23000);

% Aero params
craft.setAeroParams(0.03, 0.0100775, 0.929); % Cl, Cd, Area

% Physical params
craft.setMass(272, 0.1); % Wet mass [kg], fuel mass fraction

% Inlet params
craft.setInlet(0.2, 10, 0.4, 1/3); % Inlet diameter [m], cone half angle [deg], Pr, Mr

% RDE params
craft.setRDE(0.2, 0.17);


dt = 0.5;
i = 1;
while(craft.fuelMass > 0)
    craft.runAtTime(craft.time + dt);
    time(i) = craft.time;
    M = craft.getMach();
    fm = craft.fuelMass;
    q = craft.getDynamicPressure();
    fprintf('Time: %0.1f, Mach: %0.3f, Q: %0.3f, Alt: %0.3f km, Fuel: %0.3f kg\n', time(i), M, q, craft.altitude / 1e3, fm);
    i = i + 1;
end

