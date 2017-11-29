clear;
close all;
clc;

% Start at altitude of 23.317 km at 1634.956 m/s

craft = sim.Vehicle(1634.95646331355, 23317);

% Aero params
craft.setAeroParams(0.03, 0.0100775, 0.929); % Cl, Cd, Area

% Physical params
craft.setMass(272, 0.15); % Wet mass [kg], fuel mass fraction

% Inlet params
craft.setInlet(0.24, 10, 0.325, (1/3)); % Inlet diameter [m], cone half angle [deg], Pr, Mr

% RDE params
craft.setRDE(0.1735, 0.13);

% Flight profile params
craft.setPhis(1, 0.515); % Phi accel, phi cruise
craft.setTarget(6.5, 1500 * 47.8802589); % Target Mach, target dynamic pressure

dt = 0.5;
i = 1;
while(craft.fuelMass > 0)
    craft.runAtTime(craft.time + dt);
    time(i) = craft.time;
    M = craft.getMach();
    fm = craft.fuelMass;
    q = craft.getDynamicPressure();
    fprintf('Time: %0.1f, Mach: %0.3f, q: %0.3f, Alt: %0.3f km, Fuel: %0.3f kg\n', time(i), M, q, craft.altitude / 1e3, fm);
    i = i + 1;
end

