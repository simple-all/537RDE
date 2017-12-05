clear;
close all;
clc;

% Start at altitude of 23.317 km at 1634.956 m/s

craft = sim.Vehicle(1634.95646331355, 23317);

% Aero params
craft.setAeroParams(0.03, 0.015, 0.929); % Cl, Cd, Area

% Physical params
craft.setMass(296.3, 0.0708); % Wet mass [kg], fuel mass fraction

% Inlet params
craft.setInlet(0.24, 0.127, 10, 0.3, (1/2)); % Inlet diameter [m], inlet gap [m] cone half angle [deg], Pr, Mr

% RDE params
craft.setRDE(0.2035, 0.16);

% Flight profile params
craft.setPhis(1, 0.58); % Phi accel, phi cruise
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

