classdef Vehicle < handle
    %@VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wetMass % [kg] Initial wet mass of vehicle
        fuelMassFraction;
    end
    
    properties (SetAccess = private)
        % Sim parameters
        time = 0; % [s] Current time in simulation/flight
        lastRunTime = 0;
        
        % Movement parameters
        velocity; % [m/s] Flight velocity
        altitude; % [m] Flight altitude
        
        % Engine parameters
        detOuterDiameter; % [m] Outer diameter of engine
        detInnerDiameter; % [m] Inner diameter of engine
            
        % Physical parameters
        dryMass; % [kg] Dry mass of vehicle
        fuelMass; % [kg] Mass of fuel
        
        % Aero parameters
        Cl; % Lift coefficient
        Cd; % Drag coefficient
        liftingArea % [m^2] Lifting body area
    end
    
    properties (Constant)
        % Atmospheric parameters
        gamma_air = 1.4; % Ratio of specific heats
        R_air = 287.058; % [J/kg*K] Air gas constant
        % [Pa] Atmpsheric pressure
        altPa = 1e3 * [101.33 99.49 97.63 95.91 94.19 92.46 90.81 89.15 87.49 85.91 85.44 81.22 78.19 75.22 72.40 69.64 57.16 46.61 37.65 30.13 23.93 18.82 14.82 11.65 9.17 7.24 4.49 2.80 1.76 1.12 0.146 2.2e-2 1.09e-6 4.98e-7 4.8e-10];
        % [m] Altitude
        alth = 0.3048 * [0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 6000 7000 8000 9000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000 70000 80000 90000 100000 150000 200000 300000 500000 2000000];
        % [K] Atmospheric temperature
        altT = 273 + [15 14 13 12 11 10 9 8 7 6 5 3 1 -1 -3 -5 -14 -24 -34 -44 -54 -57 -57 -57 -57 -57 -55 -52 -59 -46 -46 -46 -46 -46 -46];
    end
    
    methods
        function obj = Vehicle()
            % Create a vehicle
        end
        
        function runAtTime(obj, time)
            obj.time = time;
            dt = obj.time - obj.lastRunTime;
            
            
            
            obj.lastRunTime = obj.time;
        end
        
        function 
        
        function setAeroParams(obj, Cl, Cd, liftingArea)
            obj.Cl = Cl;
            obj.Cd = Cd;
            obj.liftingArea = liftingArea;
        end
        
        function setMass(obj, wetMass, fuelMassFraction)
            obj.wetMass = wetMass;
            obj.fuelMassFraction = fuelMassFraction;
            
            obj.fuelMass = obj.wetMass * obj.fuelMassFraction;
            obj.dryMass = obj.wetMass - obj.fuelMass;
        end
        
        function lift = obj.getLift(obj)
            lift = obj.Cl * obj.getDynamicPressure() * obj.liftingArea;
        end
        
        function drag = obj.getDrag(obj)
            drag = obj.Cd * obj.getDynamicPressure() * obj.liftingArea;
        end
        
        function mass = getTotalMass(obj)
            mass = obj.dryMass + obj.fuelMass;
        end
        
        function q = getDynamicPressure(obj)
            P0 = obj.getLocalFlowParams();
            q = 0.5 * obj.gamma_air * P0 * obj.getMach()^2;
        end
        
        function M = getMach(obj)
            M = obj.veloctiy / obj.getSonicVelocity();
        end
        
        function a = getSonicVelocity(obj)
            a = sqrt(obj.gamma_air * obj.R_air * T0);
        end
        
        function [P0, T0] = getLocalFlowParams(obj)
            P0 = interp1(obj.alth, obj.altPa, obj.altitude, 'linear'); % Freestream static pressure
            T0 = interp1(obj.alth, obj.altT, obj.altitude, 'linear'); % Freestream static temperature
        end
    end
    
end

