classdef Vehicle < handle
    %@VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        wetMass % [kg] Initial wet mass of vehicle
        fuelMassFraction;
        phi_accel = 1.4;
        phi_cruise = 0.5;
        targetMach = 6.5;
        targetQ = 1500 * 47.8802589;
        
        % Sim parameters
        time = 0; % [s] Current time in simulation/flight
        lastRunTime = 0;
        
        % Movement parameters
        velocity; % [m/s] Flight velocity
        altitude; % [m] Flight altitude
        
        % Inlet parameters
        inletDiameter; % [m] Outer diameter of inlet
        inletConeAngle; % [deg] Half-angle of inlet cone
        Pr_inlet; % Pressure recovery factor of the inlet
        Mr_inlet = 0.33; % Mach ratio of the inlet-isolator system
        
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
        
        % Log parameters
        log;
        
        % CEA
        cea;
    end
    
    properties (Access = private)
        lastPmin = 118210; 
        lastDP = 1000;
        qLatch = 0;
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
        function obj = Vehicle(velocity, altitude)
            % Create a vehicle
            
            % Initiate log parameters
            obj.log = struct();
            obj.log.i = 1;
            
            % Create CEA runner
            obj.cea = nasa.CEARunner;
            
            obj.velocity = velocity;
            obj.altitude = altitude;
        end
        
        function runAtTime(obj, time)
            obj.time = time;
            dt = obj.time - obj.lastRunTime;
            
            [netThrust, Isp, Pmin, coneLength, q, P2, mdot_air] = obj.solveFlowpath(dt);
            
            i = obj.log.i;
            obj.log.netThrust(i) = netThrust;
            obj.log.Isp(i) = Isp;
            obj.log.Pmin(i) = Pmin;
            obj.log.coneLength(i) = coneLength;
            obj.log.q(i) = q;
            obj.log.P2(i) = P2;
            obj.log.mdot_air(i) = mdot_air;
            obj.log.M(i) = obj.getMach();
            obj.log.altitude(i) = obj.altitude;
            
            obj.lastRunTime = obj.time;
            obj.log.i = obj.log.i + 1;
        end
        
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
        
        function setInlet(obj, diameter, angle, Pr, Mr)
            obj.inletDiameter = diameter;
            obj.inletConeAngle = angle;
            obj.Pr_inlet = Pr;
            obj.Mr_inlet = Mr;
        end
        
        function setRDE(obj, od, id)
            obj.detOuterDiameter = od;
            obj.detInnerDiameter = id;
        end
        
        function [netThrust, Isp, Pmin, coneLength, q, P2, mdot_air] = solveFlowpath(obj, dt)
            % Get operating conditions
            M0 = obj.getMach();
            % Switch when we reach target mach
            if (M0 > obj.targetMach) || obj.qLatch
                phi = obj.phi_cruise;
                obj.qLatch =1;
            else
                phi = obj.phi_accel;
            end
            
            % Solve the inlet
            [M2, P2, T2, mdot_air, coneLength, q, ramDrag] = inlet.solveInlet(...
                obj.inletDiameter, obj.inletConeAngle, M0, obj.Pr_inlet,...
                obj.altitude, obj.Mr_inlet);
            
            [P0, ~] = obj.getLocalFlowParams();
            
            if (obj.fuelMass > 0)
                % Only burn if there is still fuel
                % Iterate until minimum chamber pressure is balanced
                Pmin_guess = obj.lastPmin;
                c_err = inf;
                c_maxErr = 100; % [Pa]
                c_step = obj.lastDP;
                c_lastDir = 1;
                tsteps = 1000;
                numDets = 1;
                while abs(c_err) > c_maxErr
                    % Get CEA detonation parameters
                    params = obj.cea.run('problem', 'det', 'p,atm', Pmin_guess / 101325, 't,k', T2, ...
                        'phi', phi, 'output', 'trans', 'reac', 'fuel' ,'C2H4', 'wt%', 100, 'oxid', ...
                        'Air', 'wt%', 100, 'end');
                    R = 8314 / params.output.burned.mw;
                    Pr = params.output.p_ratio;
                    Tr = params.output.t_ratio;
                    v_cj = params.output.det_vel;
                    gamma_det = params.output.burned.gamma;
                    Tmax = T2 * Tr; % [K] Temperature of burned gas
                    [Isp, Thrust, Pmin] = combustor.solveRDE(Pr, Pmin_guess, ...
                        Tmax, v_cj, R, obj.detOuterDiameter, obj.detInnerDiameter, ...
                        mdot_air, phi, gamma_det, tsteps, P0, numDets, M2);
                    
                    c_err = Pmin - Pmin_guess; % Pressure guess error
                    
                    if c_lastDir ~= sign(c_err)
                        c_step = c_step / 2;
                    end
                    
                    c_lastDir = sign(c_err);
                    
                    Pmin_guess = Pmin_guess + sign(c_err) * c_step;
                end
                obj.lastDP = max(abs(Pmin - obj.lastPmin), 10);
                obj.lastPmin = Pmin;
            else
                Thrust = 0;
            end
            
            % Update physics
            netThrust = (Thrust - ramDrag);
            % Try to increase Mach and q such that both reach their goal at
            % the same time
            if ~obj.qLatch
                currM = obj.getMach();
                currQ = obj.getDynamicPressure();
                machLeft = obj.targetMach - currM;
                qLeft = obj.targetQ - currQ;
                mqErr = inf;
                maxMqErr = 0.001; % Within 10% is close enough
                alpha = 0;
                step = 0.1; % Degrees
                lastDir =1;
                while abs(mqErr) > maxMqErr
                    acceleration = (netThrust - obj.getDrag - ...
                        (obj.getTotalMass * 9.81 * sind(alpha))) / ...
                        obj.getTotalMass;
                    nextV = obj.velocity + acceleration * dt;
                    nextAlt = obj.altitude + nextV * sind(alpha) * dt;
                    [nextQ, nextM] = obj.getTraj(nextV, nextAlt);
                    dQ = (nextQ - currQ) / qLeft;
                    dMach = (nextM - currM) / machLeft;
                    mqErr = (dMach - dQ);
                    if abs(mqErr) > maxMqErr
                        if sign(mqErr) ~= lastDir
                            step = step / 2;
                        end
                        lastDir = sign(mqErr);
                        
                        % mqErr > 1 => Increasing speed too much, should increase alpha
                        alpha = alpha - sign(mqErr) * step;
                    end
                end
            else
                alpha = 0;
            end
            acceleration = (netThrust - obj.getDrag - ...
                (obj.getTotalMass * 9.81 * sind(alpha))) / ...
                obj.getTotalMass;
            obj.velocity = obj.velocity + acceleration * dt;
            obj.altitude = obj.altitude + obj.velocity * sind(alpha) * dt;
            % Keep track of fuel
            obj.drainFuel(mdot_air, phi, dt);
        end
        
        function drainFuel(obj, mdot_air, phi, dt)
            FAR_stoich = 1/14.78733504955836;
            FAR = FAR_stoich * phi;
            spentFuel = mdot_air * FAR * dt;
            obj.fuelMass = obj.fuelMass - spentFuel;
            if (obj.fuelMass < 0)
                obj.fuelMass = 0;
            end
        end
        
        function [q, M] = getTraj(obj, vel, alt)
            % Get the desired q and M for a velocity at an altitude
            P0 = interp1(obj.alth, obj.altPa, alt, 'linear'); % Freestream static pressure
            T0 = interp1(obj.alth, obj.altT, alt, 'linear'); % Freestream static temperature
            a = sqrt(obj.gamma_air * obj.R_air * T0);
            M = vel / a;
            q = 0.5 * obj.gamma_air * P0 * M^2;
        end
        
        function lift = getLift(obj)
            lift = obj.Cl * obj.getDynamicPressure() * obj.liftingArea;
        end
        
        function drag = getDrag(obj)
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
            M = obj.velocity / obj.getSonicVelocity();
        end
        
        function a = getSonicVelocity(obj)
            [~, T0] = obj.getLocalFlowParams();
            a = sqrt(obj.gamma_air * obj.R_air * T0);
        end
        
        function [P0, T0] = getLocalFlowParams(obj)
            P0 = interp1(obj.alth, obj.altPa, obj.altitude, 'linear'); % Freestream static pressure
            T0 = interp1(obj.alth, obj.altT, obj.altitude, 'linear'); % Freestream static temperature
        end
    end
    
end

