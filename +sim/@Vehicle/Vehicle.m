classdef Vehicle < handle
    %@VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        wetMass % [kg] Initial wet mass of vehicle
        fuelMassFraction;
        phi_accel;
        phi_cruise;
        targetMach;
        targetQ;
        
        % Sim parameters
        time = 0; % [s] Current time in simulation/flight
        lastRunTime = 0;
        
        % Movement parameters
        velocity; % [m/s] Flight velocity
        altitude; % [m] Flight altitude
        
        % Inlet parameters
        inletDiameter; % [m] Outer diameter of inlet
        inletGap;
        isolatorArea;
        inletConeAngle; % [deg] Half-angle of inlet cone
        Pr_inlet; % Pressure recovery factor of the inlet
        Mr_inlet = 0.33; % Mach ratio of the inlet-isolator system
        
        % Engine parameters
        detOuterDiameter; % [m] Outer diameter of engine
        detInnerDiameter; % [m] Inner diameter of engine
        detArea;
        
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
        lastPmin = 213750; 
        lastDP = 1000;
        qLatch = 0;
        distTraveled = 0;
    end
    
    properties (Constant)
        % Atmospheric parameters
        gamma_air = 1.4; % Ratio of specific heats
        R_air = 287.058; % [J/kg*K] Air gas constant
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
            
            [netThrust, Isp, Pmin, coneLength, q, P3, mdot_air, alpha, Pavg, Tavg, Ptc] = obj.solveFlowpath(dt);
            
            i = obj.log.i;
            obj.log.netThrust(i) = netThrust;
            obj.log.Isp(i) = Isp;
            obj.log.Pmin(i) = Pmin;
            obj.log.coneLength(i) = coneLength;
            obj.log.q(i) = q;
            obj.log.P3(i) = P3;
            obj.log.mdot_air(i) = mdot_air;
            obj.log.M(i) = obj.getMach();
            obj.log.altitude(i) = obj.altitude;
            obj.log.alpha(i) = alpha;
            obj.log.distance(i) = obj.distTraveled;
            obj.log.fuelMass(i) = obj.fuelMass;
            obj.log.Tavg(i) = Tavg;
            obj.log.Pavg(i) = Pavg;
            obj.log.Ptc(i) = Ptc;
            
            obj.lastRunTime = obj.time;
            obj.log.i = obj.log.i + 1;
        end
        
        function setPhis(obj, phi_accel, phi_cruise)
            obj.phi_accel = phi_accel;
            obj.phi_cruise = phi_cruise;
        end
        
        function setTarget(obj, targetM, targetQ)
            obj.targetMach = targetM;
            obj.targetQ = targetQ;
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
        
        function setInlet(obj, diameter, inletGap, angle, Pr, Mr)
            obj.inletDiameter = diameter;
            obj.inletGap = inletGap;
            obj.isolatorArea = pi * ((obj.inletDiameter / 2)^2 - ...
                ((obj.inletDiameter - (obj.inletGap))/2)^2);
            obj.inletConeAngle = angle;
            obj.Pr_inlet = Pr;
            obj.Mr_inlet = Mr;
        end
        
        function setRDE(obj, od, id)
            obj.detOuterDiameter = od;
            obj.detInnerDiameter = id;
            obj.detArea = pi * ((od/2)^2 - (id / 2)^2);
        end
        
        function [netThrust, Isp, Pmin, coneLength, q, P3, mdot_air, alpha, Pavg, Tavg, Ptc] = solveFlowpath(obj, dt)
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
            [M2, P2, T2, mdot_air, coneLength, q, ramDrag, Pt2, Tt2] = inlet.solveInlet(...
                obj.inletDiameter, obj.inletConeAngle, M0, obj.Pr_inlet,...
                obj.altitude, obj.Mr_inlet);
            
            [P0, ~] = obj.getLocalFlowParams();
            
            if (obj.fuelMass > 0)
                % Only burn if there is still fuel
                % Iterate until minimum chamber pressure is balanced
                
                A_star_isolator = obj.isolatorArea / ...
                    aeroBox.isoBox.calcARatio(M2, obj.gamma_air);
                
                M3 = aeroBox.isoBox.machFromAreaRatio(...
                    obj.detArea / A_star_isolator, obj.gamma_air, 1);
                P3 = aeroBox.isoBox.calcStaticPressure('mach', M3, ...
                    'gamma', obj.gamma_air, 'Pt', Pt2);
                T3 = aeroBox.isoBox.calcStaticTemp('mach', M3, ...
                    'gamma', obj.gamma_air, 'Tt', Tt2);
                
                Pmin_guess = obj.lastPmin;
                c_err = inf;
                c_maxErr = 100; % [Pa]
                c_step = obj.lastDP;
                c_lastDir = 1;
                tsteps = 1000;
                numDets = 1;
                while abs(c_err) > c_maxErr && c_step > 1
                    % Get CEA detonation parameters
                    params = obj.cea.run('problem', 'det', 'p,atm', Pmin_guess / 101325, 't,k', T3, ...
                        'phi', phi, 'output', 'trans', 'reac', 'fuel' ,'C2H4', 'wt%', 100, 'oxid', ...
                        'Air', 'wt%', 100, 'end');
                    R = 8314 / params.output.burned.mw;
                    Pr = params.output.p_ratio;
                    Tr = params.output.t_ratio;
                    v_cj = params.output.det_vel;
                    gamma_det = params.output.burned.gamma;
                    Tmax = T3 * Tr; % [K] Temperature of burned gas
                    [Isp, Thrust, Pmin, Pavg, Tavg, Pmax] = combustor.solveRDE(Pr, Pmin_guess, ...
                        Tmax, v_cj, R, obj.detOuterDiameter, obj.detInnerDiameter, ...
                        mdot_air, phi, gamma_det, tsteps, P0, numDets, M3);

                    Ptc = aeroBox.isoBox.calcStagPressure('mach', M3, 'Ps', Pmax, 'gamma', gamma_det);
                    c_err = Pmin - Pmin_guess; % Pressure guess error
                    
                    if c_lastDir ~= sign(c_err)
                        c_step = c_step / 2;
                    end
                    
                    c_lastDir = sign(c_err);
                    
                    Pmin_guess = Pmin_guess + sign(c_err) * c_step;
                end
                % Speed boost on future solution generation
                obj.lastDP = max(abs(Pmin - obj.lastPmin), c_maxErr / 2);
                obj.lastPmin = Pmin;
            else
                Thrust = 0;
            end
            
            % Update physics
            netThrust = (Thrust - ramDrag);
            % Try to increase Mach and q such that both reach their goal at
            % the same time
            if ~obj.qLatch
                % Goal is to reach target mach and dynamic pressure at the
                % same time
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
                        obj.getLift*sind(alpha) - ...
                        (obj.getTotalMass * 9.81 * sind(alpha))) / ...
                        obj.getTotalMass;
                    nextV = obj.velocity + acceleration * dt;
                    nextAlt = obj.altitude + nextV * sind(alpha) * dt + ...
                        ((obj.getLift * cosd(alpha) - obj.getTotalMass * 9.81) / obj.getTotalMass) * dt;
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
                % Goal is to maintain dynamic pressure during cruise
                qErr = inf;
                maxMqErr = 1; % Within 10 Pa
                alpha = 0;
                step = 0.1; % Degrees
                lastDir =1;
                while abs(qErr) > maxMqErr
                    acceleration = (netThrust - obj.getDrag - ...
                        obj.getLift*sind(alpha) - ...
                        (obj.getTotalMass * 9.81 * sind(alpha))) / ...
                        obj.getTotalMass;
                    nextV = obj.velocity + acceleration * dt;
                    nextAlt = obj.altitude + nextV * sind(alpha) * dt + ...
                        ((obj.getLift * cosd(alpha) - obj.getTotalMass * 9.81) / obj.getTotalMass) * dt;
                    nextQ = obj.getTraj(nextV, nextAlt);
                    qErr = nextQ - obj.targetQ;
                    if abs(qErr) > maxMqErr
                        if sign(qErr) ~= lastDir
                            step = step / 2;
                        end
                        lastDir = sign(qErr);
                        
                        % qErr > 1 => Increasing speed too much, should increase alpha
                        alpha = alpha + sign(qErr) * step;
                    end
                end
            end
            acceleration = (netThrust - obj.getDrag - ...
                        obj.getLift*sind(alpha) - ...
                        (obj.getTotalMass * 9.81 * sind(alpha))) / ...
                        obj.getTotalMass;
            obj.velocity = obj.velocity + acceleration * dt;
            obj.altitude = obj.altitude + obj.velocity * sind(alpha) * dt + ...
                ((obj.getLift * cosd(alpha) - obj.getTotalMass * 9.81) / obj.getTotalMass) * dt;
            obj.distTraveled = obj.distTraveled + obj.velocity * cosd(alpha) * dt;
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
            [T0, P0] = atmosphere.atmosphere_metric(alt, 1);
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
            [T0, P0] = atmosphere.atmosphere_metric(obj.altitude, 1);
        end
    end
    
end

