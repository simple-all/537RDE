function [Isp_RDE, Thrust, Pmin, Pavg, Tavg, Pmax] = solveRDE(Pr, Pmin, Tmax, v_cj, R, D_outer, D_inner, mdot_air, phi, gamma, tsteps, P0, numDets, M)
%solveRDE Time-stepped integration solver for an RDE
% Uses the Stechmann model

FAR_stoich = 1/14.78733504955836;
FAR = FAR_stoich * phi;

A_t = pi *((D_outer / 2)^2 - (D_inner / 2)^2); % [m^2] Throat area

% Exponential decay factor
tau_c = (pi * D_outer) / (v_cj * numDets); % Time between detonation fronts
lambda = log(Pr)/tau_c;


F_gamma = ((2 * gamma^2) / (gamma - 1)) * ((2 / (gamma + 1))^((gamma + 1) / (gamma - 1)));

eta_max = (pi * (D_outer / 2)^2) / A_t; % Maximum expansion ratio

a_star_c = A_t / aeroBox.isoBox.calcARatio(M, gamma);
M_e = aeroBox.isoBox.machFromAreaRatio((pi * (D_outer / 2)^2) / a_star_c, gamma, 1);

dt = tau_c / tsteps;
err = inf;
step = 1e4;
dir = 1;
eta_nozzle = 0.93;
while abs(err) > 1e-4;
    time = 0;
    i = 1;
    while time <= tau_c
        tt(i) = time;
        P(i) = Pc(time);
        T(i) = Tc(time);
        
        % Calculate exit pressure out of the nozzle
        Pt = P(i) * ((1 + (((gamma - 1) / 2) * M^2))^(gamma / (gamma - 1)));
        Tt = T(i) * ((1 + (((gamma - 1) / 2) * M^2)));
        P_e = Pt / ((1 + (((gamma - 1) / 2) * M_e^2))^(gamma / (gamma - 1)));
        T_e = Tt / ((1 + (((gamma - 1) / 2) * M_e^2)));
        
        P_ei(i) = P_e;
        T_ei(i) = T_e;
        
        c_star(i) = sqrt(gamma * R * T(i)) / (gamma * sqrt((2 / (gamma + 1))^((gamma + 1) / (gamma - 1))));
        
        if (P_e < P0)
            cf(i) = sqrt(F_gamma * (1 - ((P0 / P(i))^((gamma - 1) / gamma))));
        else
            cf(i) = sqrt(F_gamma * (1 - ((P_e / P(i))^((gamma - 1) / gamma)))) + eta_max * ((P_e / P(i)) - (P0 / P(i)));
        end
        mdot(i) = (P(i) * A_t) / c_star(i);
        Isp(i) = (cf(i) * c_star(i) * sqrt(eta_nozzle)) / 9.81;
        F(i) = mdot(i) * Isp(i) * 9.81;
        
        i = i + 1;
        time = time + dt;
    end
    
    massFlow = mean(mdot);
    err = massFlow - (mdot_air * (FAR + 1));
    
    if sign(err) ~= dir
        step = step / 2;
    end
    dir = sign(err);
    % err > 0, Minimum pressure too high, decrease
    % err < 0, Minimum pressure too low, increase
    
    Pmin = Pmin + -sign(err) * step;
end
m_cyc = mean(mdot) * (1 + FAR) * tau_c; % [kg] Mass per cycle

vx=mean(mdot)/(0.5*((((D_outer / 2))^2)-(((D_inner / 2))^2))*combustor.K(gamma,R,Pc(0),Pr,Tmax));
if vx < 0
    warning('Axial velocity less than 0!');
end
Isp_RDE = 0;
for i = 1:numel(F)
    Isp_RDE = Isp_RDE + F(i) * dt;
end
Isp_RDE = Isp_RDE / m_cyc;

Thrust = mean(F);
Pavg = mean(P);
Tavg = mean(T);
Pmax = max(P);

    function p = Pc(time)
        p = Pr * Pmin * exp(-lambda * time);
    end

    function t = Tc(time)
        t = Tmax * (Pc(time) / Pc(0))^(1 - (1 / gamma));
    end

end

