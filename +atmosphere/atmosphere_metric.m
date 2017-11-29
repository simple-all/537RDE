function [T, P, rho, HGeo] = atmosphere_metric(Hvector, isGeo)
%ATMOSPHERE_METRIC Easy conversion of the atmosphere from SAE to SI

% Inputs:
% Hvector: [m] Vector of altitudes
% isGeo: Flag, 1 = geometric altitude, 0 = geopotential altitude

% Outputs:
% T: [k] Temperature
% P: [Pa] Pressure
% rho: [kg/m^3] Density
% HGeo: [m] Geopotential altitude vector

% Convert inputs to SAE
Hvect_SAE = Hvector / (12 * 0.0254);

% Get atmospheric data
[T, P, rho, HGeo] = atmosphere.atmosphere(Hvect_SAE, isGeo);

% Convert back to SI
T = T * (5/9);
P = P * 47.8802588888888888; 
rho = rho * 515.379;
HGeo = HGeo * (12 * 0.0254);

end

