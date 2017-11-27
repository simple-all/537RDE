%% getRange
% Calculates range
% Inputs
% 1. Free Stream Density (rho0) [slugs/ft3]
% 2. Wing Reference Area (Aero_Reference_Area) [Square inches]
% 3. Lift Coefficient (C_L) [Nondim]
% 4. Drag Coefficient (C_D) [pounds]
% 5. Gross Weight (weight_gross) [pounds]
% 6. Zero Fuel Weight (weight_zerofuel) [pounds] 


function [Range]=getRange(rho0,Aero_Ref_Area,TSFC, C_L, C_D, weight_gross, weight_zerofuel)
R_ft = 2*sqrt(2./(rho0*Aero_Ref_Area))*(1/(TSFC/3600)).*(sqrt(C_L)./C_D)*(sqrt(weight_gross)-sqrt(weight_zerofuel)); %[feet]
Range = R_ft./6076; %[Nautical Miles]
end