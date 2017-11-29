%% getRange
% Calculates range using the Breguet equation
% Inputs
% 1. Free Stream Velocity (v0) [ft/s]
% 2. Lift Coefficient (C_L) [Nondim]
% 3. Drag Coefficient (C_D) [Nondim]
% 4. Specific Impulse (I_sp) [seconds]
% 5. Gross Weight (weight_gross) [pounds]
% 6. Zero Fuel Weight (weight_zerofuel) [pounds] 

function [Range]=getRange(v0, C_L, C_D, I_sp, weight_gross, weight_zerofuel)
Range_ft = v0*(C_L/C_D)*I_sp*log(weight_gross/weight_zerofuel); %[ft]
Range = Range_ft/6076; %[Nautical Miles]
end