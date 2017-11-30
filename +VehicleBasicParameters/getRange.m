%% getRange
% Calculates range using the Breguet equationand glide range equation
% Inputs
% 1. Free Stream Velocity (v0) [ft/s]
% 2. Lift Coefficient (C_L) [Nondim]
% 3. Drag Coefficient (C_D) [Nondim]
% 4. Specific Impulse (I_sp) [seconds]
% 5. Gross Weight (weight_gross) [pounds]
% 6. Zero Fuel Weight (weight_zerofuel) [pounds] 
% 7. Altitude (altitude) [ft]

function [Range_cruise_ideal,Range_cruise_expected]=getRange(v0, LoD_ideal, LoD_expected, I_sp, weight_gross, weight_zerofuel)
% Powered Flight
Range_cruise_ft_ideal = v0*(LoD_ideal)*I_sp*log(weight_gross/weight_zerofuel); %[ft]
Range_cruise_ideal = Range_cruise_ft_ideal/6076; %[Nautical Miles]
Range_cruise_ft_expected = v0*(LoD_expected)*I_sp*log(weight_gross/weight_zerofuel); %[ft]
Range_cruise_expected = Range_cruise_ft_expected/6076; %[Nautical Miles]
% Gliding Flight
% Range_glide_ft = (C_L/C_D)*altitude; %[ft]
% Range_glide = Range_glide_ft/6076; %[Nautical Miles]
% Total Range
% Range = Range_cruise + Range_glide;
end