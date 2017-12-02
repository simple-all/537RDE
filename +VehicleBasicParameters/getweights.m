%function that determines weight of vehicle based on process Prashanth laid
%out. Inputs are...
% 1. Wing span (b) [inches]
% 2. Root Chord (c_r) [inches]
% 3. Tip Chord (c_t) [inches]
% 4. Wing Loading (wing_loading) [psf]

function [ScalingFactor, Aero_Ref_Area, weight_gross, weight_structures, weight_misc, weight_fluid, weight_zerofuel]=getweights(b, c_r, c_t, wing_loading)
%% Aero Plan form
s = b/2; %[in] half span
Aero_Ref_Area = (c_t + c_r)*s; %[sq. in]
sweep = atand((c_r - c_t)/s); %[deg] wing sweep angle

%Scaling Factor for all weight determiination
GHV_Aero_Ref_Area = 6026; %[sq.in]
ScalingFactor = Aero_Ref_Area/GHV_Aero_Ref_Area; %[Nondim]

%Determine gross weight from wing loading and aero reference area
weight_gross = (wing_loading/144)*Aero_Ref_Area; %[lbm]

structural_mass_fraction = 0.496; %[Nondim] from Prashanth
misc_mass_fraction = round(-0.0001*weight_gross + 0.3226,3); %[Nondim] Curve fit of data from Prashanth
fluid_mass_fraction = 1 - structural_mass_fraction - misc_mass_fraction; %[Nondim] Determine fluids mass fraction once structural and misc weights are known
ethylene_mass_fraction = fluid_mass_fraction*.904; %[Nondim] %Determine propellant mass fraction
other_fluids_mass_fraction = fluid_mass_fraction*.0167; %[Nondim] Determine other fluids mass fraction

weight_structures = weight_gross*structural_mass_fraction; %[lbm]
weight_misc = weight_gross*misc_mass_fraction; %[lbm]
weight_fluid = weight_gross*fluid_mass_fraction; %[lbm]
weight_propellant = round(weight_gross*ethylene_mass_fraction,1); %[lbm]
weight_zerofuel = weight_gross - weight_propellant; %[lbm]

end