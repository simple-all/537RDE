clc
clear

q = 1500; %[psf]
M = 6;
b = 24; %[in]
c_r = 90; %[in] Root chord (or length of combustion)
c_t = 30; %[in] tip chord
I_sp = 1000; %[s]
wing_loading = 30; %[psf]


[ScalingFactor, Aero_Ref_Area, weight_gross, weight_structures, weight_misc, weight_fluid, weight_zerofuel]=getweights(b, c_r, c_t, wing_loading)

[Mach,altitude,P0,T0,rho0,v0] = getOperatingConditions(q,M)

[LoD, C_L, C_D]=getLoD(q,M,Aero_Ref_Area, weight_gross)

[Range_cruise, Range_glide, Range]=getRange(v0, C_L, C_D, I_sp, weight_gross, weight_zerofuel,altitude)