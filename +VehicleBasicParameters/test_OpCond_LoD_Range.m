clc
clear

q = 1500; %[psf]
M = 6;
b = 24; %[in]
c_r = 90; %[in] Root chord (or length of combustion)
c_t = 30; %[in] tip chord
weight_gross = 400; %[lbm]
weight_zerofuel = 325; %[lbm]
I_sp = 2000; %[s]

[Mach,altitude,P0,T0,rho0,v0] = getOperatingConditions(q,M)

[LoD, C_L, C_D, Aero_Ref_Area]=getLoD(q,M,b,c_r,c_t,weight_gross)

[Range]=getRange(v0, C_L, C_D, I_sp, weight_gross, weight_zerofuel,altitude)