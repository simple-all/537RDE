%% AAE 537 Project
%  Regenerative Cooling of Combustor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear all

% Combustor Properties
Tc = 5400; %[R]
Te = 180; %[R] Temperature of ethylene
%boiling boint of ethylene at 1atm is 305 R

Pc = 188.549; %[psi]
Min = 2; %Mach number entering combustor
k = 1.1797; % gamma
T0 = Tc*(1+(k-1)/2*Min^2); %stagnation temperature of combustor
Dt = 6.83; %[in] combustor diameter
L = 10/12; %[ft] length of combustor
At_A = 1; %Area Ratio
Cp = 0.549; %[Btu/lbm.F]
Pr = 0.4882; %Prandtl Number in combustor
Rt = 1; %radius of throat curvature
r2 = 2.55; %[in]
th = 0.05; %[ft] wall thickness between coolant channels and combustor
g = 32.2*12; %[in/s^2]
K_gas = 0.00104; %[Btu/s.ft.R]
n = 50; %number of coolant channels
cstar = 4238.845*12; %[ft/s]

%figuring out channel geometry
r_in = r2 - th; % [in] [radius of bottom of channel, using 0.05 wall thickness]
circ = 2*pi*r_in; %[in] circumference;

%use a channel width of 0.1 inches
height = 0.125; %[in]
width = 0.125; %[in]
width_tot = width*n; %using 20 channels to cool
left = circ - width_tot; %[in] lefotver space to be spread between channels
space = left/n; %[in] space between each channel
A_ch = height*width; %[in^2] channel area
%A_ch = A_ch/144;

%Thermal Conductivity [Btu/ft.s.F]
%Inconel
%K = 0.00104; %[Btu/ft.s.F]
%Hastelloy
K = 0.00451; %[Btu/ft.s.F]
%Copper 101
%K = 0.0628;

%Recovery Temperature
r = Pr^(1/3); %correction factor
Tr = Tc*(1+(k-1)/2*r*Min^2);

%Adiabatic Wall Temperature
Taw = Tc*((1+r*(k-1)/2*Min^2)/(1+(k-1)/2*Min^2));

T_eth(1) = Te;

%% Loop for finding proper Twg and heat flux
for z = 1:100
% Guess Twg
Twg = 100; %[R]

% set arbitrary error value
err = 1000;

%Begin while loop
while err > 0.5

%find qdot with guess of Twg
mu = 4.7867e-6; %[lb/in.s]
B = 1+0.5*(k-1)*Min^2;
sig = 1/(0.5+0.5*(Twg/Tc)*B^0.68);
%calculate hg [Btu/in^2.s.F]
hg = (0.026*mu^0.2 * Cp*(Pc*g)^0.8)/...
    (Dt^0.2 * Pr^0.6 * cstar^0.8) * At_A^0.9 * sig;


qd = hg*(Tr - Twg);
qd = qd*144; %[Btu/s.ft^2]

%Calculate temerature of wall on coolant side
Twc = Twg - qd*(th/12)/K; %[R]


%Ethylene properties
mdot_eth = 0.154*0.43; %[lb/s] mass flowrate of ethylene
rho_eth = 35.437; %[lb/ft^3] liquid density of ethylene
Pr_eth = 2.25; %Prandtl Number
Cp_eth = 1.236; %[Btu/(lbF)]
v_eth = (mdot_eth/n)/(rho_eth*(A_ch/144)); %[ft/s] velocity through a single channel
mu_eth = 3.06e-5; %[lb/ft.s] viscosity of ethylene
K_eth = ((mu_eth*Cp_eth)/Pr_eth); %[Btu/in.s.F]
Re_eth = rho_eth*v_eth*L/mu_eth; %Reynold's Number
Dh = 4*A_ch/(2*height+2*width);

%Ethylene Nussselt Number
Nu_eth = 0.0214*Re_eth^0.8 * Pr^0.4;

%hl and heat flux for coolant side
hl = Nu_eth*K_eth/Dh; %[Btu/s.in^2.R]
qd_cool = hl*(Twc - Te); %[Btu/s.in^2]
qd_cool = qd_cool * 144; %[Btu/s.ft^2]

%find percent error between two qd values
err = abs(qd_cool - qd)/qd * 100;

Twg = Twg + 1;
end

%Calculate new temperature of the liquid
dx = 0.1;
dT = qd*((dx/12)* (width/12))/(mdot_eth*Cp_eth);
Te = Te + dT;
%ending temperature of ethylene in -177 degrees fahrenheit
%boiling point of liquid ethylene is -154 degrees fahrenheit

T_eth(z+1) = Te;
Twall(z) = Twg;


end