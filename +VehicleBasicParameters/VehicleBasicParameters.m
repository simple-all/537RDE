%% Vehicle Basic Parameters
% Andrew McClaskey

clc
clear
close all
format compact

%test

%Constants
b_max = 24; %[in]
s = b_max/2; %[in] half span
c_r = 90; %[in] Root chord (or length of combustion)
c_t = 20; %[in] tip chord
Weight_gross = 600; %[lbm]
%Inlet_Diameter = 11; %[inches]
TSFC = 0.6; %Need a real value; used in Range equation
rho_liqHyd = 4.423; %[lb/cu.ft]

%Flight Conditions
M = 6;
q = 1500; %[psf]
gamma = 1.4; %Specific Heat of Air
R = 1716; %Gas Constant [ft.lbf/R.slug]
h_max = 2400; %[100*ft]
for i = 1:length(M) %for each M vary altitude to get rho_error less than error
    for j = 1:h_max %Altitude iteration
        error = 0.001;
        h(j) = 100*j; %increase altitude by 100 feet for each iteration
        [t(j),p_psf(j),r(j),height(j)]=atmosphere(h(j),0);
        a(j) = sqrt(gamma*R*t(j)); %[ft/s]
        v(j) = M(i) * a(j); %[ft/s]
        rho(j) = (2*q) /(v(j)^2);
        rho_error(j) = (r(j)-rho(j))/r(j);
            if rho_error(j) < error
               Mach(i) = M(i);
               altitude(i) = h(j); %[ft]
               press_psf(i) = p_psf(j); %{psf]
               temp(i) = t(j); %[R]
               v0(i) = v(j); %[ft/s]
               rho0(i) = rho(j); %[slug/ft^3]
               break
            end
    end
end

%% Planform
Aero_Ref_Area = (c_t + c_r)*s; %[sq. inches]
sweep = atand((c_r - c_t)/s); %[deg] wing sweep angle
wing_x = [0 0 s s];
wing_y = [0 c_r c_t 0];
figure(1)
plot(wing_x, wing_y)
axis([0 c_r/3 0 c_r])

%% Lift, Drag, and L over D

% Drag Calculation
% Current Analysis assumes alpha = 0
% The drag coefficients were pulled from GHV and curve fit to allow for
% calculating the drag at different Mach numbers. I can expand this to use
% multiple alphas if we think that is neccesary
C_D = 0.00023*M^2 - 0.00282*M + 0.01869;
D = (q/144)*Aero_Ref_Area*C_D

% Lift Calculation
% L = W; determine the lift coefficient required to support the weight at
% the specified dynamic pressure
C_L = Weight_gross/((q/144)*Aero_Ref_Area);

%Calculate L/D using C_L and C_D
L_D = C_L/C_D

%% Tank Geometry
% Split sphere and cylinder "pill capsle" style tank. d is diameter, r is
% radius, l is total length of tank, a is length of center cylinder
% section.
% Currently sized for three tanks, one long center, two short (each side of
% long tank). Tanks are parrallell to longitudinal axis.
DesiredTankPressure = 1.3; %[MPa] %Target Hydrogen Fuel Tank Pressure
P_tank = DesiredTankPressure*(145.038) %[psi]
d = [7 5 5]; %[in] Diameter of tank
r = d./2; %[in]
l = [60 30 30]; %[in] 3 tanks
a = l - (2.*d) %[in]
V_tank = pi.*r.^2.*((4/3).*r+a) %[cu.in]
W_tank = (3*P_tank.*V_tank)./(1*10^6) %[lbm] %Burst pressure is 3 times max operating
W_liqHyd = (rho_liqHyd/1728).*sum(V_tank)
Weight_zerofuel = Weight_gross - W_liqHyd; %[lbm] 

%% Range
% Assumed the initial weight and zero fuel weight as specified in the
% "constrants" section. Need to discuss the zero fuel weight.
R_ft = 2*sqrt(2/(rho0*Aero_Ref_Area))*(1/(TSFC/3600))*(sqrt(C_L)/C_D)*(sqrt(Weight_gross)-sqrt(Weight_zerofuel));
R_nm = R_ft/6076
%% Maybe Junk
%Read in Excel file of Coefficient of Drag data
%Currently uses phi = 0.5 data. This can be updated once we know the real
%equivalence ratio
% xlsread('DragCoefficient.xlsx');
% AoA = xlsread('DragCoefficient.xlsx', 'B4:B9');
% C_D_M4 = xlsread('DragCoefficient.xlsx', 'C4:C9');
% C_D_M45 = xlsread('DragCoefficient.xlsx', 'D4:D9');
% C_D_M5 = xlsread('DragCoefficient.xlsx', 'E4:E9');
% C_D_M55 = xlsread('DragCoefficient.xlsx', 'F4:F9');
% C_D_M6 = xlsread('DragCoefficient.xlsx', 'G4:G9');
% C_D_M65 = xlsread('DragCoefficient.xlsx', 'H4:H9');
% C_D_M7 = xlsread('DragCoefficient.xlsx', 'I4:I9');
% % Calculate Drag for range of alpha for each mach number 
% % D = q*S*C_D
% D4 = (q/144)*Aero_Ref_Area.*C_D_M4; %[lbf}
% D45 = (q/144)*Aero_Ref_Area.*C_D_M45; %[lbf}
% D5 = (q/144)*Aero_Ref_Area.*C_D_M5; %[lbf}
% D55 = (q/144)*Aero_Ref_Area.*C_D_M55; %[lbf}
% D6 = (q/144)*Aero_Ref_Area.*C_D_M6 %[lbf}
% D65 = (q/144)*Aero_Ref_Area.*C_D_M65; %[lbf}
% D7 = (q/144)*Aero_Ref_Area.*C_D_M7; %[lbf}
% 
% output_drag_table = [D4 D45 D5 D55 D6 D65 D7];
% xlswrite('drag.xlsx',output_drag_table)
