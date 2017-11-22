%% Vehicle Basic Parameters
% Andrew McClaskey

clc
clear
close all
format compact


%Constants
b_max = 24; %[in]
s = b_max/2; %[in] half span
c_r = 90; %[in] Root chord (or length of combustion)
c_t = 30; %[in] tip chord
Weight_gross = 600; %[lbm]
%Inlet_Diameter = 11; %[inches]
TSFC = 0.6; %Need a real value; used in Range equation
rho_liqHyd = 4.423; %[lb/cu.ft]
gamma = 1.4; %Specific Heat of Air
R = 1716; %Gas Constant [ft.lbf/R.slug]


%% Tank Geometry
% Split sphere and cylinder "pill capsule" style tank. d is diameter, r is
% radius, l is total length of tank, a is length of center cylinder
% section.
% Currently sized for three tanks, one long center, two short (each side of
% long tank). Tanks are parrallell to longitudinal axis.
DesiredTankPressure = 0.17; %[MPa] %Target Hydrogen Fuel Tank Pressure
P_tank = DesiredTankPressure*(145.038) %[psi]
d = [7 5 5]; %[in] Diameter of tank
r = d./2; %[in]
l = [55 30 30]; %[in] 3 tanks
a = l - (2.*d); %[in]
V_tank = pi.*r.^2.*((4/3).*r+a) %[cu.in]
W_tank = (3*P_tank.*V_tank)./(1*10^6) %[lbm] %Burst pressure is 3 times max normal operating
W_liqHyd = (rho_liqHyd/1728).*sum(V_tank)
Weight_zerofuel = Weight_gross - W_liqHyd; %[lbm]
tank1_x = [0  0   d(1) d(1)]; %x-coordinate for tank1
tank1_y = [0 l(1) l(1)  0]; %y-coordinate for tank1
tank2_x = [0  0   -d(2) -d(2)]; %x-coordinate for tank2
tank2_y = [0 l(2) l(2)  0]; %y-coordinate for tank2

%% Planform
Aero_Ref_Area = (c_t + c_r)*s; %[sq. inches]
sweep = atand((c_r - c_t)/s); %[deg] wing sweep angle
wing_x = [0 0 s s];
wing_y = [0 c_r c_t 0];
figure(1)
plot(wing_x, wing_y,'k', -wing_x, wing_y,'k', tank1_x, tank1_y, tank2_x, tank2_y)
axis([-c_r/3 c_r/3 0 c_r])
title('Planform Area')
xlabel('x-axis')
ylabel('y-axis')


%% Flight Conditions
q = [1000 1500 2000]; %[psf]
for k = 1:length(q)
M = [4:0.1:10];
h_max = 2400; %[100*ft]
for m = 1:length(M) %for each M vary altitude to get rho_error less than error
    for j = 1:h_max %Altitude iteration
        error = 0.001;
        h(j) = 100*j; %increase altitude by 100 feet for each iteration
        [t(j),p_psf(j),r(j),height(j)]=atmosphere.atmosphere(h(j),0);
        a(j) = sqrt(gamma*R*t(j)); %[ft/s]
        v(j) = M(m) * a(j); %[ft/s]
        rho(j) = (2*q(k)) /(v(j)^2);
        rho_error(j) = (r(j)-rho(j))/r(j);
            if rho_error(j) < error
               Mach(k,m) = M(m);
               altitude(k,m) = h(j); %[ft]
               P0_psf(k,m) = p_psf(j); %{psf]
               T0(k,m) = t(j); %[R]
               v0(k,m) = v(j); %[ft/s]
               rho0(k,m) = rho(j); %[slug/ft^3]
               break
            end
    end
end


%% Lift, Drag, and L over D
% Drag Calculation
% Current Analysis assumes alpha = 0
% The drag coefficients were pulled from GHV and curve fit to allow for
% calculating the drag at different Mach numbers. I can expand this to use
% multiple alphas if we think that is neccesary. Data that was curve fit
% was for M = 4 to 6.
C_D = 0.00023*M.^2 - 0.00282.*M + 0.01869;
D = (q(k)/144)*Aero_Ref_Area.*C_D;
figure(2)
plot(M,D)
hold on
xlabel('Mach No.')
ylabel('Drag')
title('Drag Polar')
legend(sprintf('q = %i',q(1)),sprintf('q = %i',q(2)),sprintf('q = %i',q(3)) )

% Lift Calculation
% L = W; determine the lift coefficient required to support the weight at
% the specified dynamic pressure
C_L(k) = Weight_gross/((q(k)/144)*Aero_Ref_Area);
% Calculate L/D using C_L and C_D
L_D = C_L(k)./C_D;
figure(3)
plot(M, L_D)
hold on
xlabel('Mach No.')
ylabel('L/D')
title('Lift to Drag Polar')
legend(sprintf('q = %i',q(1)),sprintf('q = %i',q(2)),sprintf('q = %i',q(3)) )


end
hold off
hold off

% %% Range
% % Assumed the initial weight and zero fuel weight as specified in the
% % "constrants" section. Need to discuss the zero fuel weight.
% R_ft = 2*sqrt(2./(rho0*Aero_Ref_Area))*(1/(TSFC/3600)).*(sqrt(C_L)./C_D)*(sqrt(Weight_gross)-sqrt(Weight_zerofuel));
% R_nm = R_ft./6076
% %% Maybe Junk
% %Read in Excel file of Coefficient of Drag data
% %Currently uses phi = 0.5 data. This can be updated once we know the real
% %equivalence ratio
% % xlsread('DragCoefficient.xlsx');
% % AoA = xlsread('DragCoefficient.xlsx', 'B4:B9');
% % C_D_M4 = xlsread('DragCoefficient.xlsx', 'C4:C9');
% % C_D_M45 = xlsread('DragCoefficient.xlsx', 'D4:D9');
% % C_D_M5 = xlsread('DragCoefficient.xlsx', 'E4:E9');
% % C_D_M55 = xlsread('DragCoefficient.xlsx', 'F4:F9');
% % C_D_M6 = xlsread('DragCoefficient.xlsx', 'G4:G9');
% % C_D_M65 = xlsread('DragCoefficient.xlsx', 'H4:H9');
% % C_D_M7 = xlsread('DragCoefficient.xlsx', 'I4:I9');
% % % Calculate Drag for range of alpha for each mach number 
% % % D = q*S*C_D
% % D4 = (q/144)*Aero_Ref_Area.*C_D_M4; %[lbf}
% % D45 = (q/144)*Aero_Ref_Area.*C_D_M45; %[lbf}
% % D5 = (q/144)*Aero_Ref_Area.*C_D_M5; %[lbf}
% % D55 = (q/144)*Aero_Ref_Area.*C_D_M55; %[lbf}
% % D6 = (q/144)*Aero_Ref_Area.*C_D_M6 %[lbf}
% % D65 = (q/144)*Aero_Ref_Area.*C_D_M65; %[lbf}
% % D7 = (q/144)*Aero_Ref_Area.*C_D_M7; %[lbf}
% % 
% % output_drag_table = [D4 D45 D5 D55 D6 D65 D7];
% % xlswrite('drag.xlsx',output_drag_table)