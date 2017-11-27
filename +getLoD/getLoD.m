%% getLoD
% function that calcualtes the L/D given the following inputs...
% 1. Dynamic Pressure in Pounds per Square Foot (q)
% 2. Mach Number (M)
% 3. Wing span (b)
% 4. Root Chord (c_r)
% 5. Tip Chord (c_t)
% 6. Gross Weight (weight_gross)

function [LoD, C_L, C_D, Aero_Ref_Area]=getLoD(q,M,b,c_r,c_t,weight_gross)
%% Aero Plan form
s = b/2; %[in] half span
Aero_Ref_Area = (c_t + c_r)*s; %[sq. inches]
sweep = atand((c_r - c_t)/s); %[deg] wing sweep angle

%% Lift, Drag, and L over D
% Drag Calculation
% Current Analysis assumes alpha = 0
% The drag coefficients were pulled from GHV and curve fit to allow for
% calculating the drag at different Mach numbers. I can expand this to use
% multiple alphas if we think that is neccesary. Data that was curve fit
% was for M = 4 to 6.
C_D = 0.00023*M^2 - 0.00282*M + 0.01869;
D = (q/144)*Aero_Ref_Area.*C_D;

% Lift Calculation
% L = W; determine the lift coefficient required to support the weight at
% the specified dynamic pressure
C_L = weight_gross/((q/144)*Aero_Ref_Area);
% Calculate L/D using C_L and C_D
LoD = C_L/C_D;
end