%% getLoD
% function that calcualtes the L/D given the following inputs...
% 1. Dynamic Pressure in Pounds per Square Foot (q) [psf]
% 2. Mach Number (M) [Nondim]
% 3. Aero Reference Area )Areo_Ref_Area [sq.in]

function [LoD, C_L, C_D]=getLoD(q,M,Aero_Ref_Area, weight_gross)

%% Lift, Drag, and L over D
% Drag Calculation
% Current Analysis assumes alpha = 0
% The drag coefficients were pulled from GHV and curve fit to allow for
% calculating the drag at different Mach numbers. I can expand this to use
% multiple alphas if we think that is neccesary. Data that was curve fit
% was for M = 4 to 6.
C_D = 0.00023*M^2 - 0.00282*M + 0.01869; %[Nondim]
D = (q/144)*Aero_Ref_Area.*C_D; %[pounds]

% Lift Calculation
% L = W; determine the lift coefficient required to support the weight at
% the specified dynamic pressure
C_L = weight_gross/((q/144)*Aero_Ref_Area); %[Nondim]
% Calculate L/D using C_L and C_D
LoD = C_L/C_D; %[Nondim]
end