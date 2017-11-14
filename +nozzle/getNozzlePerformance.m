function [Me, Ve ] = getNozzlePerformance( D_outer, D_inner, GAMMA, Ts_throat, MW, Mt )
%getNozzlePerformance Uses basic isentropic flow equations to get the exit
%   performance of an aerospike nozzle 
%   
%   Assuming axisymmetric aerospike nozzle that converges to an
%       infinitesimally small final point.
%   Does NOT generate a contour for the aerospike nozzle, only gives final
%       performance.
%
% INPUTS:
%   D_outer: outer diameter [unit independent] of the combustor, essentially where the
%       top of the aerospike starts
%   D_inner: inner diameter [unit independent] of the combustor
%   GAMMA: [] burned value of the ratio of specific heats (products),
%       assumed constant
%   Ts_throat: [K] throat static temp, used for velocity 
%   MW: molecular weight of escaping gas 
%   Mt: mach number at the throat of supersonic engine
%
% OUPUTS:
%   Eps_max: [] maximum expansion ratio for the aerospike nozzle
%   M_exit: [] exit mach at the end of the aerospike
%   Ve: [m/s] gas velocity at the end of the aerospike
% 
% ASSUMPTIONS: 
%   Throat area is the area of the ring formed by the outer and inner
%       diameters (which are based on gap size)
%   Constant GAMMA for the calculations
%   Performance of the nozzle is based on the maximum expansion ratio
%   Isentropic flow across the nozzle

%% Calculate Nozzle Performance
Ru = 8314;

throat_area = pi * D_outer^2 /4 - pi * D_inner^2 /4;    %xx^2 units

Eps_max =   pi*D_outer^2/4 / throat_area;      


M_guess = 6;                                            %supersonic guess

func = @(Me) Me/Mt * ((1 + (GAMMA-1)/2*Mt^2)/(1 + (GAMMA-1)/2*Me^2))^...
    (1/(GAMMA-1)) * sqrt((1+(GAMMA-1)/2*Mt^2)/(1+(GAMMA-1)/2*Me^2)) -...
    1/Eps_max;

options = optimset('Display', 'off');
Me = fsolve(func, M_guess, options);


Ve = sqrt(2* GAMMA * (Ru/MW) * Ts_throat / (GAMMA-1) * ...
    (1 - 1/(1+(GAMMA+1)/2*Me^2)));                      %m/s escape vel



end

