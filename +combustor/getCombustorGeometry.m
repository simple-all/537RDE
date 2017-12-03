function [ lambda, D_outer, gap_width, D_inner,do_range, di_range, L] = getCombustorGeometry( fuel_cellsize )
%getCombustorGeometry Uses basic empirical relations for sizing an RDE
%combustor
%   Everything taken from Bykovskii paper:
%   Assume a cell size is an input given a certain fuel and pressure. any 
%       units inputted are what are outputted
%   Assuming the relations from Bykovskii for sizing (choosing the base
%       values without any +_ accounted for)

%% Calculate Chamber Sizes

lambda = 0.7 * fuel_cellsize;

D_outer = 40 * lambda;

gap_width = 0.2 * (12 * fuel_cellsize);         %minimum gap width required

D_inner = D_outer - 2*gap_width;                %inner diameter assuming the minimum gap width is all that's needed to run

% Full range of possible values

h = 12 * lambda;
K = [5 9];
do_range = (h * K) ./ pi;
delta = 0.2 * h;
di_range = do_range - 2 * delta;
L = [2, 4] * h;




end

