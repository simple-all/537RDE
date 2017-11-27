function [ D,dD ] = D_dD( M,gamma )
% D(M)
D=M/((1+((gamma-1)/2)*M^2)^((gamma+1)/(2*(gamma-1))));

% dD/dM = dD(M)
dD=(D/M)*(1-M^2)/(1+0.5*(gamma-1)*M^2);
end

