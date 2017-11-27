%% getRange
% Calculates range



function [Range]=getRange(rho0,Aero_Ref_Area,TSFC, C_L, C_D, weight_gross, weight_zerofuel)
R_ft = 2*sqrt(2./(rho0*Aero_Ref_Area))*(1/(TSFC/3600)).*(sqrt(C_L)./C_D)*(sqrt(weight_gross)-sqrt(weight_zerofuel));
Range = R_ft./6076;
end