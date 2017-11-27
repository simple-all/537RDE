function cp=cp_air(T)
% Curve fit for Cp of the air
% Units T~[K]; cp~[J/(kg K)]
cp=(1.9327e-10)*(T^4)-(7.9999e-7)*(T^3)+(1.1407e-3)*(T^2)-(4.4890e-1)*T+1.0575e3;
end