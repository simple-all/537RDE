function pc = p_det(theta,prcj,pmin)
pc=prcj*pmin*exp(-(log(prcj))*(theta/(2*pi)));
end