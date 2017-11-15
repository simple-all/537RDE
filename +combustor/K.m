function k = K (gamma,R,pmax,pr,Tmax)
k = (pmax*2*pi*gamma/(R*Tmax*log(pr)))*(exp((log(pr))/gamma)-1);
end