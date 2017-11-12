
function coneang=find_cone_angle(M,thetas,g) % find_cone_angle.m;
% Taylor-Maccoll conical flow solution using equations given by
% Anderson, 1990, "Modern Compressible Flow (with Historical % Perspective)"
% 2nd edition, McGraw-Hill, p301-303. %
% M - Mach number upstream of shock
% thetas - shock wave angle (in degrees)
% g - ratio of specific heats %
% Limiting cone angle is set to 0.1 degrees %
% David BT Sercombe, 18 December 2003
% David R Buttsworth, 22 December 2003
global gamma
m1=M; % mach no.
gamma=g; % ratio of specific heats
% starting shock angle
thetas=thetas.*(pi)/180;
% stream deflection angle just behind the wave
delta=atan(2.*cot(thetas).*(((m1.^2).*(sin(thetas).^2)-1)./...
    ((m1.^2).*(gamma+cos(2*thetas))+2)));
% normal component of the free stream Mach No.
mn1=m1.*sin(thetas);
% normal component of the post shock Mach no.
mn2=sqrt(((mn1.^2)+(2/(gamma-1)))./((2*gamma./...
    (gamma-1)).*(mn1.^2)-1));
% calculation of post shock Mach no.
m2=mn2./sin(thetas-delta);
% initial dimensionless component of the post shock velocity
v_in=(2/((gamma-1).*(m2.^2))+1).^(-.5);
% radial and normal components of the velocity calculated % as boundary conditions
v_rin=v_in.*cos(thetas-delta); v_thin=v_in.*sin(thetas-delta);
% lower bound for cone angle (limiting case) 
endtheta=0.1.*(pi)/180;
% 'options' command switches on event detection (crossing point % detection) and refines the final solution by a factor of 4
options=odeset('Events','on','Refine',4);
% 'ode23' solver returns the numerical solution for the
% current shock angle
[theta,v]=ode23('conical.taymaceqn',[thetas, endtheta],[v_rin, v_thin], options);
% converts the angle values to degrees
theta2=theta.*(180/(pi));
% returns the cone angle from the numerical values found
coneang=theta2(length(theta2));
end
% EOF