function thetas=find_cone_shock_angle(M,thetac,g);
% find_cone_shock_angle(M,thetac,g);
%
% Taylor-Maccoll conical flow solution using equations given
% Anderson, 1990, "Modern Compressible Flow (with Historical % Perspective)"
% 2nd edition, McGraw-Hill, p301-303. %
% M - Mach number upstream of shock
% thetac - cone angle (in degrees)
% g - ratio of specific heats %
% David BT Sercombe and David R Buttsworth, 22 December 2003
% Function iterates through a series of shock angles using a
% iteration technique until the found cone angle = the known
% sets the maximum number of iterations
nitermax=50;
% sets the iteration tolerance
tol=1e-9;
% returns the accuracy (initial value for the ?while? loop)
accuracy=tol*2;
% initial value for the ?while
niter=0;
% sets the maximum shock angle
thetasmax=60;
% finds the corresponding cone
coneangmax=conical.find_cone_angle(M,thetasmax,g);
% sets the minimum shock angle
thetasmin=asin(1/M)*180/pi+0.2;
% finds the corresponding cone angle
coneangmin=conical.find_cone_angle(M,thetasmin,g);
% ?while? loop continues the iteration until the solution % converges to within the tolerance or exceeds the maximum
% number of iterations allowable
while (niter<nitermax)&&(accuracy>tol)
    % finds the mean of the shock angle
    thetasmean=(thetasmax+thetasmin)/2;
    % finds the corresponding cone angle
    coneangmean=conical.find_cone_angle(M,thetasmean,g);
    % ?if? command determines whether the found value % forms the new upper or lower bound
    if coneangmean>thetac
        thetasmax=thetasmean;
    else
        thetasmin=thetasmean;
    end
    % advances the iteration counter
    niter=niter+1;
    % determines the new difference between the % found cone angle and the known cone angle
    accuracy=abs(coneangmean-thetac)/thetac;
end
% returns the found conical shock angle
thetas=thetasmean;
% EOF