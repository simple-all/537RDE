
function [value, isterminal, direction] = taymaceqn(t,g,flag); % taymaceqn.m
% Taylor-Maccoll conical flow solution using equations given by
% Anderson, 1990, "Modern Compressible Flow (with Historical
% Perspective)"
% 2nd edition, McGraw-Hill, p301-303.
%
% Function contains the differential equations in state-variable
% format for use with the 'ode23' numerical differential equation % solver.
% David BT Sercombe and David R Buttsworth, 20 December 2003
% specification of global variables
global gamma
% state variable form of the differential equation governing % conical flow
gdot(1)=g(2);
a=(gamma-1)/2; c=a.*(2*((g(1)).^3-g(1))+(g(2)-g(2).*g(1).^2-g(2).^3)...
    .*cot(t)-2*g(1).*g(2).^2)-(g(1).^2.*g(2).^2); d=a.*((g(2)).^2+(g(1)).^2 - 1) + (g(2)).^2; gdot(2)=c./d;
% events detection flagging routine set up to detect when g(2) % (angular velocity) is zero, stop the iteration from proceeding % and provide an estimation of where the crossing point occurs
if (nargin<3) || isempty(flag);
    value=[gdot(1);gdot(2)];
elseif flag == 'events';
    value=g; % values to check
    isterminal=[0;1]; % terminal on g(2)
    direction=[0;-1]; % detects falling slope
else
    error(['unknown flag''', flag]); % error flag
end
end% EOF