clear;
clc;
close all;

Skin = [212.3 504.7 893.6 1352.9 1486.2];
Flaps = [47.5 129.3 242.4 366.1 514.9];
Tails = [17.5 47.1 87.1 131.3 183.8];
TPS = [33.7 88.5 161.9 242.4 337.9];
Spars = [15.8 31 46.6 61.6 77.2];
Inlet = [28.2 54.4 80.2 105.8 131.4];
Isolator = [37.7 74.1 110.2 146.2 182.2];
Nozzle = [119.5 223.2 344 441.7 554.5];
Ballast = [70 175 300 400 500];

flowRate = [1 2 3 4 5];


ppSkin = polyfit(flowRate, Skin, 2);
ppFlaps = polyfit(flowRate, Flaps, 2);
ppTails = polyfit(flowRate, Tails, 2);
ppTPS = polyfit(flowRate, TPS, 2);
ppSpars = polyfit(flowRate, Spars, 2);
ppInlet = polyfit(flowRate, Inlet, 2);
ppIsolator = polyfit(flowRate, Isolator, 2);
ppNozzle = polyfit(flowRate, Nozzle, 2);
ppBalast = polyfit(flowRate, Ballast, 2);

mdot = 0.1:0.1:6;
skinW = polyval(ppSkin, mdot);
flapsW = polyval(ppFlaps, mdot);
tailsW = polyval(ppTails, mdot);
tpsW = polyval(ppTPS, mdot);
sparsW = polyval(ppSpars, mdot);
refArea = 0.77;

% figure;
% hold on;
% 
% plot(flowRate, Skin, 'r*');
% plot(flowRate, Flaps, 'g*');
% plot(flowRate, Tails, 'b*');
% plot(flowRate, TPS, 'k*');
% plot(flowRate, Spars, 'c*');
% 
% plot(mdot, skinW, 'r-');
% plot(mdot, flapsW, 'g-');
% plot(mdot, tailsW, 'b-');
% plot(mdot, tpsW, 'k-');
% plot(mdot, sparsW, 'c-');

%legend('Skin', 'Flaps', 'Tails', 'TPS', 'Spars');

skinWeight = polyval(ppSkin, refArea);
flapWeight = polyval(ppFlaps, refArea);
tailWeight = polyval(ppTails, refArea);
tpsWeight = polyval(ppTPS, refArea);
sparsWeight = polyval(ppSpars, refArea);
aeroTotalWeight = skinWeight + flapWeight + tailWeight + tpsWeight + sparsWeight;

fuel = 21 * 2.2;
hp = 1.33;

Avionics = 294;

P_operating = 14.7;
Volume_tank = 2257.8;
tankWeight = (300*P_operating*Volume_tank)/(1*10^6);


isRho = 4506; % Density of titanium [kg/m^3]


% Cone
coneRad = 0.113;
coneLength = 30.25 * 0.0254;
SACone = pi * coneRad * (coneRad + sqrt(coneLength^2 + coneRad^2));
coneWall = 0.005;
VCone = SACone * coneWall;
coneWeight = VCone * isRho * 2.2;

% Isolator
isOd = 0.24;
isId = 0.113;
length = 45.34 * 0.0254 + 0.1;
isWall = 0.005;
isOuterV = pi * (((isOd + isWall * 2) / 2)^2 - (isOd / 2)^2) * length;
isInnerV = pi * ((isId / 2)^2 - ((isId + isWall * 2) / 2)^2) * length;
isWeight = (isInnerV + isOuterV * isRho) * 2.2; % lb

% Combustor
cRho = 7940; % kg/m^3
cOd = 0.2035;
cId = 0.16;
cLength = 10 * 0.0254;
tWall = 0.05 * 0.0254;
tCoolant = 0.125 * 0.0254;
tOuterWall = 0.005;

AOuter = pi * (((cOd + tWall * 2) / 2)^2 - (cOd / 2)^2) +...
    pi * 0.5 * (((cOd + (tWall + tCoolant) * 2)/2)^2 - ((cOd + tWall * 2) / 2)^2) + ...
    pi * (((cOd + (tWall + tCoolant + tOuterWall) * 2) / 2)^2 - ((cOd + (tWall + tCoolant) * 2)/2)^2);

AInner = pi * ((cId / 2)^2 - ((cId + tWall * 2) / 2)^2) +...
    pi * 0.5 * (((cId + tWall * 2) / 2)^2 - ((cId - (tWall + tCoolant) * 2)/2)^2) + ...
    pi * (((cId - (tWall + tCoolant) * 2)/2)^2 - ((cId - (tWall + tCoolant + tOuterWall / 3) * 2) / 2)^2);

VComb = (AOuter + AInner) * cLength;
combWeight = VComb * cRho * 2.2; % lb

% Nozzle
nozzleWeight = 1.55 * 2.2;

% Ballast
ballastWeight = polyval(ppBalast, refArea);

% Turbopump
turboWeight = 11.7;

totalWeight =aeroTotalWeight + fuel + Avionics + tankWeight + coneWeight + isWeight + nozzleWeight + combWeight + ballastWeight + turboWeight;

w = [100 300 600];
M = [9.4 7.8 5.6];

ppMach = polyfit(w, M, 1);