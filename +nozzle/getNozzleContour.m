function [ x,r ] = getNozzleContour( gamma, P_alt, P_chamber, D_outer, D_inner, Ts_throat, MW, Mt) 
%getNozzleContour Creates an ideal aerospike nozzle given certain inputs
%   
%% Original Author
% David Imbaratto, Cal Poly
% September 04, 2008
% This program creates an aerospike profile based on fluid properties % and
% design NPR using a Prandtl-Meyer expansion technique.

%% EDITS
% Drew Sherman, Purdue University
% 11/14/2017
% Functionalized the code to make a more general application that could be
%   looped and used to design for a variety of inputs. 
%
% INPUTS
%   gamma: [] specific heat ratio of combusted products
%   P_alt: [Pa] atmospheric pressure that the nozzle will expand to 
%   P_chamber: [Pa] average chamber pressure of the engine 
%   D_outer: [mm] outer diameter of the combustor, essentially where the
%       top of the aerospike starts
%   D_inner: [mm] inner diameter of the combustor
%   Ts_throat: [K] throat static temp, used for velocity 
%   MW: molecular weight of escaping gas 
%   Mt: mach number at the throat of supersonic engine

% OUTPUTS
%   x: [mm] length along the aerospike nozzle
%   r: [mm] radius of the aerospike nozzle along xS


%% Calculations
% Designate an NPR (Pc/Pe)
NPR = P_chamber/P_alt;

% Finding Nozzle Performance Parameters
[Me, Ve] = getNozzlePerformance(D_outer, D_inner, gamma, Ts_throat, MW, Mt);
M = linspace(1, Me, 100);

P_c_P_x = (1+(gamma-1)/2.*M.^2).^(gamma/(gamma-1));

P_x = NPR*P_alt ./ P_c_P_x;                         %[Pa]

% A_ratio = A / At
A_ratio = 1./M.*(2/(gamma+1)*(1+(gamma-1)/2.*M.^2)).^((gamma+1)/(2*(gamma-1)));

% Designate an exit radius
% Use r_e = 1 for nozzle scaled to unity
r_e = D_outer/2;                                    %[mm]
A_e = pi() * r_e^2;                                 %[mm^2]
A_t = A_e / A_ratio(end);                           %[mm^2]

% Prandtl-Meyer function
K = ((gamma+1)/(gamma-1))^(1/2);
v = K .* atan((M.^2-1).^(1/2)./K)-atan((M.^2-1).^(1/2));
theta_t = v(end);
r_t = (r_e^2 - A_t*cos(theta_t)/pi())^(1/2);
theta = theta_t - v;
mu = asin(1./M);

% Radial Coordinate of the spike contour
r = real((r_e^2 - (r_e^2 - r_t^2).*A_ratio.*sin(mu +theta)./(sin(mu).*cos(theta_t))).^(1/2));

% Axial Coordinate of the spike contour
x = real((r_e - r)./(tan(mu + theta)));

%% Plotting

% Plot the spike contour
cont_fig = figure;
figure(cont_fig);
hold on;
plot(x,r);
title(['Aerospike Contour for NPR = ', num2str(NPR), ', \gamma = ' num2str(gamma)]);
xlabel('x [mm]');
ylabel('r [mm]');

% Plot the nozzle throat
plot([x(1) 0],[r_t r_e],'r');
box on;
grid on;
axis equal;
hold off;

% Plot the area of spike contour
area_fig = figure;
hold on;
area(x,r);

% Plot the cowl
area([0 x(1)],[r_e (r(2)-r(1))/(x(2)-x(1))*x(1)+r_e], 'BaseValue', (r_e + ...
    r_e- r_t)*1.4);
title(['Aerospike Area Contour for NPR = ', num2str(NPR), ', \gamma = ' num2str(gamma)]);
colormap summer;
xlabel('x / r_e');
ylabel('r / r_e');
box on;
grid on;
axis equal;
alpha(0.5);
hold off;

% Plot 3-D Nozzle
ThreeD_fig = figure;
phi = linspace(0, 2*pi(), 100);
s_y = zeros(size(r,2), length(phi));
s_x = zeros(size(r,2), length(phi));
x_x = zeros(size(r,2), length(phi));
for k = 1:size(r,2);
    for j = 1:100;
        s_y(j,k) = r(k).*sin(phi(j));
        s_x(j,k) = r(k).*cos(phi(j));
        x_x(j,k) = x(k);
    end;
end;

surf(x_x, s_y, s_x,'FaceColor','interp','EdgeColor','none',...
 'FaceLighting','phong');
title(['Aerospike 3-D Contour for NPR = ', num2str(NPR), ', \gamma = ' num2str(gamma)]);
camlight right;
colormap pink;
axis equal;
xlabel('x');
ylabel('r');
zlabel('r');
view([23,19]);


end

