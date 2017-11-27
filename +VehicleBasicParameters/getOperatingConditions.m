%% getOperatingConditions
% Finds operating conditions based on dynamic pressure (q) [psf] and Mach
% Mumber (M)

function [Mach,altitude,P0,T0,rho0,v0] = getOperatingConditions(q,M)
%% Find operating conditions
gamma = 1.4; % Ratio of specific heats
R = 1716; %Gas Constant [ft.lbf/R.slug]
h_max = 2400; %[100*ft]
for m = 1:length(M) %for each M vary altitude to get rho_error less than error
    for j = 1:h_max %Altitude iteration
        error = 0.001;
        h(j) = 100*j; %increase altitude by 100 feet for each iteration
        [t(j),p_psf(j),r(j),height(j)]=atmosphere(h(j),0);
        a(j) = sqrt(gamma*R*t(j)); %[ft/s]
        v(j) = M(m) * a(j); %[ft/s]
        rho(j) = (2*q) /(v(j)^2);
        rho_error(j) = (r(j)-rho(j))/r(j);
            if rho_error(j) < error
               Mach(m) = M(m);
               altitude(m) = h(j); %[ft]
               P0(m) = p_psf(j); %{psf]
               T0(m) = t(j); %[R]
               v0(m) = v(j); %[ft/s]
               rho0(m) = rho(j); %[slug/ft^3]
               break
            end
    end
end
end