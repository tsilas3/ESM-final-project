clear all
close all

%% Forward Euler Diffusion

% Varying diffusivity stuff
% Equation caluculated in Google Sheets from values in R. Mills (1973)

kappa_m = @(T_k) (16.4 - 0.151*T_k + 3.48e-4*T_k^2)*1e-9;

%% Parameters and constants

dz = 0.1; %depth grid spacing (m)
zf = 1; % lake depth in meters
zs = 0:dz:zf; % vector of depths in the lake

dt = 86400; % time grid spacing = 1 day
tf = 3600*24*365*3; % 3 years
ts = 0:dt:tf;

conductivity = 0.5918; % W m^-1 K^-1 - slightly less than the value at 10C

% Net incoming shortwave radiation:

Kstar = @(t) 500+100*sin(2*pi*t/(365*86400));
% assume all radiation absorbed is right at the surface: e^(eta*z)=0

% Net longwave radiation:
epsilon = 0.97;
sigma = 5.67e-8;
% assume longwave radiation out with constant temp 8C
L_u = epsilon*sigma* (8 + 273.15)^4;

% Assume flux of sensible and latent heat are zero
Q_e = 0;
Q_h = 0;

%Total heat balance
heat = @(t) Kstar(t) - L_u + Q_e + Q_h;

%% Pre-allocate

T = 8 .* ones(length(zs), 1); %initial T throughout lake is annual avg 8C

T_all = nan(length(zs), length(ts));
T_all(:,1) = T;

known_vec = zeros(length(zs), 1);

%% Loop

for k = 1:length(ts)-1

    D = kappa_m(T(1));

    C_D = D*dt/(dz)^2;

    % Forward Euler method

    M = sparse(length(zs),length(zs));

    % create matrix

    for i = 1:length(zs)
        for j = 1:length(zs)
            if i==j
                M(i,j) = 1-2*C_D;
            elseif i-1==j
                M(i,j) = C_D;
            elseif i+1==j
                M(i,j) = C_D;
            end
        end
    end

    M(1,1) = 1-C_D;
    M(1,2) = C_D;
    M(end,end-1)=C_D;
    M(end,end) = 1-C_D;

    % Forward Euler Diffusion

    
    F_sfc = heat(ts(k))/conductivity;
    
    known_vec(1) = -C_D*dz*F_sfc;
    
    Tnew = M*T + known_vec;
    
    T_all(:,k+1) = Tnew;
    T = Tnew;
end

% plot

figure(7);
plot(zs, T_all(:,length(ts))) % end
title('Temperature profile of the whole lake with depth at final time')
xlabel('Depth')
ylabel('Temperature (C)')

figure(8)
plot(ts/86400/365, T_all(1, :))
title('Temperature profile of the lake surface with time')
xlabel('Time (year)')
ylabel('Temperature (C)')

figure(9)
plot(ts/86400/365, T_all(end, :))
title('Temperature profile of the lake bottom with time')
xlabel('Time (year)')
ylabel('Temperature (C)')
