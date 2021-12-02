clear all
close all

% 
% % constants
% A = 2.331e7; %set surface area of lake equal to 2.331e7 m^2
% u_2 = 5; %set 2 meter wind equal to 5 m/s
% phi = deg2rad(43.2); % latitude
% g = -9.8; % gravity, m/s^2
% 
% % Eddy Diffusivity
% 
% k = 0.4; %von Karman constant
% wstar = (1.2e-3) * u_2; %surface value of friction velocity
% p_0 = 1; %Prandtl number
% kstar = 6.6*(sin(phi))^0.5 * u_2^-1.84; % parameter of Ekman profile
% 
% rho = 1000; % kg/m^3, density of water
% % rho = @(Tk) (1 - 1.9549e-5*(abs(Tk-277))^1.68)*10^3;
% 
% drhodz = 0;
% %drhodz = ?
% 
% N = 1; %assume this is constant to simplify
% %N = @(Tk) (-g/rho(Tk))*drhodz;
% 
% Ri = @(z) (-1+((1+40*N^2*k^2*z^2/(wstar^2*exp(-2*kstar*z))))^0.5)/20;
% 
% %Kzt = @(z) (k*wstar*z/p_0)*exp(-kstar*z)*(1+37*Ri(z)^2)^(-1);
% Kzt = 1; %make this constant
% 
% % Heat source term
% 
% beta = 0.4; %shortwave radiation absorbed in surface layer
% 
% K_d = 1366; % W/m^2 incoming global shortwave radiation
% % used solar constant; this differs from the paper
% 
% A_s = 0.3;% shortwave albedo of water surface
% % this is an assumption we made that differs from the paper
% 
% 
% Kstar = K_d * (1-A_s); % net shortwave radiation at the surface
% 
% eta = @(z) 1.1925*z^-0.424; % light extinction coefficient
% % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6559290/
% 
% 
% Phi = @(z) (1-beta) * Kstar * exp(-eta(z)*z);
% 
% 
% % Set up math for surface boundary condition
% 
% %epsilon_a = 
% 
% sigma = 5.67e-8; % W m^-2 K^-4; Stefan-Boltzmann constant
% %T_a = 
% 
% L_d = epsilon_a * sigma * (T_a+273)^4;
% 
% A_1 = 0.3; % longwave albedo of lake surface
% L_dstar = (1-A_1)*L_d;
% 
% epsilon = 0.97;
% L_u = epsilon*sigma*(T_s+273)^4;
% 
% %E = ?
% 
% %T_s = %what is this temperature? how do we get it?
% L_e = 1.91846e6*((T_s+273)/(T_s+329.09))^2; %latent heat of vaporization
% Q_e = 1*L_e*E;
% 
% 
% R = gamma*((T_s-T_a)/e_0 - e_a);
% Q_h = R*Q_e;

%% Forward Euler Diffusion

% Code goes here

dz = 0.1; %depth grid spacing (m)
zf = 1; % lake depth in meters
zs = 0:dz:zf; % vector of depths in the lake

dt = 1; % time grid spacing = 1 second
tf = 3600; % 1 hr
ts = 0:dt:tf;

D = dz^2/4*dt; %hard code

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

% Pre-allocate

T_all = nan(length(zs), length(ts));

T = 10 .* ones(length(zs), 1);
T(1) = 15;
T(end) = 5;

T_all(:,1) = T;

% Forward Euler Diffusion

for k = 1:length(ts)-1
    Tnew = M*T;
    
    T_all(:,k+1) = Tnew;
    T = Tnew;
end

% plot

figure(1);
plot(zs, T_all(:,length(ts))) % end
