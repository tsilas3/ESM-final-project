
dz = 0.1; %depth grid spacing (m)
zf = 1; % lake depth in meters
zs = 0:dz:zf; % vector of depths in the lake

dt = 1; % time grid spacing = 1 day
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

M(1,1) = 1-C_D;
M(1,2) = 0;
M(end,end-1)=C_D;
M(end,end) = 1-C_D;


% Pre-allocate

T_all = nan(length(zs), length(ts));

T = 8 .* ones(length(zs), 1); %initial T throughout lake is annual avg 8C

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
