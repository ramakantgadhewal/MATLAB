function CFD01_heat1D
% 1D heat transfer equation solver using FDM
% input the data
x1 = 0;             % left endpoints of x
x2 = 1;             % right endpoints of x
gamma_0 = 1;
                    % heat transfer coefficient
sigma_0 = 1.0;
                    % the characteristic footstep:sigma = gamma *dt/dx2
Nx = 500;
                    % the number of x mesh points
T  = [0, 125, 250, 375, 500, 750, 1000];
                    % the ordinal number of interesting time points
                    % requirement: monotone increasing
Nt = max(T);        % the number of t mesh points
ind = 2;
                    % the index of the solver: 1 for FTCS, 2 for BTCS-dir

%% discrete the solution domain
% the discrete variable
u  = zeros(1,Nx+1);
x  = zeros(1,Nx+1);
t  = zeros(1,Nt+1);
dx = (x2-x1)/Nx;
dt = sigma_0*dx^2/gamma_0;
for i = 1:Nx+1
    x(i) = x1 + dx*(i-1);
    u(i) = IC(x(i));
end
for i = 1:Nt+1
    t(i) = i*dt;
end
% the plot variable
leng_T = length(T);
u_plot = zeros(Nx+1,leng_T);
d_T = zeros(leng_T-1);
d_T(1) = T(1);
for i = 1:leng_T-1
    d_T(i+1) = T(i+1) - T(i);
end
%% iteration
tic
if ind == 1         % FTCS
    u0 = zeros(1,Nx+1);
    for i = 1:leng_T
        for j = 1:d_T(i)        % the core iteration
            t_ord = j + T(i-1);
            [u(1,1),u(1,Nx+1)] = BC(t(t_ord));
            for k = 2:Nx
                u0(1,k) = sigma_0*u(1,k+1) + (1-2*sigma_0)*u(1,k) +sigma_0*u(1,k-1);
            end
            u(1,2:Nx) = u0(1,2:Nx);
            u0 = 0;
        end
        u_plot(:,i) = u(1,:);
    end
elseif ind == 2     % BTCS
    u0 = zeros(1,Nx+1);
    Q  = zeros(Nx+1,Nx+1);
    Q(1,1)       = 1;
    Q(Nx+1,Nx+1) = 1;
    for k = 2:Nx
        Q(k,k-1) = sigma_0;
        Q(k,k)   = - 1 - 2*sigma_0;
        Q(k,k+1) = sigma_0;
    end
    for i = 1:leng_T
        for j = 1:d_T(i)
            t_ord = j + T(i-1);
            u0 = -u;
            [u0(1,1),u0(1,Nx+1)] = BC(t(t_ord));
            u(1,:) = u0(1,:)/Q';
        end
        u_plot(:,i) = u(1,:);
    end
end
toc
%% plot the interesting time points
for i = 1:leng_T
    plot(x(:),u_plot(:,i));
    hold on;
end
%% write to the workspace
TIME = sigma_0/(Nx^2)*T;
end

function [f] = IC(x)
% input the initial condition
if x < 0.3
    f = 0;
elseif x < 0.6
    f = 1;
elseif x <= 1.0
    f = 1 + 2.5*(x-0.6);
end
end

function [bd_l,bd_r] = BC(t)
% input the boundary condition
% % steady
% bd_l = 0;
% bd_r = 2;
% unsteady
bd_l = 0;
bd_r = 2 + sin(pi*500*t);
end