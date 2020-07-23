function CFD01_heat2D
% 1D heat transfer equation solver using FDM
%% include the global flags and input the data
 include_flags_CFD01_heat2D;
 input_file_CFD01_heat2D;

%% discrete the solution domain
% the discrete variable
u  = zeros(Nx+1,Ny+1);
x  = zeros(1,Nx+1);
y  = zeros(1,Ny+1);
t  = zeros(1,Nt+1);
dx = (x2-x1)/Nx;
dy = (y2-y1)/Ny;
dt = min(sigma_1*dx^2,sigma_2*dy^2)/gamma_0;
                        % in order to satisfy the stability conditions
for i = 1:Nx+1
    x(i) = x1 + dx*(i-1);
    for j = 1:Ny+1
        y(j) = y1 + dy*(j-1);
        u(i,j) = IC(x(i),y(j));
    end
end
for i = 1:Nt+1
    t(i) = i*dt;
end
% the plot variable
T = zeros(1,Np);
for i = 2:Np
    T(i) = (i-1)*Nt/(Np-1);
end
T(1) = 0;
leng_T = length(T);
d_T = zeros(leng_T-1);
d_T(1) = T(1);
for i = 1:leng_T-1
    d_T(i+1) = T(i+1) - T(i);
end
%% iteration
tic
if ind_sol == 1         % FTCS
    for i = 1:leng_T
        for ii = 1:d_T(i)
            u0 = u;
            % calculate (i+1)-th step
            for k = 2:Ny
                for j = 2:Nx
                    u(j,k) = sigma_1 * u0(j,k+1) + sigma_2 * u0(j+1,k) ...
                           + (1 - 2*sigma_1 - 2*sigma_2) * u0(j,k) ...
                           + sigma_1 * u0(j,k-1) + sigma_2 * u0(j-1,k);
                end
            end
            % update the boundary conditions
            t_ord = ii + T(i-1);
            for m = 1:Nx+1
                u(m,1)    = BC(x(m),y1,t(t_ord));
                u(m,Ny+1) = BC(x(m),y2,t(t_ord));
            end
            for n = 1:Ny+1
                u(1,n)    = BC(x1,y(n),t(t_ord));
                u(Nx+1,n) = BC(x2,y(n),t(t_ord));
            end
        end
        % plot the interesting time points
        subplot(Np_h,Np_c,i)
        hold on;
        mesh(x,y,u');
    end
end
toc
%% plot the interesting time points
% if ind_plot == 1
%     mesh(x,y,u');
% elseif ind_plot == 2
%     surf(x,y,u');
%     shading interp;
% end
% if ind_view == 1
%     view(90,-90);
% end
%% write to the workspace
TIME_x = sigma_1*Nt/(Nx^2)
end

function [f] = IC(x,y)
% input the initial condition
f = 2*x*y;
if (x < 0.75) && (x > 0.25) && (y < 0.75) && (y > 0.25)
    f = 1;
end
end

function [f] = BC(x,y,t)
% input the boundary condition
% % steady
% if (x == 0) || (y == 0)
%     f = 0;
% elseif x == 1
%     f = 2*y;
% elseif y == 1
%     f = 2*x;
% end
% unsteady
if (x == 0) || (y == 0)
    f = 0;
elseif x == 1
    f = 2*y + sin(pi*80*t);
elseif y == 1
    f = 2*x + sin(pi*80*t);
end
end