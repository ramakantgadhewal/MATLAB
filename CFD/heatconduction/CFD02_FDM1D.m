function CFD02_FDM1D
% 1D heat transfer equation solver using FDM
%% include the global flags and input the data
 include_flags_CFD02_FDM1D;
 input_file_CFD02_FDM1D;

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
TIME = sigma_0/(Nx^2)*T
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

function x = followup(A,b)
% Thomas solver of the linear algebraic equation Ax = b,
%       where Q is a tri-diagonal matrix.
%   Sometimes, the scale of A is too big for internal 
%       storage. Under this circumstance, A is required
%       to be stored as a built-in matrix.
n = rank(A);
for i = 1:n
    if(A(i,i)==0)
        disp('Error: There is a diagonal element valued zero!');
        return;
    end
end
d = ones(n,1);
a = ones(n-1,1);
c = ones(n-1);
for i = 1:n-1
    a(i,1) = A(i+1,i  );
    c(i,1) = A(i  ,i+1);
    d(i,1) = A(i  ,i  );
end
d(n,1) = A(n,n);
% solve Ly = b
for i = 2:n
    ad = a(i-1,1)/d(i-1,1);
    d(i,1) = d(i,1) - ad*c(i-1,1);
    b(i,1) = b(i,1) - ad*b(i-1,1);
end
% solve Ux = y
x(n,1) = b(n,1)/d(n,1);
for i = (n-1):(-1):1
    x(i,1) = (b(i,1)-c(i,1)*x(i+1,1))/d(i,1);
end
end