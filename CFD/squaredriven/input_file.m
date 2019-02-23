%% % Input file for CFD05squaredriven2D

%% include flags
include_flags;

%% domain definition
x1 = 0;
x2 = 1;
y1 = 0;
y2 = 1;
Nx = 128;
Ny = 128;
dx = (x2-x1)/Nx;
dy = (y2-y1)/Ny;
CFL = 0.1;
                % CFL number determining the time-propelling step length
dt = CFL * min(dx,dy);
                % uv_max = 1;
Re = 1000;
                % Re number of the problem
                % Re := U*L/nu
                %   L = 0.1 m, nu = 1e-6 m2/s.

%% Initialization
x = x1:dx:x2;
y = y1:dy:y2;
u = zeros(Nx+1,Nx+1);   % x-velocity
v = zeros(Nx+1,Nx+1);   % y-velocity
w = zeros(Nx+1,Nx+1);   % vorticity
pf = zeros(Nx+1,Nx+1);  % stream function

%% Boundary condition
u(:,1)    = 0;
u(2:Nx,Nx+1) = 1;
u(1,:)    = 0;
u(Nx+1,:) = 0;
v(:,1)    = 0;
v(:,Nx+1) = 0;
v(1,:)    = 0;
v(Nx+1,:) = 0;
pf(:,1)    = 0;
pf(:,Nx+1) = 0;
pf(1,:)    = 0;
pf(Nx+1,:) = 0;

%% Iteration and p-calculation parameters
err_steady_sup = 1.17e-1;
err_steady = 0;
err_pfConv_sup = 1e-6;
err_pfConv = 0;
cal_p_flag = 0;
UVmax_flag = 0;
ind = 3;
                % ind = 1: 1-order upwind scheme
                % ind = 2: 2-order upwind scheme
                % ind = 3: FTCS scheme