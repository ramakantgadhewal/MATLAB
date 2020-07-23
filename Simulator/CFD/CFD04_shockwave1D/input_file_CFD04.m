%% input the domain of the 1D shockwave problem
include_flags_CFD04;

%% Flags marking the discrete scheme
%   Rusanov:            1
%   Jameson:            2
%   1-order FVS:        3
%   2-order FVS:        4
%   2-order FVS(Runge): 5
%   1-order Roe(FDM):   6
%   1-order Roe(FVM):   7
%   2-order TVD(FVM):   8
scheme_flag = 8;

%% define the domain constants
xleft = -0.5;   % the left endpoint of the domain
xright = 0.5;   % the right endpoint of the domain
T = 0.25;       % the end of the time calculation
Time = 0.20;    % the timepoint we are interested in
Nx = 3200;       % number of grid across x-axis
Nt = 4000;      % number of grid across t-axis
gam = 1.4;      % adiabatic exponent
fixEnt = 1e-1;  % entropy fixing factor
% calculate the grid size
dt = T/Nt;
dx = (xright - xleft)/Nx;

%% Initialize variables
rho  = zeros(Nx+1, Nt+1);
m    = zeros(Nx+1, Nt+1);
epsi = zeros(Nx+1, Nt+1);
u    = zeros(Nx+1, Nt+1);
p    = zeros(Nx+1, Nt+1);

% initial condition
rho(1:Nx/2+1,1) = 1;
rho(Nx/2+2:Nx+1,1) = 0.125;
u(1:Nx/2+1,1) = 0.75;
u(Nx/2+2:Nx+1,1) = 0;
p(1:Nx/2+1,1) = 1;
p(Nx/2+2:Nx+1,1) = 0.1;

% boundary condition
rho(1,:) = 1;
rho(Nx+1,:) = 0.125;
u(1,:) = 0.75;
u(Nx+1,:) = 0;
p(1,:) = 1;
p(Nx+1,:) = 0.1;

% initialize the conservative vector (rho, m, epsi)
m(:,:) = rho(:,:).*u(:,:);
epsi(:,:) = p(:,:)./(gam - 1) + 0.5.*rho(:,:).*u(:,:).^2;