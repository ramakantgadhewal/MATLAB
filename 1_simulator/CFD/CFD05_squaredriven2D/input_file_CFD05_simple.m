%% Input file for CFD05_incompr2D_simple
%% include flags
include_flags_CFD05_simple;

%% domain definition
% - cfl = uMag*dt/ds = dt/ds
% - sgm = dt/Re/ds/ds
% - ds = min(dx,dy)
%  from sgm < 1/2 ==> dt < Re*ds*ds/2 ... (1)
%  from cfl^2 < 2*sgm ==> dt < 2/Re   ... (2)
% (1),(2) ==> dt = cap*min(Re*ds*ds/2,2/Re), cap<1
Re = 10000;
cap = 2;
Nx = 128;
Ny = 128;

% output parameters
outputName = 'CFDtest_Re10000_SIMPLE';
outputInt = 200;

% iteration parameters
nIter_max = 1e6;
err_pres = 0;
err_pres_thre  = 1e-6;
err_uv_thre = 1e-6;
err_dpres_thre  = 1e-6;

% algorithm parameters
ind_uiter = 1;
                % 1: 1-order upwind scheme
                % 2: 2-order upwind scheme
ind_piter = 1;
                % 1: Gauss-Seidel scheme

% derived quantities
uMag = 1;          % non-dim
lMag = 1;          % non-dim
Lx = lMag;
Ly = lMag;
dx = Lx/Nx;
dy = Ly/Ny;
ds = min(dx,dy);
dt = cap * min(Re*ds*ds/2,2/Re);

%% For post-processing
x_p = (dx:dx:Lx)-dx/2;
y_p = (dy:dy:Ly)-dy/2;
x_u = 0:dx:Lx;
y_u = 0:dy:Ly;

%% Initialization
%  y
%  |
%   -- A，---
%  |         |
% B，  O，  ，D
%  |         |
%  0 - C， -- --- x
%
%  O: p(1,1)
%  A: v(1,2)
%  C: v(1,1)
%  B: u(1,1)
%  D: u(2,1)
%
% physical paras
u = zeros(Nx+1,Ny);
v = zeros(Nx,Ny+1);
p = zeros(Nx,Ny);
% Numerical paras
dp = zeros(Nx,Ny);
a = zeros(Nx-1,Ny);
a_up = zeros(Nx,Ny);
a_un = zeros(Nx,Ny);
a_vp = zeros(Nx-1,Ny+1);
a_vn = zeros(Nx-1,Ny+1);
b = zeros(Nx,Ny-1);
b_vp = zeros(Nx,Ny);
b_vn = zeros(Nx,Ny);
b_up = zeros(Nx+1,Ny-1);
b_un = zeros(Nx+1,Ny-1);
c = zeros(Nx+1,Ny+1);

%% Boundary condition
% u(:,Ny+1/2) = uMag;
u(1,:) = 0;
u(Nx+1,:) = 0;
v(:,1) = 0;
v(:,Ny+1) = 0;