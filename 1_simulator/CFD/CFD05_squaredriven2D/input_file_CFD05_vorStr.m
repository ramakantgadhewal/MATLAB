%% Input file for CFD05_incompr2D_vorStr
%% include flags
include_flags_CFD05_vorStr;

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
outputName = 'CFDtest_Re10000_FTCS';
outputInt = 200;

% iteration parameters
nIter_max = 1e6;
err_uMag_thre = 1e-6;
err_psf_thre  = 1e-6;

% algorithm parameters
ind_ps2omg = 1;
                % 1: 2-order center scheme
                % 2: 3-order scheme
ind_uv2omg = 1;
                % 1: FTCS scheme
                % 2: 1-order upwind scheme
                % 3: 2-order upwind scheme
ind_psifun = 1;
                % 1: Gauss-Seidel scheme
ind_ps2uv = 1;
                % 1: 2-order center scheme
                % 2: 3-order scheme

% derived quantities
uMag = 1;          % non-dim
lMag = 1;          % non-dim
Lx = lMag;
Ly = lMag;
dx = Lx/Nx;
dy = Ly/Ny;
ds = min(dx,dy);
dt = cap * min(Re*ds*ds/2,2/Re);

%% Initialization
x = 0:dx:Lx;
y = 0:dy:Ly;
U = zeros(Nx+1,Ny+1);     % x-velocity
V = zeros(Nx+1,Ny+1);     % y-velocity
ix = 2:Nx;
iy = 2:Ny;
inx = 3:Nx-1;
iny = 3:Ny-1;

%% Boundary condition
U(2:Nx,Ny+1) = uMag;