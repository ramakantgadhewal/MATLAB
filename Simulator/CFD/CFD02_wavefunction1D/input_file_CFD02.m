%% Input file for CFD02

% include the global flags
 include_flags_CFD02

% input the general data
x1 = 0;             % left endpoints of x
x2 = 1;             % right endpoints of x
a = 1;
                    % wave velocity
CFL = 0.1;
                    % CFL number
Nx = 256;
                    % the number of x mesh points
dx = (x2-x1)/Nx;
                    % calculation of dx
dt = CFL*dx/a;
                    % calculation of dt
TIME  = [0 0.1];
                    % the interesting physical time points
                    % requirement: monotone increasing, integer times of dt
T  = TIME/dt;
                    % the ordinal number of interesting time points
                    % requirement: monotone increasing
Nt = max(T);        % the number of the mesh timepoints

% input the initial value parameters
% initialCondition = 1:
epsil = 0.1;        % small number for initial condition
Nphi = Nx/4;        % the number of different wavenumbers
phik = rand(Nphi,1);% random phase for each single wave with different wavenumbers
k0 = 5;             % parameter indicates those waves with high wavenumbers

% input the solver indicator
initialCondition_flag = 2;
timepointsPlot_flag = 1;
calError_flag = 0;
l = 4;
r = 3;
prL = 0;
prR = 1;
e_res = 0.05;
                    % parameter controlling the resolution ratio
nu = 8;
                    % parameter controlling the weight when optimizing
mv = 0.004;
                    % parameter controlling the movement due to dispersion
ind = [12];
name_ind = {'1-order Upwind'; 'Lax-Wendroff'; 'Warming-Beam'; 'Up2Down1'; 'Up1Down2F\_M';...
    ''; ''; ''; ''; '';...
    'Original UprDownlF'; 'Optimized Integration UprDownlF'; 'Optimized Taylor-expansion UprDownlF'};
                    % the index of the differential scheme (periodic boundary condition):
                    %   1 for Up1Down0F scheme
                    %   2 for Lax-Wendroff scheme
                    %   3 for Warming-Beam scheme
                    %   4 for Up2Down1 scheme
                    %   5 for Up2Down1F_M scheme
                    %   11 for UplDownrF_origin scheme
                    %   12 for UplDownrF_opt scheme (using Integration method)
                    %   13 for UplDownrF_opt scheme (using Taylor expansion method) (unfinished)
