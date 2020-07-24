%% Input file for CFD01_heat2D

% include the global flags
 include_flags_CFD01_heat2D;

 % input the data
x1 = 0;             % left endpoints of x
x2 = 1;             % right endpoints of x
y1 = 0;             % left endpoints of y
y2 = 1;             % right endpoints of y
gamma_0 = 1;
                    % heat transfer coefficient
sigma_1 = 0.25;
                    % the characteristic footstep:sigma = gamma *dt/dx2
sigma_2 = 0.25;
                    % the characteristic footstep:sigma = gamma *dt/dy2
Nx = 100;
                    % the number of x mesh points
Ny = 100;
                    % the number of y mesh points
Nt = 1000;        
                    % the number of t mesh points
Np  = 9;
Np_c = 3;
Np_h = 3;
                    % the number of subplots: (Np-1)|Nt && Np = Np_h * Np_c
                    % Generally, Np = length(T), where T is the time points
                    %   we are interested in.
ind_sol = 1;
                    % choose the solver: 1 for FTCS
% ind_plot = 1;
%                     % choose the way of plot:
%                     %       1:      mesh
%                     %       2:      surf_interp
% ind_view = 0;
%                     % choose the angle of view:
%                     %       0:      (-37.5,30) (default)
%                     %       1:      (90,-90)