%% Input file for CFD01_heat1D

% include the global flags
include_flags_CFD01_heat1D;

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