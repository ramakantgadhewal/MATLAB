%% Input Data for Problem 2.7
%% for 2D case
nsd = 2;	      % Number of space dimensions 
ndof= 2;     	  % Number of degrees-of-freedom per node
nnp = 12;    	  % Number of nodal points
nel = 22;       % Number of elements
nen = 2;     	  % Number of nodes per element
nd  = 4;        % Number of prescribed displacement degrees-of-freedom

%% for 3D case
% nsd = 3;	      % Number of space dimensions 
% ndof= 3;     	  % Number of degrees-of-freedom per node
% nnp = 12;    	  % Number of nodal points
% nel = 22;         % Number of elements
% nen = 2;     	  % Number of nodes per element
% nd  = 16;          % Number of prescribed displacement degrees-of-freedom

%% BOTH are OK
neq = ndof*nnp;	  % Number of equations

f 	= zeros(neq,1);   % Initialize force vector(local)
d 	= zeros(neq,1);   % Initialize displacement matrix(local)
K 	= zeros(neq);     % Initialize stiffness matrix(local)

% element properties
CArea   = ones(1,22);      	% Normalized elements area
E       = ones(1,22) .* 1.5;% Normalized Young's modulus
CArea_ex= -2;               % Normalizing coefficient for elements area
E_ex    = 11;               % Normalizing coefficient for Young's modulus
f_ex    = 3;                % Normalizing coefficient for external force
P_ex    = 0;                % Normalizing coefficient for the coordinate of node points
mag_ex  = 3.5;              % Magnification coefficient for visualization

% solver choosing
solver  = 2;           % 1 for reduction, while 2 for penalty.
crtcl_value = 1e6;     % convergence critical value
epsl    = 1e-4;        % convergence accuracy
pnl0    = 1e4;         % initial penalty factor
pnl     = 1e2;         % step-length of penalty factor

%% for 2D case
% prescribed displacements
% -node num(glo)- -node deg-
d0 = [1             1       % d0 is a local variable.
      1             2
      11            1
      11            2]';
% prescribed forces
f(10) = -7;                 % Normalized force

%% for 3D case
% % prescribed displacements
% % -node num(glo)- -node deg-
% d0 = [1             1       % d0 is a local variable.
%       1             2
%       1             3
%       2             3
%       3             3
%       4             3
%       5             3
%       6             3
%       7             3
%       8             3
%       9             3
%       10            3
%       11            1
%       11            2
%       11            3
%       12            3]';
% % prescribed forces
% f(14) = -7;              % Normalized force

%% BOTH are OK
% output plots
plot_truss 	= 'yes';
plot_nod	= 'yes';
plot_disp   = 'yes';

% mesh generation
truss_mesh_prb2_7;