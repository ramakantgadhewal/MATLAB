%% Input Data for Problem 2.4
nsd = 2;	      % Number of space dimensions 
ndof= 2;     	  % Number of degrees-of-freedom per node
nnp = 4;    	  % Number of nodal points
nel = 3;     	  % Number of elements
nen = 2;     	  % Number of nodes per element
nd  = 5;          % Number of prescribed displacement degrees-of-freedom

neq = ndof*nnp;	  % Number of equations

f 	= zeros(neq,1);         % Initialize force vector(local)
d 	= zeros(neq,1);         % Initialize displacement matrix(local)
K 	= zeros(neq);           % Initialize stiffness matrix(local)

% element properties
CArea 	= [1   2   1];   	% Normalized elements area
E     	= [1   1   1];   	% Normalized Young's modulus
CArea_ex= -2;               % Normalizing coefficient for elements area
E_ex    = 11;               % Normalizing coefficient for Young's modulus
f_ex    = 3;                % Normalizing coefficient for external force
P_ex    = 0;                % Normalizing coefficient for the coordinate of node points
mag_ex  = 4;                % Magnification coefficient for visualization

% solver choosing
solver  = 2;           % 1 for reduction, while 2 for penalty.
crtcl_value = 1e6;     % convergence critical value
epsl    = 1e-4;       % convergence accuracy
pnl0    = 1e4;         % initial penalty factor
pnl     = 1e2;         % step-length of penalty factor

% prescribed displacements
% -node num(glo)- -node deg-
d0 = [1             1       % d0 is a local variable.
      1             2
      2             1
      2             2
      3             2]';

% prescribed forces
f(7)    = 1;        % Normalized force

% output plots
plot_truss 	= 'yes';
plot_nod	= 'yes';
plot_disp   = 'yes';

% mesh generation
truss_mesh_prb2_4;