%% Input Data for Example 2.2 
nsd = 2;	      % Number of space dimensions 
ndof= 2;     	  % Number of degrees-of-freedom per node
nnp = 3;    	  % Number of nodal points
nel = 2;     	  % Number of elements
nen = 2;     	  % Number of element nodes
nd  = 4;          % Number of prescribed displacement degrees-of-freedom

neq = ndof*nnp;	  % Number of equations
 
f 	= zeros(neq,1);   % Initialize force vector(local)
d 	= zeros(neq,1);   % Initialize displacement matrix(local)
K 	= zeros(neq);     % Initialize stiffness matrix(local)
 
% element properties
CArea 	= [1       1   ];   	% Elements area  
E     	= [1000    1000];   	% Young's Modulus
CArea_ex= 0;               % Normalizing coefficient for elements area
E_ex    = 0;               % Normalizing coefficient for Young's modulus
f_ex    = 0;               % Normalizing coefficient for external force
P_ex    = 0;               % Normalizing coefficient for the coordinate of node points
mag_ex  = 1;               % Magnification coefficient

% solver choosing
solver  = 2;           % 1 for reduction, while 2 for penalty.
crtcl_value = 1e6;     % convergence critical value
epsl    = 1e-15;       % convergence accuracy
pnl0    = 1e10;        % initial penalty factor
pnl     = 1e8;         % step-length of penalty factor

% prescribed displacements
% -node num- -node deg-
d0 = [1        1       % d0 is a local variable.
      1        2
      2        1
      2        2]';

% prescribed forces
f(5)	= 10;	   % External force

% output plots
plot_truss 	= 'yes';
plot_nod	= 'yes';
plot_disp   = 'yes';

% mesh generation and geometry calculation
truss_mesh_ex2_2;