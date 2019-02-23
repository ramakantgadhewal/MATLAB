%% Input Data from Figure 2.8
nsd	= 1;     % Number of space dimensions 
ndof= 1;     % Number of degrees-of-freedom per node
nnp	= 3;     % Total number of global nodes
nel	= 2;     % Total number of elements
nen	= 2;     % Number of nodes in each element
nd	= 1;     % Number of prescribed displacement degrees-of-freedom
 
neq	= ndof*nnp;	% Number of equations
 
f	= zeros(neq,1);	% Initialize force vector
d	= zeros(neq,1);	% Initialize displacement vector
K	= zeros(neq); 	% Initialize stiffness matrix
 
% element properties 
CArea	= [.5       1];   	% Elements cross-sectional area  
E       = [1        1];   	% Young's Modulus
CArea_ex= 0;                % Normalizing coefficient for elements area
E_ex    = 0;                % Normalizing coefficient for Young's modulus
f_ex    = 0;                % Normalizing coefficient for external force
P_ex    = 0;                % Normalizing coefficient for the coordinate of node points
mag_ex  = 1;                % Magnification coefficient

% solver choosing
solver  = 2;           % 1 for reduction, while 2 for penalty.
crtcl_value = 1e6;     % convergence critical value
epsl    = 1e-15;       % convergence accuracy
pnl0    = 1e10;        % initial penalty factor
pnl     = 1e8;         % step-length of penalty factor

% prescribed displacements
% -node num- -node deg-
d0 = [1        1]'; % d0 is a local variable.

% prescribed forces
f(3)	= 1;      	% Force at node 3 in the x-direction

% output plots
plot_truss 	= 'yes';
plot_nod	= 'yes';
plot_disp   = 'yes';
 
% mesh generation and geometry calculation
truss_mesh_ex2_8;