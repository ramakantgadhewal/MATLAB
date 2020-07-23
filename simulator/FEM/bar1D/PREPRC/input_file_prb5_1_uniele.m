% Input Data for Problem 5.2(uniform elements)

nsd = 1;  % number of space dimensions
ndof = 1; % number of degrees-of-freedom per node
nnp = 3;  % number of nodal points
nel = 1;  % number of elements
nen = 3;  % number of element nodes

neq = ndof*nnp;    % number of equations

f = zeros(neq,1);  % initialize nodal force vector
d = zeros(neq,1);  % initialize nodal displacement vector
K = zeros(neq);    % initialize stiffness matrix

flags = zeros(neq,1);  % initialize flag vector
e_bc = zeros(neq,1);   % initialize vector of essential boundary condition
n_bc = zeros(neq,1);   % initialize vector of natural boundary condition

% element and material data (given at the element nodes)
E = 5*ones(nnp,1);       % nodal values Young's modulus
body = 100*ones(nnp,1);    % nodal values body forces
CArea = ones(nnp,1); % nodal values of cross-sectional area

% gauss integration
ngp = 2;  % number of gauss points

% flags(i) = 1 or 2  : node i is located on the natural or essential boundary
% essential boundary conditions
flags(1) = 2;  % flags to mark nodes located on the essential boundary
e_bc(1) = 0;   % value of essential B.C
nd = 1;        % number of nodes on the essential boundary

% natural boundary conditions
flags(nnp) = 1; % flags to mark nodes located on the natural boundary
n_bc(nnp) = 0;  % value of natural B.C

% point forces
P = 0;  % array of point forces
xp = 0;  % array of coordinates where point forces are applied
np = 0;  % number of point forces

% output plots
plot_bar = 'no';
plot_nod = 'no';
plot_type = 'conduct';  % choose the type of plot
plot_exact = 1;        % choose the exact solution to plot
nplot = nnp*10;         % number of points in the element to plot displacements and stresses

% mesh generation
bar_mesh_prb5_1_uniele;