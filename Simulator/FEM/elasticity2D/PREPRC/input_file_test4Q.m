% Input Data for STAPPpp test4Q 
include_flags;

%% material properties 
E  = 1000;     % Young's modulus  
ne = 0.3;      % Poisson ratio   
D  = E/(1-ne^2) * [1    ne     0           
                   ne    1     0   
                   0     0     (1-ne)/2]; 

%% mesh specifications 
nsd  = 2;         % number of space dimensions 
nnp  = 8;         % number of nodal nodes 
nel  = 5;         % number of elements 
nen  = 4;         % number of element nodes
nee  = 4;         % number of element edges
ndof = 2;         % degrees of freedom per node 
neq  = nnp*ndof;  % number of equations 

f = zeros(neq,1);      % initialize nodal force vector 
d = zeros(neq,1);      % initialize nodal displacement matrix 
K = zeros(neq);        % initialize stiffness matrix 

counter    = zeros(nnp,1);  % counter of nodes for stress plots 
nodestress = zeros(nnp,3);  % stresses at nodes for the stress plots [sxx syy sxy] 

P     = zeros(neq,1);          % point forces applied at nodes 
b     = zeros(nen*ndof,nel);   % body forces defined at nodes 
P(1) = -2;
P(3) = 3;
P(5) = 2;
P(7) = -3;
%% boundary conditions and force array
% array to set B.C flags 
% dof:    1x    1y    2x     2y    3x    3y    4x     4y  
flags = zeros(16,1);
flags(1) = 2;
flags(2) = 2;
flags(4) = 2;
flags(7) = 2;
nd    = 4;             % number of essential boundary conditions (x and y)
ngp   = 2;
% essential B.C array
e_bc  = zeros(neq,1);
% natural B.C. array  
nbe   = 0;             %  number of edges on the boundary

%% plot flag
plot_mesh      = 'no'; 
plot_nod       = 'no'; 
plot_disp      = 'no'; 
compute_stress = 'yes'; 
plot_stress_xx = 'no'; 
plot_mises     = 'no'; 
fact           = 10;      % factor for scaled displacements plot  

%% mesh generation 
% node:  1    2    3    4    5    6    7    8
x   =  [0.0  2.0  2.0  0.0  0.4  1.4  1.5  0.3];     % X coordinate 
y   =  [0.0  0.0  3.0  2.0  0.4  0.6  2.0  1.6];     % Y coordinate 
% the order is required to be clockwise
IEN =  [1  2  6  5;
        2  3  7  6;
        3  4  8  7;
        4  1  5  8
        5  6  7  8]';     % connectivity array 
 
% function to plot the mesh 
plotmesh; 