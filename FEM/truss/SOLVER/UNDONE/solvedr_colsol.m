% partition and solve the system of equations [Using Reduction Technique~Active Column Solution]
function [d,f_E] = solvedr_colsol(K,f)
include_flags;

% partition the matrix K, vectors f and d
K_E	 = K(1:nd,1:nd);                     	  % Extract K_E matrix 
K_F	 = K(nd+1:neq,nd+1:neq);                  % Extract K_E matrix
K_EF = K(1:nd,nd+1:neq);                      % Extract K_EF matrix
d_E  = zeros(nd,1);                           % Extract d_E vector
f_F  = f(nd+1,1);                             % Extract f_F vector

% construct the skyline matrix m[r] of K_F
% //TODO:complete the function
m   = makeskyline(K_F);

% use colsol to solve for d_F
r   = neq - nd;                 % the actual degree of freedom of the system
R   = f_F;
[IERR, K_F, R] = COLSOL(r, m, K_F, R);
d_F = R;
 
% reconstruct the global displacement d
d = [d_E             
     d_F];                
 
% compute the reaction f_E
f_E = K_E*d_E+K_EF*d_F;
 
% write to the workspace
solution_vector_d	= d   .* 10^d_ex
reactions_vector 	= f_E .* 10^f_ex
