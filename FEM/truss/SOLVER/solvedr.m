% partition and solve the system of equations [Using Reduction Technique]
% ONLY effective when the prescribed displacements take the first few global degrees of freedom 
function [d,f_E] = solvedr(K,f)
include_flags;

% partition the matrix K, vectors f and d
K_E	 = K(1:nd,1:nd);                     	  % Extract K_E matrix 
K_F	 = K(nd+1:neq,nd+1:neq);                  % Extract K_E matrix
K_EF = K(1:nd,nd+1:neq);                      % Extract K_EF matrix
f_F  = f(nd+1:neq);                           % Extract f_F vector
d_E  = zeros(nd,1);                           % Extract d_E vector

% solve for d_F
d_F	=K_F\( f_F - K_EF'* d_E);
 
% reconstruct the global displacement d
d = [d_E             
     d_F];                
 
% compute the reaction f_E
f_E = K_E*d_E+K_EF*d_F;
 
% write to the workspace
cond_K_F            = cond(K_F)
solution_vector_d	= d   .* 10^d_ex
reactions_vector 	= f_E .* 10^f_ex