% partition and solve the system of equations [Using Penalty Method]
% NO requirement to the index of the prescribed displacements
function [d,f_E] = solvedr_pnl(K,f,d0);
%% include global flags
include_flags;

%% initialization
err   = 0;          % error of prescribed displacements
c_ind = 0;          % cycle index
I     = zeros(nd);  % global degrees of the prescribed displacements

%% loop of penalty for d
K1 = K;             % temporary global stiffness matrix
for i = 1:nd
    I(i) = (d0(1,i)-1)*ndof + d0(2,i);
end
for i = 1:nd
    K1(I(i),I(i)) = pnl0;
end
d = K1\f;
for i = 1:nd
    err = err + d(I(i))^2;
end
K1_0 = K1;          % initial penalty matrix (for output)
% falling into the loop
while err >= epsl
    err   = 0;      % initial the error
    for i = 1:nd
        K1(I(i),I(i)) = K1(I(i),I(i)) + pnl;
    end
    d = K1\f;
    for i = 1:nd
        err = err + d(I(i))^2;
    end
    c_ind = c_ind + 1;
    while c_ind >= crtcl_value
        error('The loop of penalty may not converge!');
    end
end

%% compute the reaction f_E
for i = 1:nd
    f_E(i,1) = K(I(i),:)*d(:);
end

%% write to the workspace
cond_K1_0           = cond(K1_0)
cond_K1             = cond(K1)
cycle_index         = c_ind
solution_vector_d	= d   .* 10^d_ex
reactions_vector 	= f_E .* 10^f_ex
