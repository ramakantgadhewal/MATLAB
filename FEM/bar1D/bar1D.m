%%%%%%%%%%%%%%%%%%
% 1D FEM Program (Chapter 5) %
% Haim Waisman, Rensselaer %
%%%%%%%%%%%%%%%%%%
clear all;
close all;

% include global variables
include_flags;

% Preprocessing
[K,f,d] = preprocessor;

% Element matrix computations and assembly
for e = 1:nel
    [ke,fe] = barelem(e);
    [K, f] = assembly(K,f,e,ke,fe);
end

% Add nodal boundary force vector
f = naturalBC(f);

% Partition and solution
[d,f_E] = solvedr(K,f,d);

% Postprocessing and calculate the N.B.C.error
postprocessor(d);

% plot the exact solution
if plot_exact == 1
    ExactSolution_conduct;
else
    ExactSolution;
end