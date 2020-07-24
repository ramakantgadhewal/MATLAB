%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Truss (Chapter 2)        %
% Haim Waisman, Rensselaer    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all; 
 
% include global variables
include_flags;  

% Preprocessor Phase 
[K,f,d0]	= preprocessor;

% Calculation and assembly of element matrices
for e = 1:nel
    ke      = trusselem(e);
    K       = assembly(K,e,ke);
end

% Solution Phase
if solver == 1
    [d,f_E]	= solvedr(K,f);
elseif solver == 2
    [d,f_E]	= solvedr_pnl(K,f,d0);
elseif solver == 3
    [d,f_E] = solvedr_colsol(K,f);
end

% Postprocessor Phase 
postprocessor(d);