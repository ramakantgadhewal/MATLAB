function [d] = transbkd(d1,d0);
% //TODO:complete the function
I = zeros(nd);  % global degrees of the prescribed displacements
for i = 1:nd
    I(i) = (d0(1,i)-1)*ndof + d0(2,i);
end