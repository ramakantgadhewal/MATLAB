% generate the element stiffness matrix for each element
function [ke, be] = trusselem(e);
include_flags;

const = CArea(e)*E(e)*leng(e);                  % constant coefficient for each truss element

if nsd == 1
    be(e,:) = [-1, 1]/leng(e);
    ke = const * (be(e,:))' * be(e,:);                    % stiffness matrix for 1D
elseif nsd == 2
    c  = phi(1,e);
    s  = phi(2,e);
    be(e,:) = [-c,-s, c, s]/leng(e);
    ke = const * (be(e,:))' * be(e,:);                    % stiffness matrix for 2D
elseif nsd == 3
    c1 = phi(1,e);
    c2 = phi(2,e);
    c3 = phi(3,e);
    be(e,:) = [-c1,-c2,-c3, c1, c2, c3]/leng(e);
    ke = const * (be(e,:))' * be(e,:);                    % stiffness matrix for 3D
end
