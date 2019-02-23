function [ke, fe] = elast2Delem(e) 
include_flags; 
 
ke  = zeros(nen*ndof,nen*ndof);     % initialize element stiffness 
fe  = zeros(nen*ndof,1);            % initialize element force vector 
 
% get coordinates of element nodes  
je = IEN(:,e);   
C  = [x(je); y(je)]; 

if nen == 4
    [w,gp] = gauss(ngp);   % get gauss points and weights 
    C = C';
    % compute element stiffness matrix and element nodal force vector
    for i=1:ngp 
       for j=1:ngp 
           eta = gp(i);             
           psi = gp(j); 
           N             = NmatElast2D(eta,psi);       % shape functions matrix   
           [B, detJ]     = BmatElast2D(eta,psi,C);     % derivative of the shape functions
           ke = ke + w(i)*w(j)*B'*D*B*detJ;   % element conductance matrix 
           be = N*b(:,e);                     % interpolate body forces using element shape functions 
           fe = fe + w(i)*w(j)*N'*be*detJ;    % element nodal force vector  
       end        
    end
elseif nen == 3
    % compute element stiffness matrix and element nodal force vector
    [B,detJ] = BmatElast2D(0,0,C);     % derivative of the shape functions 
    ke = 0.5*detJ*B'*D*B;
        % assume that b is stored in (x1 y1 x2 y2 x3 y3)
    bx = b(1,e) + b(3,e) + b(5,e);
    by = b(2,e) + b(4,e) + b(6,e);
    be = detJ/24*[bx+b(1,e) by+b(2,e) bx+b(3,e) by+b(4,e) bx+b(5,e) by+b(6,e)]';
    fe = fe + be;
end
end