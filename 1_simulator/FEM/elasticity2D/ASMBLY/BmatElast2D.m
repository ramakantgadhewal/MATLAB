% B matrix function for 2D elasticity  
function [B, detJ] = BmatElast2D(eta,psi,C) 
    include_flags;
    if nen == 4
      % Calculate the Grad(N) matrix 
        GN    = 0.25 * [eta-1  1-eta   1+eta   -eta-1; 
                        psi-1  -psi-1  1+psi    1-psi]; 
        J     = GN*C;        % Compute Jacobian matrix  
         
        detJ  = det(J);     % Jacobian 
       
        BB     = J\GN;       % compute the derivative of the shape functions 
         
        B1x     = BB(1,1); 
        B2x     = BB(1,2); 
        B3x     = BB(1,3); 
        B4x     = BB(1,4); 
        B1y     = BB(2,1); 
        B2y     = BB(2,2); 
        B3y     = BB(2,3); 
        B4y     = BB(2,4); 
        
        B = [ B1x      0     B2x     0      B3x    0      B4x     0  ; 
                0     B1y     0     B2y      0     B3y     0      B4y;  
              B1y     B1x    B2y    B2x     B3y    B3x    B4y     B4x];
    elseif nen == 3
        Ae = det([1 C(1,1) C(2,1);
                  1 C(1,2) C(2,2);
                  1 C(1,3) C(2,3)]) * 0.5;
        x21 = C(1,2) - C(1,1);
        x32 = C(1,3) - C(1,2);
        x13 = C(1,1) - C(1,3);
        y23 = C(2,2) - C(2,3);
        y31 = C(2,3) - C(2,1);
        y12 = C(2,1) - C(2,2);
        B = 0.5/Ae* [y23  0  y31  0  y12  0;
                     0   x32  0  x13  0  x21;
                     x32 y23 x13 y31 x21 y12];
        detJ = 2*Ae;
    end
end