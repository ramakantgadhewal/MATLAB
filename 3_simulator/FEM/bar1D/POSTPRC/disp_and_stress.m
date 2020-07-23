% compute stresses and displacements in the current element  
function [err] = disp_and_stress(e,d,err) 
include_flags; 

de = d(LM(:,e));              % extract element nodal displacements 
IENe = IEN(:,e);              % extract element connectivity information 
xe = x(IENe);                 % extract element coordinates 
J = (xe(nen) - xe(1))/2;      % Jacobian  
[w , gp] = gauss(ngp);        % Gauss points and weights

% compute stresses at Gauss points 
for i = 1:ngp   
    xt  = 0.5*(xe(1)+xe(nen))+J*gp(i);     % Gauss point in the physical coordinates 
    gauss_pt(i) = xt;                      % store gauss point information   
         
    N  = Nmatrix1D(xt,xe);        % extract shape functions  
    B  = Bmatrix1D(xt,xe);        % extract derivative of shape functions  
  
    Ee = N*E(IENe);                  % Young's modulus at element Gauss points    
    stress_gauss(i)  = Ee*B*de;      % compute stresses at Gauss points 
end

% print stresses at element gauss points 
fprintf(1,'%d\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\n',e,gauss_pt(1),gauss_pt(2),stress_gauss(1),stress_gauss(2)); 
     
% compute displacements and stresses  
xplot  = linspace(xe(1),xe(nen),nplot);      % equally distributed coordinate within an element                   

for i = 1:nplot 
    xi   = xplot(i);               % current coordinate   
    N    = Nmatrix1D(xi,xe);       % shape functions  
    B    = Bmatrix1D(xi,xe);       % derivative of shape functions    

    Ee   = N*E(IENe);              % Young's modulus   
    displacement(i) = N*de;        % displacement output 
    stress(i)       = Ee*B*de;     % stress output 
end 

% plot displacements and stresses
figure(2) 
subplot(2,1,1);  
plot(xplot,displacement); hold on; 
ylabel('displacement');  title('FE analysis of 1D bar'); 

subplot(2,1,2); plot(xplot,stress); hold on; 
ylabel('stress'); xlabel('x');

% compute the error in the natural boundary (only valid if only one N.B.C.)
for i = 1:nen
    if flags(IEN(e,i)) == 1        % check if natural boundary
        xi0 = x(IEN(e,i));
        N0  = Nmatrix1D(xi0,xe);    % shape functions  
        B0  = Bmatrix1D(xi0,xe);    % derivative of shape functions    
        Ee0 = N0*E(IENe);          % Young's modulus
        stress0 = Ee0*B0*de;       % stress on the natural boundary
        err = err + abs(stress0 - n_bc(IEN(e,i)));
    end
end