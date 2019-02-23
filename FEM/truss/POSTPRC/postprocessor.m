% Postprocessing function
function postprocesser(d);
include_flags;

% prints the element number and corresponding stresses
 fprintf(1,'element\t\t\tstress\n');
% Compute stress for each element
 stress = zeros(nel);
 for e=1:nel
     de        = d(LM(:,e));             % nodal displacements for each element
     const     = E(e);                   % constant coefficient for each element
     stress(e) = const*be(e,:)*de;       % computes stress(using strain matrix)
     stress(e) = stress(e) * 10^(f_ex - CArea_ex);
     fprintf(1,'%d\t\t\t%f\n',e,stress(e));
 end
% plot the deformed structure of the truss
 plotdisp(d);