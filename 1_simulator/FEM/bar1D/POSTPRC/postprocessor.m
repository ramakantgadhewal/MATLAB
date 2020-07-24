% postprocessing
function postprocessor(d)
include_flags;
if strcmpi(plot_type,'conduct')==1;
    fprintf(1,'\n Print heat flux at the Gauss points \n')
    fprintf(1,'Element\t\t x(gauss1) \t\t x(gauss2) \t\t heat flux(gauss1) \t\t heat flux(gauss2)\n')
    fprintf(1,'--------------------------------------------------------------------------------- \n')
    err = 0;
    for e = 1:nel
       % compute heat flux(\nabla T) and temperature for the current element
       [err] = htflx_and_tmpr(e,d,err);
    end
    fprintf(1,'Compute the error in the natural boundary (only valid if only one N.B.C.)')
    % compute the error in the natural boundary (only valid if only one N.B.C.)
    n_bc_error = err
else
    fprintf(1,'\n Print stresses at the Gauss points \n')
    fprintf(1,'Element\t\t x(gauss1) \t\t x(gauss2) \t\t stress(gauss1) \t\t stress(gauss2)\n')
    fprintf(1,'--------------------------------------------------------------------------------- \n')
    err = 0;
    for e = 1:nel
       % compute stresses and displacements for the current element
       [err] = disp_and_stress(e,d,err);
    end
    fprintf(1,'Compute the error in the natural boundary (only valid if only one N.B.C.)')
    % compute the error in the natural boundary (only valid if only one N.B.C.)
    n_bc_error = err
end