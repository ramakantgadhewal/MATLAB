function plottruss;
include_flags;

% compute the actual scale of coordinate
x_atl = x .* 10^P_ex;
y_atl = y .* 10^P_ex;
if nsd == 3
    z_atl = z .* 10^P_ex;
end

% check if truss plot is requested
if strcmpi(plot_truss,'yes')==1;  
    for i = 1:nel
        XX = [x_atl(IEN(1,i)) x_atl(IEN(2,i)) ];
        YY = [y_atl(IEN(1,i)) y_atl(IEN(2,i)) ];
        if nsd == 3;
            ZZ = [z_atl(IEN(1,i)) z_atl(IEN(2,i))];
            plot3(XX,YY,ZZ,'b');hold on;
        else
            plot(XX,YY,'b');hold on;
        end

        % check if node numbering is requested
        if strcmpi(plot_nod,'yes')==1; 
            if nsd == 3;
                text(XX(1),YY(1),ZZ(1),sprintf('%0.5g',IEN(1,i)));
                text(XX(2),YY(2),ZZ(2),sprintf('%0.5g',IEN(2,i)));
            else
                text(XX(1),YY(1),sprintf('%0.5g',IEN(1,i)));
                text(XX(2),YY(2),sprintf('%0.5g',IEN(2,i)));
            end
        end
    end
    title('Truss Plot');
end

% print mesh parameters
fprintf(1,'\tTruss Params \n');
fprintf(1,'------------------- \n');
fprintf(1,'No. of Elements  %d \n',nel);
fprintf(1,'No. of Nodes     %d \n',nnp);
fprintf(1,'No. of Equations %d \n\n',neq);
