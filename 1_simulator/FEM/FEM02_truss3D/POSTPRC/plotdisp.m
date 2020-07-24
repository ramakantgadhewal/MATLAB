function plotdisp(d);
%% Include global flags
include_flags;

%% compute the actual scale of coordinate
x_atl = x .* 10^P_ex;
y_atl = y .* 10^P_ex;
if nsd == 3
    z_atl = z .* 10^P_ex;
end

%% correct the scale of displacement for visualization
d_cor = d .* 10^(d_ex + mag_ex);

%% plot the deformed structure of the truss
% check if displacement plot is requested
if strcmpi(plot_disp,'yes')==1;
    if nsd == 1
        for i = 1:nnp
            x_atl(i) = x_atl(i) + d_cor((i-1)*ndof+1);   
            y_atl(i) = 0;
        end
    elseif nsd == 2
        for i = 1:nnp
            x_atl(i) = x_atl(i) + d_cor((i-1)*ndof+1);   
            y_atl(i) = y_atl(i) + d_cor((i-1)*ndof+2);
        end
    elseif nsd == 3
        for i = 1:nnp
            x_atl(i) = x_atl(i) + d_cor((i-1)*ndof+1);   
            y_atl(i) = y_atl(i) + d_cor((i-1)*ndof+2);
            z_atl(i) = z_atl(i) + d_cor((i-1)*ndof+3);
        end
    end
    for i = 1:nel
        XX = [x_atl(IEN(1,i)) x_atl(IEN(2,i)) ];
        YY = [y_atl(IEN(1,i)) y_atl(IEN(2,i)) ];
        if nsd == 3;
            ZZ = [z_atl(IEN(1,i)) z_atl(IEN(2,i))];
            plot3(XX,YY,ZZ,'r');hold on;
        else
            plot(XX,YY,'r');hold on;
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
    title('Deformed Truss Plot');
end