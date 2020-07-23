function [x,y] = optimization(x,y)
    % calculate the source term
    cal_source;

    % iteration to refine the grids
    exam = 1;
    iter = 0;
    while exam == 1 && iter < 1e4
        [x,y,exam] = iteration(x,y);
        iter = iter + 1;
    end
end

function cal_source
    include_flags;
    Psour = zeros(2*m+num+1,n);
    Qsour = zeros(2*m+num+1,n);
    % NOTICE: eta-line[j=const,1:n] is defined along the airfoil,
    % - while ksi-line[i=const,1:2*m+num+1] is orthonogal to it.
    % Mesh width is required to be 1 along ksi&eta axis in compute coordinate;
    % - in other words, dksi = deta = 1.
    for i = 1:2*m+num+1
        for j = 1:n
        P = 0;
        Q = 0;
            for mm = 1:attr_num
                Psour(i,j) = P - attr(mm,5)*sign(i-attr(mm,2))*exp(-attr(mm,6)...
                    *sqrt(attr(mm,1)*(i-attr(mm,2))^2+attr(mm,3)*(j-attr(mm,4))^2));
                
                Qsour(i,j) = Q - attr(mm,5)*sign(j-attr(mm,4))*exp(-attr(mm,6)...
                    *sqrt(attr(mm,1)*(i-attr(mm,2))^2+attr(mm,3)*(j-attr(mm,4))^2));
            end
        end
    end
end