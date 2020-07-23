function [x,y,exam] = iteration(x,y)
 % There are two ADDitional options:
 % - 1.to determine whether to iterate the trailing line or not
 % - 2.to determine whether to iterate the trailing area or not
 % - To make the grids more smooth, it is better to iterate over the whole
 % - - domain.
    include_flags;

    a = x;      % temp x
    b = y;      % temp y
    err = 0;    % iteration error

    % NOTICE: eta-line[j=const,1:n] is defined along the airfoil,
    % - while ksi-line[i=const,1:2*m+num+1] is orthonogal to it.
    % Mesh width is required to be 1 along ksi&eta axis in compute coordinate;
    % - in other words, dksi = deta = 1.
    if iterflag_mode == 1
        LeftEndPoint = 2;
        RightEndPoint = 2*m+num;
    else
        LeftEndPoint = m+1;
        RightEndPoint = m+num+1;
    end
    for i = LeftEndPoint:RightEndPoint
        for j=2:n
            % approximation of derivatives
            xksi = (x(i+1,j)-x(i-1,j))/2;
            xeta = (x(i,j+1)-x(i,j-1))/2;
            yksi = (y(i+1,j)-y(i-1,j))/2;
            yeta = (y(i,j+1)-y(i,j-1))/2;
            xksieta = (x(i+1,j+1)-x(i+1,j-1)-x(i-1,j+1)+x(i-1,j-1))/4;
            yksieta = (y(i+1,j+1)-y(i+1,j-1)-y(i-1,j+1)+y(i-1,j-1))/4;
            tempAlpha = xeta^2 + yeta^2;
            tempBeta  = xksi*xeta + yksi*yeta;
            tempGamma = xksi^2 + yksi^2;
            Jacobi  = xksi*yeta - yksi*xeta;
            % construct the difference scheme of Laplace equation
            bwe  = tempAlpha;
            bsn  = tempGamma;
            bp   = 2*(bwe+bsn);
            cpx  = - 2*tempBeta * xksieta;
            cpy  = - 2*tempBeta * yksieta;
            Pbar = Jacobi^2*(Psour(i,j)*xksi + Qsour(i,j)*xeta);
            Qbar = Jacobi^2*(Psour(i,j)*yksi + Qsour(i,j)*yeta);
            a(i,j) = (bwe*x(i-1,j)+bwe*x(i+1,j)+bsn*x(i,j-1)+bsn*x(i,j+1)+cpx+Pbar)/bp;
            b(i,j) = (bwe*y(i-1,j)+bwe*y(i+1,j)+bsn*y(i,j-1)+bsn*y(i,j+1)+cpy+Qbar)/bp;
            dist = (a(i,j)-x(i,j))^2+(b(i,j)-y(i,j)^2);
            % find the maximum of square error among elements
            if dist > err
                err = dist;
            end
        end
    end
    if iterflag_glueline == 1
        for i = 2:m
            iacross = 2*m+num+2-i;
            % approximation of derivatives
            xksi = (x(i+1,1)-x(i-1,1))/2;
            xeta = (x(i,2)-x(iacross,2))/2;
            yksi = (y(i+1,1)-y(i-1,1))/2;
            yeta = (y(i,2)-y(iacross,2))/2;
            xksieta = (x(i+1,2)-x(iacross-1,2)-x(i-1,2)+x(iacross+1,2))/4;
            yksieta = (y(i+1,2)-y(iacross-1,2)-y(i-1,2)+y(iacross+1,2))/4;
            tempAlpha = xeta^2 + yeta^2;
            tempBeta  = xksi*xeta + yksi*yeta;
            tempGamma = xksi^2 + yksi^2;
            Jacobi  = xksi*yeta - yksi*xeta;
            % construct the difference scheme of Laplace equation
            bwe  = tempAlpha;
            bsn  = tempGamma;
            bp   = 2*(bwe+bsn);
            cpx  = - 2*tempBeta * xksieta;
            cpy  = - 2*tempBeta * yksieta;
            Pbar = Jacobi^2*(Psour(i,j)*xksi + Qsour(i,j)*xeta);
            Qbar = Jacobi^2*(Psour(i,j)*yksi + Qsour(i,j)*yeta);
            a(i,1) = (bwe*x(i-1,1)+bwe*x(i+1,1)+bsn*x(iacross,2)+bsn*x(i,2)+cpx+Pbar)/bp;
            b(i,1) = (bwe*y(i-1,1)+bwe*y(i+1,1)+bsn*y(iacross,2)+bsn*y(i,2)+cpy+Qbar)/bp;
            a(iacross,1) = a(i,1);
            b(iacross,1) = b(i,1);
            dist = (a(i,1)-x(i,1))^2+(b(i,1)-y(i,1)^2);
            % find the maximum of square error among elements
            if dist > err
                err = dist;
            end
        end
    end
    
    if err > qual
        exam = 1;
    else
        exam = 0;
    end
    x = a;      % return new x
    y = b;      % return new y
end