function [u_plot] = cal_UprDownlF(u_plot,u,leng_T,d_T,l,r,b)
include_flags_CFD03_FDM;
% construct Q s.t. Q*u = \Sigma(b*u)
    I = eye(Nx+1);
    Q = zeros(Nx+1);
    for i = 1:l
        N0 = l + 1 - i;
        b0 = b(i);
        Q0 = b0*[I(1:Nx+1,N0+1:Nx+1) I(1:Nx+1,1:N0)];
        Q = Q + Q0;
    end
    Q = Q + b(l+1)*I;
    for i = 1:r
        N0 = Nx + 1 - i;
        b0 = b(l + i + 1);
        Q0 = b0*[I(1:Nx+1,N0+1:Nx+1) I(1:Nx+1,1:N0)];
        Q = Q + Q0;
    end
% the core iteration (4-order Runge-Kutta method)
    wait_flag = 0;
    hwait = waitbar(0,'Please wait...');
    for i = 1:leng_T
        for j = 1:d_T(i)
            u0 = u;
            u1 = u0 - CFL/4*Q * u0;
            u2 = u0 - CFL/3*Q * u1;
            u3 = u0 - CFL/2*Q * u2;
            u  = u0 - CFL  *Q * u3;
            wait_flag = wait_flag + 1/Nt;
            waitbar(wait_flag)
        end
        u_plot(:,i) = u(:,1);
    end
    close(hwait);
end