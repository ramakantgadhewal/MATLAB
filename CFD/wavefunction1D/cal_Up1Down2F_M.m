function [u_plot] = cal_Up1Down2F_M(u_plot,u,leng_T,d_T)
include_flags_CFD03_FDM;

    I  = eye(Nx+1);
    Q1 = -CFL*(CFL-1)*(CFL-2)*[I(1:Nx+1,2:Nx+1) I(1:Nx+1,1)];
    Q2 = 3*(CFL+1)*(CFL-1)*(CFL-2)*I;
    Q3 = -3*(CFL+1)*CFL*(CFL-2)*[I(1:Nx+1,Nx+1) I(1:Nx+1,1:Nx)];
    Q4 = (CFL+1)*CFL*(CFL-1)*[I(1:Nx+1,Nx:Nx+1) I(1:Nx+1,1:Nx-1)];
    Q  = 1/6*(Q1 + Q2 + Q3 + Q4);
    for i = 1:leng_T
        for j = 1:d_T(i)        % the core iteration
            u = Q * u;
        end
        u_plot(:,i) = u(:,1);
    end
end