function [u_plot] = cal_Up1Down0(u_plot,u,leng_T,d_T)
include_flags_CFD03_FDM;

    I = eye(Nx+1);
    Q1 = (1-CFL)*I;
    Q2 = CFL*[I(1:Nx+1,2:Nx+1) I(1:Nx+1,1)];
    Q  = Q1 + Q2;
    for i = 1:leng_T
        for j = 1:d_T(i)        % the core iteration
            u = Q * u;
        end
        u_plot(:,i) = u(:,1);
    end
end