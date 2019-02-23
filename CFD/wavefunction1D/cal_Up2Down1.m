function [u_plot] = cal_Up2Down1(u_plot,u,leng_T,d_T)
include_flags_CFD03_FDM;

    I  = eye(Nx+1);
    Q1 = (2-CFL)/3*[I(1:Nx+1,2:Nx+1) I(1:Nx+1,1)];
    Q2 = (1+CFL)/3*[I(1:Nx+1,3:Nx+1) I(1:Nx+1,1:2)];
    Q  = Q1 + Q2;
    for i = 1:leng_T
        for j = 1:d_T(i)        % the core iteration
            u = Q * u;
        end
        u_plot(:,i) = u(:,1);
    end
end