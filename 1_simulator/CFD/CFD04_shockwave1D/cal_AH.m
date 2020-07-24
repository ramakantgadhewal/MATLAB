function [A] = cal_AH(u, H, gam)
    A = [0                     1                       0
        (gam-3)*u*u/2         (3-gam)*u                gam-1
        0.5*(gam-1)*u^3-u*H    H-(gam-1)*u*u           gam*u];
end