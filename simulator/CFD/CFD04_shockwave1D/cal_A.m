function [A] = cal_A(u, E, gam)
    A = [0                     1                       0
        (gam-3)*u*u/2         (3-gam)*u                gam-1
        (gam-1)*u^3-gam*u*E   -1.5*(gam-1)*u*u+gam*E   gam*u];
end