function [R] = cal_R(u, c, H)
    R = [1 1 1
        u-c u u+c
        H-u*c 0.5*u*u H+u*c];
end