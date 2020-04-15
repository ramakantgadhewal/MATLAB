function [L] = cal_L(u, c, gam)
    coeff = (gam - 1) / c / c;
    L = [0.5*(0.5*coeff * u * u + u / c)     -0.5*(coeff * u + 1 / c)    0.5*coeff
        1 - 0.5*coeff * u * u                coeff * u                   -coeff
        0.5*(0.5*coeff * u * u - u / c)      -0.5*(coeff * u - 1 / c)    0.5*coeff];
end