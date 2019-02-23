function [ R_PAL ] = SM02_Palresist(R)
    n = length(R);
    R_PAL_0 = 0;
    for i = 1:n
        R_PAL_0 = R_PAL_0 + sum(1/R(i));
    end
    R_PAL = 1 / R_PAL_0;
end