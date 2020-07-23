function [b_opt,k0_opt] = calScheme_optTaylor(l,r,prL,prR,nu,e_res)
    syms k a y y0;
    Num = l + r - 3;
    y = sym(zeros(Num));
    y0 = sym(zeros(Num));
    for i = 1:Num
        y(i) = sin((-l+i-1)*k);
        y0(i) = taylor(y(i),k,'Order',Num);
    end
    % TODO: finish the construction of b_opt
    
end