function nls3_Newton(x0)
eps = 1e-7;
M = 1e+3;
n = 1;
x = x0;
a(1) = x;
x1 = x - (x^3 + 2*x^2 + 10*x -20)/(3*x^2 + 4*x +10);
while abs(x1 - x) > eps
    x = x1;
    x1 = x - (x^3 + 2*x^2 + 10*x -20)/(3*x^2 + 4*x +10); 
    n = n + 1
    a(n) = x;
    a(n)
    if n > M
        error('迭代可能不收敛！');
    end
end
end

