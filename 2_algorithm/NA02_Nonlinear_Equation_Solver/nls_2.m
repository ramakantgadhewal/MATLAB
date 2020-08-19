function nls_2(x0)
eps = 1e-7;
M = 1e+3;
n = 1;
x = x0;
a(1) = x;
x1 = nthroot(20- 10*x - 2*x^2,3);
while abs(x1 - x) > eps
    x = x1;
    x1 = nthroot(20- 10*x - 2*x^2,3); 
    n = n + 1
    a(n) = x;
    a(n)
    if n > M
        error('迭代可能不收敛！');
    end
end
end

