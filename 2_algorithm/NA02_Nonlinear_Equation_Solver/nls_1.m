function nls_1(x0)
eps = 1e-7;
M = 300;
n = 1;
x = x0;
a(1) = x;
x1 = (20- 2*x^2 - x^3)/10;
while abs(x1 - x) > eps
    x = x1;
    x1 = (20- 2*x^2 - x^3)/10; 
    n = n + 1
    a(n) = x;
    a(n)
    if n > M
        error('迭代可能不收敛！');
    end
end
end


