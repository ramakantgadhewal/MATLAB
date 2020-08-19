function nls1_Stff(x0)
global a;
eps = 1e-7;
M = 1e+3;
n = 1;
x = x0;
a(1) = x;
y = (20- 2*x^2 - x^3)/10;
z = (20- 2*y^2 - y^3)/10;
x1 = x - (y - x)^2 / (z - 2*y + x);
while abs(x1 - x) > eps
    x = x1;
    y = (20- 2*x^2 - x^3)/10;
    z = (20- 2*y^2 - y^3)/10;
    x1 = x - (y - x)^2 / (z - 2*y + x);
    n = n + 1
    a(n) = x;
    a(n)
    if n > M
        error('迭代可能不收敛！');
    end
end
end

