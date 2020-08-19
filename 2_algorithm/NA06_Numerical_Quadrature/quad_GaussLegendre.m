function [I,I0,err,exp] = quad_GaussLegendre(f,a,b,n,l)
% 功能：利用n点Gauss-Legendre求积公式的l次复化形式计算任意给定函数在给定区间上的积分，并给出误差
% 输入变量：被积函数f，被积区间[a,b]，节点数量n(默认为5)，复化个数l(默认为4)
% 输出变量：数值积分值I，精确积分值I0，积分计算误差err.[程序运行时间"tic-toc"]
% 变量设置说明：
%   [1]符号变量：
%       函数f(x)的自变量     x
%       替换后的变量         t
%       变量替换式x(t)       xt
%       替换后的函数f(x(t))  ft
%   [2]数值变量：
%       Legendre多项式零点   x0
%       P'_n在零点处的取值   y1
%       P_n-1在零点处的取值  y2
%       求积系数             A
%       子区间端点          [a1,b1]
%       f(x)在零点处的取值   f_x0
%       数值积分值           I
%       精确积分值           I0
%       积分计算误差         err
tic
if nargin == 4
    l = 4;
elseif nargin == 3
    n = 5;
    l = 4;
end
% 预分配内存
A = zeros(n,1);
f_x0 = zeros(n,1);
I = 0;
% 计算Gauss点和求积系数
[x0,y1,y2] = Gauss(n);
for j = 1:1:n
   A(j) = 2 /(n * y1(j) * y2(j));
end
% 复化n点Gauss-Legendre方法求积
for i = 1:1:l
    % 将[a,b]分为l个子区间
    a1 = a + (i-1)/l * (b-a);
    b1 = a + i    /l * (b-a);
    % 通过变量替换，将积分区间[a,b]标准化为[-1,1]
    syms x xt t;
    x = symvar(sym(f));
    k = (b1-a1)/2;
    m = (b1+a1)/2;
    xt = k * t + m;
    ft = subs(f,x,xt);
    % 计算Gauss点处被积函数值
    for j = 1:1:n
        f_x0(j) = subs(ft, t ,x0(j));
    end
    % 给出数值积分结果，并与精确值比较
    I = I + k * dot(A,f_x0);
end
I0  = vpa(int(f,symvar(sym(f)),a,b));
err = vpa(abs(I - I0));
toc
exp = log10(err);
end

function [x0,Pn_prime,Pn_minus_1] = Gauss(n)
% 预分配内存
Pn_prime   = zeros(n,1);
Pn_minus_1 = zeros(n,1);
% 给出n-1,n阶Legendre多项式
an1  = 1/(2^(n-1) * factorial(n-1));
an   = an1 /(2 * n);
syms x;
Pnm1 = an1* diff((x^2-1)^(n-1),x,n-1);
Pn   = an * diff((x^2-1)^n,    x,n  );
Pnp  = an * diff((x^2-1)^n,    x,n+1);
% 计算n阶Legendre多项式零点及所需函数值
[x0] = solve(Pn);
for i = 1:1:n
    Pn_prime(i)   = subs(Pnp, x,x0(i));
    Pn_minus_1(i) = subs(Pnm1,x,x0(i));
end
end
