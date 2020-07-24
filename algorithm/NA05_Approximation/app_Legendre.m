function app_Legendre(f,a,b,n)
% 主程序
% 分配节点及预分配内存
N = fix(100*(b-a));
x = zeros(N+1);
y = zeros(N+1);
t_ft = zeros(N+1);
L = zeros(N+1);
err = zeros(N+1);
% 将任意区间[a,b]变换至[-1,1]
syms x_f xt t;
x_f = symvar(sym(f));
k = (a+b)/2;
if k == 0 
    xt = (b-a)/2*t;
else
    xt = (a+b)/2 * t + (b-a)/2;
end
% 计算复合函数f(x(t))在[-1,1]上的n次Legendre截断多项式函数
ft = subs(f,x_f,xt);
l = Legendre(ft,n);
for i = 1:1:(N+1)
    x(i) = a + (i-1)/N*(b-a);
    t_ft(i) = -1 + (i-1)/N*2;
    y(i) = f(x(i));
    L(i) = subs(l,t,t_ft(i));
    err(i) = L(i) - y(i); 
end
plot(x,y,'b',x,L,'g',x,err,'r');
end

function [l] = Legendre(f,k)
% 功能：计算f(x)的Legendre级数的系数c(0)~c(k)
syms t;
g = subs(f,symvar(sym(f)),sym('t'));
% 预分配内存
P(1:k+1) = t;
c(1:k+1) = 0.0;
% 递归计算系数
P(1) = 1;
P(2) = t;
c(1)=int(g * P(1),t,-1,1)/2;
c(2)=int(g * P(2),t,-1,1)*(2*2-1)/2;
l = c(1)* P(1) + c(2)* P(2);
for i =3:k+1
    P(i) = ((2*i-3)*P(i-1)*t - (i-2)*P(i-2))/(i-1);
    c(i) = int(g * P(i),t,-1,1)*(2*i-1)/2;
    l = l + c(i)* P(i);
end
end