function app_Tchebychev(f,a,b,n)
% 主程序
% 分配节点及预分配内存
N = fix(100*(b-a));
x = zeros(N+1);
y = zeros(N+1);
t_ft = zeros(N+1);
T = zeros(N+1);
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
% 计算复合函数f(x(t))在[-1,1]上的n次Tchebychev截断多项式函数
ft = subs(f,x_f,xt);
l = Tchebychev(ft,n);
for i = 1:1:(N+1)
    x(i) = a + (i-1)/N*(b-a);
    t_ft(i) = -1 + (i-1)/N*2;
    y(i) = f(x(i));
    T(i) = subs(l,t,t_ft(i));
    err(i) = T(i) - y(i); 
end
plot(x,y,'b',x,T,'g',x,err,'r');
end

function [l] = Tchebychev(f,k)
% 功能：算法Ⅱ(变量替换)计算f(x)的Tchebychev级数的系数c(0)~c(k)
syms t;
g = subs(f,symvar(sym(f)),cos(sym('t')));
% 预分配内存
P(1:k+1) = t;
c(1:k+1) = 0.0;
% 递归计算系数
P(1) = 1;
P(2) = t;
c(1)  = int(g,t,0,pi)/pi;
c(2)  = int(g*cos(t),t,0,pi)*2/pi;
l = c(1)* P(1) + c(2)* P(2);
for i =3:k+1
    P(i)  = 2*t*P(i-1)-P(i-2);
    c(i)  = int(g*cos((i-1)*t),t,0,pi)*2/pi;
    l = l + c(i)* P(i);
end
end