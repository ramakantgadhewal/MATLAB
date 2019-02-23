function interp_spline(a,b,n)
% 等距三次样条插值主程序
% 功能：A.插值节点赋值[subroutine f(x)]
%       B.节点补充条件[subroutine coffi(y,n)]
%       C.描点作图[subroutine spline(t,x,y,n)]

% A
x = zeros(n+1,1);
y = zeros(n+1,1);
for j = 1:1:(n+1)
    x(j) = a + (b-a)*(j-1)/n;
    y(j) = f(x(j));
end
% B
global m;
m = coffi(y,n);
% C
N = 10 * n;
x_0 = zeros(N+1);
F = zeros(N+1);
S = zeros(N+1);
err = zeros(N+1);
for i = 1:1:(N+1)
    x_0(i) = a + (b-a)*(i-1)/N;
    F(i) = f(x_0(i)); 
    S(i) = spline(x_0(i),x,y,n);
    err(i) = F(i) - S(i);
end
plot(x_0,F,'r',x_0,S,'b',x_0,err,'g');
end


function [y] = f(x)
% 被插值函数
y = 1/(1+25*x^2);
end


function [m] = coffi(y,n)
% 求取补充条件：样条函数节点处一阶导数值（自然边界条件）
% subsubroutine:followup(A,g)

a = zeros(n,1); % save lambda_1 to n
b = zeros(n,1); % save mu_0 to n-1
c = zeros(n+1,1); % save 非对角系数
g = zeros(n+1,1); % save 右端向量
a(:) = 1/2;
b(:) = 1/2;
c(:) = 2;
for i = 2:1:n
    g(i) = 3 * a(1) * ((y(i)-y(i-1)) * n/2 + (y(i+1)-y(i)) * n/2);
end
g(1) = 3 * (y(2)-y(1)) * n/2;
g(n+1) = 3 * (y(n+1)-y(n)) * n/2;
a(n) = 1;
b(1) = 1; 

A1 = diag(c);
A2 = diag(b,1);
A3 = diag(a,-1);
A = A1 + A2 + A3;
m = followup(A,g);
end
function [x] = followup(A,g)
% 追赶法求解严格三对角

n = length(g) - 1;
x = zeros(n+1,1);
c = diag(A);
a = diag(A,-1);
b = diag(A,1);
for i = 2:1:(n+1)
    c(i,1) = c(i,1) - (a(i-1,1)/c(i-1,1))*b(i-1,1);
    g(i,1) = g(i,1) - (a(i-1,1)/c(i-1,1))*g(i-1,1);
end
x(n+1) = g(n+1,1)/c(n+1,1);
for i = n:-1:1
    x(i) = (g(i,1)-b(i,1)*x(i+1,1))/c(i,1);
end
end


function [l] = spline(t,x,y,n)
% 三次样条插值函数
% subsubroutine:alpha(t,x,k,n),beta(t,x,k,n)

global m;
for k = 1:1:(n+1)
    y(k) = y(k) * alpha(t,x,k,n)+ m(k) * beta(t,x,k,n);
end
l = sum(y);
end
function [p] = alpha(t,x,k,n)
% 分段Hermite插值基Ⅰ
p = 0;
if (k == 1) && (t <= x(2))
    p = (1 + 2*(t-x(1))/(x(2)-x(1)))*((t-x(2))/(x(1)-x(2)))^2;
end
if (k == n + 1) && (t >= x(n))
    p = (1 + 2*(t-x(n+1))/(x(n)-x(n+1)))*((t-x(n))/(x(n+1)-x(n)))^2;
end
if (k >= 2) && (k <= n)
    if (t >= x(k-1)) && (t <= x(k))
        p = (1 + 2*(t-x(k))/(x(k-1)-x(k)))*((t-x(k-1))/(x(k)-x(k-1)))^2;
    else
        if (t >= x(k)) && (t <= x(k+1))
            p = (1 + 2*(t-x(k))/(x(k+1)-x(k)))*((t-x(k+1))/(x(k)-x(k+1)))^2;
        end
    end
end
end
function [p] = beta(t,x,k,n)
% 分段Hermite插值基Ⅱ
p = 0;
if (k == 1) && (t <= x(2))
    p = (t-x(1))*((t-x(2))/(x(1)-x(2)))^2;
end
if (k == n + 1) && (t >= x(n))
    p = (t-x(n+1))*((t-x(n))/(x(n+1)-x(n)))^2;
end
if (k >= 2) && (k <= n)
    if (t >= x(k-1)) && (t <= x(k))
        p = (t-x(k))*((t-x(k-1))/(x(k)-x(k-1)))^2;
    else
        if (t >= x(k)) && (t <= x(k+1))
            p = (t-x(k))*((t-x(k+1))/(x(k)-x(k+1)))^2;
        end
    end
end
end