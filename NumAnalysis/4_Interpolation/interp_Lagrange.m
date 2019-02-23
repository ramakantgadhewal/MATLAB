function interp_Lagrange(a,b,n)
% Lagrange插值主程序
% 插值节点赋值
x = zeros(n+1,1);
y = zeros(n+1,1);
for j = 1:1:(n+1)
    x(j) = a + (b-a)*(j-1)/n;
    y(j) = f(x(j));
end
% 描点作图
N = 10 * n;
x_0 = zeros(N+1);
F = zeros(N+1);
L = zeros(N+1);
err = zeros(N+1);
for i = 1:1:(N+1)
    x_0(i) = a + (b-a)*(i-1)/N;
    F(i) = f(x_0(i)); 
    L(i) = Lagrange(x_0(i),x,y,n);
    err(i) = F(i) - L(i);
end
plot(x_0,F,'r',x_0,L,'b',x_0,err,'g');
end


function [y] = f(x)
% 被插值函数
y = 1/(1+25*x^2);
end


function [l] = Lagrange(t,x,y,n)
% Lagrange插值函数（内部已耦合Lagrange插值基）
for k = 1:1:(n+1)
    for m = 1:1:(n+1)
        if m ~= k
            y(k) = y(k)*(t-x(m))/(x(k)-x(m));
        end
    end
end
l = sum(y);
end
