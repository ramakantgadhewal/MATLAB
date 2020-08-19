function interp_linear(a,b,n)
% Piecewise linear interpolation

% assign the interpolation points
x = zeros(n+1,1);
y = zeros(n+1,1);
for j = 1:1:(n+1)
    x(j) = a + (b-a)*(j-1)/n;
    y(j) = f(x(j));
end
% plot the figure
N = 10 * n;
x_0 = zeros(N+1);
F = zeros(N+1);
L = zeros(N+1);
err = zeros(N+1);
for i = 1:1:(N+1)
    x_0(i) = a + (b-a)*(i-1)/N;
    F(i) = f(x_0(i)); 
    L(i) = linear(x_0(i),x,y,n);
    err(i) = F(i) - L(i);
end
plot(x_0,F,'r',x_0,L,'b',x_0,err,'g');
end


function [y] = f(x)
% function to be interpolated
y = 1/(1+25*x^2);
end


function [l] = linear(t,x,y,n)
% linear interpolation function
for k = 1:1:(n+1)
    y(k) = y(k) * phi(t,x,k,n);
end
l = sum(y);
end

function [p] = phi(t,x,k,n)
% linear basis
p = 0;
if (k == 1) && (t <= x(2))
    p = (t-x(2))/(x(1)-x(2));
end
if (k == n + 1) && (t >= x(n))
    p = (t-x(n))/(x(n+1)-x(n));
end
if (k >= 2) && (k <= n)
    if (t >= x(k-1)) && (t <= x(k))
        p = (t-x(k-1))/(x(k)-x(k-1));
    else
        if (t >= x(k)) && (t <= x(k+1))
            p = (t-x(k+1))/(x(k)-x(k+1));
        end
    end
end
end


