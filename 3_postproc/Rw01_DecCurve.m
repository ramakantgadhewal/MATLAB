function Rw01_DecCurve(t,w)
% 拟合公式：w  = (w0 + b)*exp(-a*t) - b
%          -w' = a*(w0 + b)*exp(-a*t) = A*exp(B*t)
%         {A  =  a*(w0 + b); {a = -B       ;
%         {B  = -a         ; {b = -A/B - w0;

n = length(t);
t = t - 1;

% 一种数值微分公式：两点法
w_p1 = zeros(n-1,1);
t_p1 = zeros(n-1,1);
for i=1:1:n-1
    w_p1(i) = log(w(i)-w(i+1));
    t_p1(i) = t(i);
end
p1 = polyfit(t_p1,w_p1,1);
B1 = p1(1);
A1 = exp(p1(2));
wp1 = zeros(n-2,1);
for i=1:1:n-1
    wp1(i) = log(A1) + B1*t(i);
end
% plot(t_p1,w_p1,t_p1,wp1);

a1 = -B1
b1 = -A1/B1 - w(1)
w1 = zeros(n,1);
for i=1:1:n
    w1(i) = (w(1)+b1)*exp(-a1*t(i)) - b1;
end
w_a = w - w1;
plot(t,w,t,w1);


% 另一种数值微分公式：间隔点法
w_p2 = zeros(n-2,1);
t_p2 = zeros(n-2,1);
for i=1:1:n-2
    w_p2(i) = log((w(i)-w(i+2))/2);
    t_p2(i) = t(i);
end
p2 = polyfit(t_p2,w_p2,1);
B2 = p2(1);
A2 = exp(p2(2));
wp2 = zeros(n-2,1);
for i=1:1:n-2
    wp2(i) = log(A2) + B2*t(i);
end
% plot(t_p2,w_p2,t_p2,wp2);

a2 = -B2
b2 = -A2/B2 - w(1)
w2 = zeros(n,1);
for i=1:1:n
    w2(i) = (w(1)+b2)*exp(-a2*t(i)) - b2;
end
w_b = w - w2;
% plot(t,w,t,w2);
% plot(t,w_a,'k',t,w_b);
end

