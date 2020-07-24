function FM01_check_stability(b,d1,d2,n)
% 分析无粘近似下边界层流动的稳定性：
%        { 2by,               0 <= y <= 1/2 (Ⅰ) 
% U(y) = { 2(1-b)y+2b-1,     1/2 < y <= 1   (Ⅱ)
%        { 1,                      y > 1    (Ⅲ)
% 参数说明：
%   d1、d2       区间端点
%   n            作图节点数
%   b            待定参数（b<0.5对应有拐点速度剖面，b>0.5对应无拐点速度剖面）
l   = length(b);
a   = zeros(n,1);
a0  =    0;
b1  =  1.5;
b2  = -1.5;
sym = zeros(n,l);
for j = 1:l
    for i = 1:n
        a(i) = d1 + i/n*(d2-d1);
        K1   = (2*b(j)-1)/a(i);
        K2   = exp(a(i))/(1-exp(a(i)));
        K3    = 1 - b(j);
        a1   = a(i)*K2/K3;
        a2   = (-exp(-a(i))+a(i)-1)*K2+exp(-a(i));
        a3   = ( exp(-a(i))+a(i)-1)*K1;
        delta = (a2)^2 - 4*a1*a3;
        if delta >= 0
            sym(i,j)  =  1 - j/(2*l);
        else
            sym(i,j)  = -1 + j/(2*l);
        end
    end
end
plot(a0,b1,a0,b2);
hold on
for k = 1:l
    if sym(:,k) > 0
        p(k) = plot(a,sym(:,k),'k');
    else
        p(k) = plot(a,sym(:,k),'g');
    end
end
legend([p],...
    ['b(1) = ',num2str(b(1))], ...
    ['b(2) = ',num2str(b(2))]);
% legend([p(1),p(2)],['b(1) = ',num2str(b(1))],['b(2) = ',num2str(b(2))]);
title({'二维平板边界层线性稳定性';'(Note.sym>0为稳定区间)'});
xlabel('\alpha');
ylabel('sym');
end