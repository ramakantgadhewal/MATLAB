function Conjgrad(m,eps)
if nargin == 1
    eps = 1.0e-8;
end
T = zeros(m,1);
M = zeros(m,1);
for n =6:1:m
    H = ones(n);
    L = zeros(n);
    x_0 = ones(n,1);
    x0  = zeros(n,1);
    for i=1:1:n
        for j=1:1:n   
            H(i,j)=1/(i+j-1);
        end
    end
    b = H * x_0;
% without preprocessing
    r1 = b-H*x0;
    p1 = r1;
    d  = dot(r1,r1)/dot(p1,H*p1);
    x  = x0+d*p1;
    r2 = r1-d*H*p1;
    f  = dot(r2,r2)/dot(r1,r1);
    p2 = r2+f*p1;
    n_1= 1;
    for i=1:(n-1)
        x0 = x;
        p1 = p2;
        r1 = r2;
        d  = dot(r1,r1)/dot(p1,H*p1);
        x  = x0+d*p1;
        r2 = r1-d*H*p1;
        f  = dot(r2,r2)/dot(r1,r1);
        p2 = r2+f*p1;
        n_1= n_1 + 1;
    end
    d  = dot(r1,r1)/dot(p1,H*p1);
    x  = x+d*p2;
    err1 = norm(x-x_0);
% with preprocessing
    F = zeros(n);
    for i=1:n
        F(i,i) = sqrt(2*i-1);
    end
    H = F*(H*F');
    b = F*b;
    r1 = b-H*x0;
    p1 = r1;
    d  = dot(r1,r1)/dot(p1,H*p1);
    x  = x0+d*p1;
    r2 = r1-d*H*p1;
    f  = dot(r2,r2)/dot(r1,r1);
    p2 = r2+f*p1;
    n_2= 1;
    for i=1:(n-1)
        x0 = x;
        p1 = p2;
        r1 = r2;
        d  = dot(r1,r1)/dot(p1,H*p1);
        x  = x0+d*p1;
        r2 = r1-d*H*p1;
        f  = dot(r2,r2)/dot(r1,r1);
        p2 = r2+f*p1;
        n_2= n_2 + 1;
    end
    d  = dot(r1,r1)/dot(p1,H*p1);
    x  = x+d*p2;
    x  = F'*x;
    err2 = norm(x-x_0);
    T(n,1)=log(max(err1/err2,err2/err1));
    M(n,1)=n/m;
end
plot(M,T);
xlabel('方程组阶数m /实验最高阶数n');
ylabel('ln(较大误差/较小误差)');
end

