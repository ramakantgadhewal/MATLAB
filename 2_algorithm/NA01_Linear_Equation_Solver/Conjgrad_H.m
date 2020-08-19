function Conjgrad_H(n,eps)
H = ones(n);
x_0 = ones(n,1);
x0  = zeros(n,1);
for i=1:1:n
    for j=1:1:n   
        H(i,j)=1/(i+j-1);
    end
end
b = H * x_0;
% without preprocessing
if nargin == 1
    eps = 1.0e-8;
end
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
n_1
err1 = norm(x-x_0)
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
n_2
err2 = norm(x-x_0)
end