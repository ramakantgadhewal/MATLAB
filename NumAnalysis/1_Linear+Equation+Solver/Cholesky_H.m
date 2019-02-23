function Cholesky_H(n,a,ind)
% ind = 0: without Tikhonov (n<13)
%     = 1: with Tikhonov
H = ones(n);
I = zeros(n);
x0 = ones(n,1);
for i=1:1:n
    I(i,i) = 1;
    for j=1:1:n   
        H(i,j)=1/(i+j-1);
    end
end
b = H * x0;
T = H'* H;
H1 = a * I + T;
b1 = T * x0;
if ind == 0
    R = chol(H);
    x = R\(R'\b);
    e = x - x0;
    err = norm(e,2)
elseif ind == 1
    R1 = chol(H1);
    x1 = R1\(R1'\b1)
    e1 = x1 - x0;
    err1 = norm(e1,2)
end

