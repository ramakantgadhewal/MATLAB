function Gauss_H(n,a)
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
[~,B,~]=svd(H);
B(1,1)
B(n,n)
H1 = a * I + T;
b1 = T * x0;
[~,B1,~]=svd(H);
B1(1,1)
B1(n,n)
cd = cond(H)
cd1 = cond(H1)
    %without tikhonov
    [L,U] = lu(H);
    x = U\(L\b);
    e = x - x0;
    err = norm(e,2)
    %with tikhonov
    [L1,U1] = lu(H1);
    x1 = U1\(L1\b1);
    e1 = x1 - x0;
    err1 = norm(e1,2)
end

