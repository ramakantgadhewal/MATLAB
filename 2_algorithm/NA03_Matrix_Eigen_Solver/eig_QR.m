function eig_QR(n,eps)
A = zeros(n);
E = A;
for i=1:1:n-1
    A(i,i) = 2;
    E(i,i+1) = -1;
end
A(n,n) = 2;
A = A + E + E';
l = ones(n,1);
err = 2 * eps;
tic = 0;
while err > eps
    l1 = l;
    [q,r] = qr(A);
    A = r*q;
    l = diag(A);
    err = norm(l-l1,inf);
    tic = tic + 1;
    if tic > 1e+4
        error('可能出错！');
    end
end
l
tic
end

