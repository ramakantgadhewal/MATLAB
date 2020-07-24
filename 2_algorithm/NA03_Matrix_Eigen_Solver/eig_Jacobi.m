function eig_Jacobi(n,ubound,eps)
% n   为输入矩阵的维数.
% ubound 为非对角元的上界.
% err 为对角元即特征值要求的精度.
global A
A = zeros(n);
E = A;
for i=1:1:n-1
    A(i,i) = 2;
    E(i,i+1) = -1;
end
A(n,n) = 2;
A = A + E + E';
delta = ((2*(n-1))^0.5)/n;
err = 2*eps;
% sum  为输出非对角元的平方和.
% lamda为输出特征值列阵.
% tic  为最终扫描次数.
sum = 0;
lamda = zeros(n,1);
tic = 0;
while (delta > ubound) || (err > eps)
   B = A;
   delta = delta / n;
   for i=1:1:n
       for j=1:1:i-1
           if (abs(A(i,j)) > delta)
                Givens(n,i,j);
           end
       end
       for j=i+1:1:n
           if (abs(A(i,j)) > delta)      
                Givens(n,i,j);
           end
       end
   end
   for i=2:1:n
       err = max(abs(B(i-1,i-1)-A(i-1,i-1)),abs(B(i,i)-A(i,i)));
   end
   tic = tic + 1;
   if tic > 1e+2
       error('可能出错！');
   end
end
for i=1:1:n
    for j=1:1:i-1
        sum = sum + A(i,j)*A(i,j);
    end
    for j=i+1:1:n
        sum = sum + A(i,j)*A(i,j);
    end
    lamda(i,1) = A(i,i);
end
sum
lamda
A
tic
end

function Givens(n,k,l)
global A
B = A;
m = 10^(-10);                           % m 为近似计算的预置阈值.   
if B(k,k)==B(l,l)
    A(k,k) = B(k,k)/2 + B(k,l) + B(l,l)/2;
    A(l,l) = B(k,k)/2 - B(k,l) + B(l,l)/2;
    A(k,l) = 0;
    A(l,k) = 0;
    for p=1:1:n
        if (p~=k) && (p~=l)
            A(p,k) = (B(p,k) + B(p,l))*((0.5)^0.5);
            A(p,l) = (-B(p,k) + B(p,l))*((0.5)^0.5);
            A(k,p) = A(p,k);
            A(l,p) = A(p,l);
        end
    end
else
    tt = 2 * B(k,l) / (B(k,k)-B(l,l));   % tt  = tan2O
    if abs(tt) < m                       
        t  = 1 + tt/2;                   % t = tanO
    else
        t  = ((1+tt^2)^0.5-1)/tt;
    end
    if abs(t) < m                       
        t0 = 1 + t/2;                   % t0 = tanO/2
    else
        t0 = ((1+t^2)^0.5-1)/t;
    end
    c  = (1-t0*t0)/(1+t0*t0);           % c  = cosO
    ss = 2*t/(1+t*t);                   % ss = sin2O    
    s  = 2*t0/(1+t0*t0);                % s  = sinO
    A(k,k) = B(k,k)*(c^2) + B(k,l)* ss + B(l,l)*(s^2);
    A(l,l) = B(k,k)*(s^2) - B(k,l)* ss + B(l,l)*(c^2);
    A(k,l) = 0;
    A(l,k) = 0;
    for p=1:1:n
        if (p~=k) && (p~=l)
            A(p,k) =  B(p,k)* c + B(p,l)* s;
            A(p,l) = -B(p,k)* s + B(p,l)* c;
            A(k,p) = A(p,k);
            A(l,p) = A(p,l);
        end
    end
end
end