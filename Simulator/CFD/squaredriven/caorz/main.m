function main(Re,u,n,nt,dif0)

[psi,omg]=base(n);

dif=100;
m=0;
while dif>dif0
    [psi,omg,dif]=operate(Re,psi,omg,n,nt,u);
    m=m+1;
end
m

x=0:1/n:1;
y=x;
[Y,X]=meshgrid(y,x);
figure(1)
contour(X,Y,psi,[psi(1:n,n/2).',psi(n/2,1:n)])
end