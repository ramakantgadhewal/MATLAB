function res=upwind(vort,u,v,Dt,Dx,Dy,Re)
upl=(u+abs(u))/2;
umi=(u-abs(u))/2;
vpl=(v+abs(v))/2;
vmi=(v-abs(v))/2;
na=size(vort);
n=na(1);
dxvort=vort(:,2:n)-vort(:,1:n-1);
dyvort=vort(1:n-1,:)-vort(2:n,:);
d2xvort=zeros(n,n);
d2xvort(:,2:n-1)=dxvort(:,2:n-1)-dxvort(:,1:n-2);
d2xvort(:,1)=dxvort(:,2);
d2xvort(:,n)=dxvort(:,n-1);
d2yvort=zeros(n,n);
d2yvort(2:n-1,:)=dyvort(1:n-2,:)-dyvort(2:n-1,:);
d2yvort(1,:)=d2yvort(2,:);
d2yvort(n,:)=d2yvort(n-1,:);
alpha0=1;
res=vort(2:n-1,2:n-1)-Dt/Dx*upl(2:n-1,2:n-1).* ...
    (dxvort(2:n-1,1:n-2)+alpha0/2*d2xvort(1:n-2)) ...
    -Dt/Dx*umi(2:n-1,2:n-1).* ...
    (dxvort(2:n-1,2:n-1)-alpha0/2*d2xvort(2:n-1,3:n))- ...
    Dt/Dy*vpl(2:n-1,2:n-1).* ...
    (dyvort(2:n-1,2:n-1)+ alpha0/2*d2yvort(3:n,2:n-1))- ...
    Dt/Dy*vmi(2:n-1,2:n-1).* ...
    (dyvort(1:n-2,2:n-1)-alpha0/2*d2yvort(1:n-2,2:n-1))+ ...
    Dt/Re*(d2xvort(2:n-1,2:n-1)/Dx^2+d2yvort(2:n-1,2:n-1)/Dy^2);
end