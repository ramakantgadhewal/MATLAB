function [psi,omg,dif]=operate(Re,psi0,omg0,n,nt,u0)

dx=1/n;
dy=1/n;
dt=1/nt;

psi1=zeros(n+5);
psi1(3:n+3,3:n+3)=psi0(:,:);
psi1(4:n+2,1)=psi0(2:n,2);
psi1(4:n+2,2)=psi0(2:n,2);
psi1(4:n+2,n+4)=psi0(2:n,n)+2*dy*u0;
psi1(4:n+2,n+5)=psi0(2:n,n)+3*dy*u0;
psi1(n+4,4:n+2)=psi0(n,2:n);
psi1(n+5,4:n+2)=psi0(n,2:n);
psi1(1,4:n+2)=psi0(2,2:n);
psi1(2,4:n+2)=psi0(2,2:n);

u=zeros(n+1);
v=zeros(n+1);
for i=1:n+1
    for j=1:n+1
        u(i,j)=(psi1(i+2,j+3)-psi1(i+2,j+1))/(2*dy);
        v(i,j)=-(psi1(i+3,j+2)-psi1(i+1,j+2))/(2*dx);
    end
end

omg1=zeros(n+3);
omg1(2:n+2,2:n+2)=omg0(:,:);
for i=3:n+1
    omg1(1,i)=(psi1(3,i)-2*psi1(2,i)+psi1(1,i))/(dx^2);
    omg1(n+3,i)=(psi1(n+5,i)-2*psi1(n+4,i)+psi1(n+3,i))/(dx^2);
    omg1(i,1)=(psi1(i,3)-2*psi1(i,2)+psi1(i,1))/(dy^2);
    omg1(i,n+3)=(psi1(i,n+5)-2*psi1(i,n+4)+psi1(i,n+3))/(dy^2);
end
omg=zeros(n+1);
for i=1:n+1
    for j=1:n+1
        omg(i,j)=omg0(i,j)-u(i,j)*(dt/(2*dx))*(omg1(i+2,j+1)-omg1(i,j+1))-v(i,j)*(dt/(2*dy))*(omg1(i+1,j+2)-omg1(i+1,j))+(1/Re)*((dt/(dx^2))*(omg1(i+2,j+1)+omg1(i,j+1)-2*omg1(i+1,j+1))+(dt/(dy^2))*(omg1(i+1,j+2)+omg1(i+1,j)-2*omg1(i+1,j+1)));
    end
end

dif=100;
psik=psi0;


m=0;
while dif>1 && m<50
    psit=zeros(n+3);
    psit(2:n+2,2:n+2)=psik(:,:);
    psit(3:n+1,1)=psik(2:n,2);
    psit(3:n+1,n+3)=psik(2:n,n)+2*dy*u0;
    psit(1,3:n+1)=psik(2,2:n);
    psit(n+3,3:n+1)=psik(n,2:n);
    dif=0;
    for i=1:n+1
        for j=1:n+1
            psik(i,j)=0.25*(psit(i+2,j+1)+psit(i,j+1)+psit(i+1,j+2)+psit(i+1,j)-dx*dy*omg(i,j));
            dif=dif+(psik(i,j)-psit(i+1,j+1))^2;
        end
    end
    m=m+1;
end
psi=psik;

dif=0;
for i=1:n+1
    for j=1:n+1
        dif=dif+(psi(i,j)-psi0(i,j))^2;
    end
end
end