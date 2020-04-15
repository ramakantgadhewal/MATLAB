%vort-flow
vort=zeros(129,129);
flow=zeros(129,129);
U=zeros(129,129);
V=zeros(129,129);
Dx=1/128;
Dy=1/128;
Re=400;
Dt=1/1024;


%initialize
U(1,2:128)=1;
vort(1,2:128)=2/Dx;
vort(1,1)=1/Dx;
vort(1,129)=1/Dx;

ind=1;
iterator=0;
while (ind>1e-4)
    iterator=iterator+1;
    Utemp=U;
    Vtemp=V;
    vorttemp=vort;
    vort(2:128,2:128)=vort(2:128,2:128)-Dt/2/Dx*U(2:128,2:128).*(vort(2:128,3:129)-vort(2:128,1:127))+Dt/2/Dy*V(2:128,2:128).*(vort(1:127,2:128)-vort(3:129,2:128))+Dt/Re*((vort(2:128,1:127)-2*vort(2:128,2:128)+vort(2:128,3:129))/Dx^2+(vort(1:127,2:128)-2*vort(2:128,2:128)+vort(3:129,2:128))/Dy^2);
%     vort(2:128,2:128)=upwind(vort,U,V,Dt,Dx,Dy,Re);
    phiind = 1;
    while (phiind>1e-6)
        flowtemp=flow;
        flow(2:128,2:128)=0.25*(flow(2:128,1:127)+flow(2:128,3:129)+flow(1:127,2:128)+flow(3:129,2:128)-Dx^2*vort(2:128,2:128));
        flowtemp=abs(flow-flowtemp);
        phiind=max(max(flowtemp));
    end
    
    U(2:128,2:128)=(flow(1:127,2:128)-flow(3:129,2:128))/(2*Dy);
    V(2:128,2:128) = -(flow(2:128,3:129)-flow(2:128,1:127))/(2*Dx);
    
    vort(2:128,1)=2*(flow(2:128,2)-flow(2:128,1))/Dx^2;
    vort(2:128,129)=2*(flow(2:128,128)-flow(2:128,129))/Dx^2;
    vort(129,2:128)=2*(flow(128,2:128)-flow(129,2:128))/Dy^2;
    vort(1,2:128)=2*(flow(2,2:128)-flow(1,2:128)+1*Dy)/Dy^2;
    vort(1,1)=(vort(1,2)+vort(2,1))/2;
    vort(1,129)=(vort(1,128)+vort(2,129))/2;
    vort(129,129)=(vort(128,129)+vort(129,128))/2;
    vort(129,1)=(vort(128,1)+vort(129,2))/2;
    
    ind=max(max(max(U-Utemp)),max(max(V-Vtemp)));
    if(mod(iterator,100)==0)
    fprintf('%s%i%s%f\n','the ',iterator,'th step get accuracy =',ind);
    end
    
    x = 0:1/128:1;
    y = 0:1/128:1;
    [X,Y] = meshgrid(y,x);
    contour(X,Y,flow,30);

end

