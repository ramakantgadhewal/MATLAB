function plot_shockwave

    include_flags;

    x = zeros(Nx+1,1);
    t = zeros(Nt+1,1);
    for i = 1:Nx+1
        x(i,1) = xleft + dx*(i-1);
    end
    for j = 1:Nt+1
        t(j,1) = dt*(j-1);
    end

    figure(1);
    mesh(t,x,rho);
    xlabel('t');
    ylabel('x');
%     
%     figure(2);
%     mesh(t,x,u);
%     
%     figure(3);
%     mesh(t,x,p);
    
    figure(4);
    plot(x(:,1),rho(:,Time/dt),'r-',x(:,1),u(:,Time/dt),'b-',x(:,1),p(:,Time/dt),'k-');
    xlabel('x');
    ylabel('Conservative variables');
    legend('\rho','u','p');

end