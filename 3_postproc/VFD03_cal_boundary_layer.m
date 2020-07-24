function VFD03_cal_boundary_layer(X,f20)
% Solve the ODE
    options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
    [X, Y] = ode45(@bdlayer, [0 X], [0 0 f20], options);
    figure(1);
    plot(X,Y(:,1),'-',X,Y(:,2),'-.',X,Y(:,3),'--');
    legend('f', 'f^\prime', 'f^\prime^\prime')
% Calculate the coefficients
    T_flag = 0;
    if T_flag == 1      % T boundary layer (Prob 10.2)
        Pr = 1;
        dx = X(2,1) - X(1,1);
        [m,n] = size(X);
        g = zeros(m,1);
        A = 0;
        for i = 1:m
            A = A + (Y(i,3)^Pr) * dx;
            g(i) = A;
        end
        g = g./A;
        g = [0;g];
        figure(2);
        plot(X,g(1:m),'r-');
        legend('g(x)');
    end
end

function df = bdlayer(x,f)
% f(1) = f
% f(2) = f'
% f(3) = f''
df = zeros(3,1);
df(1) = f(2);
df(2) = f(3);
df(3) = -f(3)*f(1) + f(2)*f(2) - 1;
end