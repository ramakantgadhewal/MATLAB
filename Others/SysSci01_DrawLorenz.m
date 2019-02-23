function SysSci01_DrawLorenz
% ============================
% -- 描述: 绘画洛伦兹方程的相图
% ============================

    global A B C Re
    % (C + lambda)[lambda^2 + lambda*(A + 1) + A*(1 - B)] = 0.
    % 节点 B = 0; 1; 1.3456; 13.926; 24.06; 24.7368.
    A = 10;
    B = 26;
    C = 8/3;
    Re = 1;
    x0 = 0.1;
    y0 = 0.1;
    z0 = 0.1;
%     
%     global mu
%     mu = 1;
%     x0 = 4;
%     y0 = 6;

%     global lambda1 lambda2
%     lambda1 = -2;
%     lambda2 = -3;
%     x0 = -1;

    global n h
    % 相空间维数
    nd = 3;
    % 迭代步数
    n = 1e4;
    % 迭代步长
    h = 1e-2;
    tic
    if nd == 1
        x = zeros(n,1);
        t = zeros(n,1);
        x(1,1) = x0;
        for k=2:n
            dx = dxdt_Lorenz1(x(k-1,1));
            x(k,1) = x(k-1,1) + dx;
            t(k,1) = t(k-1,1) + h;
        end
        plot(t(:,1),x(:,1));
        grid on;
    elseif nd == 2
        x = zeros(n,2);
        x(1,:) = [x0,y0];
        for k=2:n
            [dx,dy] = dxdt_Lorenz2(x(k-1,1),x(k-1,2));
            x(k,1) = x(k-1,1) + dx;
            x(k,2) = x(k-1,2) + dy;
        end
        plot(x(:,1),x(:,2));
        grid on;
    elseif nd == 3
        x = zeros(n,3);
        x(1,:) = [x0,y0,z0];
        for k=2:n
            [dx,dy,dz] = dxdt_Lorenz3(x(k-1,1),x(k-1,2),x(k-1,3));
            x(k,1) = x(k-1,1) + dx;
            x(k,2) = x(k-1,2) + dy;
            x(k,3) = x(k-1,3) + dz;
        end
        plot3(x(:,1),x(:,2),x(:,3));
        grid on;
    end
    toc
end

%% 4-order Runge-Kutta method
function [dx] = dxdt_Lorenz1(x)
global n h

K1 = f11(x);
K2 = f11(x + h*K1/2);
K3 = f11(x + h*K2/2);
K4 = f11(x + h*K3);

dx = (K1 + 2*K2 + 2*K3 + K4)*h/6;
end
function [dx,dy] = dxdt_Lorenz2(x,y)
global n h

K1 = f21(x,y);
K2 = f21(x + h*K1/2,y + h*K1/2);
K3 = f21(x + h*K2/2,y + h*K2/2);
K4 = f21(x + h*K3,y + h*K3);

L1 = f22(x,y);
L2 = f22(x + h*L1/2,y + h*L1/2);
L3 = f22(x + h*L2/2,y + h*L2/2);
L4 = f22(x + h*L3,y + h*L3);

dx = (K1 + 2*K2 + 2*K3 + K4)*h/6;
dy = (L1 + 2*L2 + 2*L3 + L4)*h/6;
end
function [dx,dy,dz] = dxdt_Lorenz3(x,y,z)
global n h

K1 = f31(x,y,z);
K2 = f31(x + h*K1/2,y + h*K1/2,z + h*K1/2);
K3 = f31(x + h*K2/2,y + h*K2/2,z + h*K2/2);
K4 = f31(x + h*K3,y + h*K3, z + h*K3);

L1 = f32(x,y,z);
L2 = f32(x + h*L1/2,y + h*L1/2,z + h*L1/2);
L3 = f32(x + h*L2/2,y + h*L2/2,z + h*L2/2);
L4 = f32(x + h*L3,y + h*L3, z + h*L3);

M1 = f33(x,y,z);
M2 = f33(x + h*M1/2,y + h*M1/2,z + h*M1/2);
M3 = f33(x + h*M2/2,y + h*M2/2,z + h*M2/2);
M4 = f33(x + h*M3,y + h*M3, z + h*M3);

dx = (K1 + 2*K2 + 2*K3 + K4)*h/6;
dy = (L1 + 2*L2 + 2*L3 + L4)*h/6;
dz = (M1 + 2*M2 + 2*M3 + M4)*h/6;
end

%% definition of the time-invariant system
function g=f11(x)
global lambda1 lambda2
g = lambda1 + lambda2*x - x*x*x;
end

function g=f21(x,y)
global mu
g = x*(1 + 3*x - x*x - y);
end
function g=f22(x,y)
global mu
g = y*(x-1);
end

function g=f31(x,y,z)
global A B C Re
g = - A* x + A * y;
end
function g=f32(x,y,z)
global A B C Re
g = B * x - y - x*z;
end
function g=f33(x,y,z)
global A B C Re
g = -C * z + x*y;
end