function VFD02_cal_seperate_pt
dx = 0.0001;
Nx = 2000;
x = zeros(Nx+1);
Z = zeros(Nx+1);
lamb = zeros(Nx+1);
for i = 1:Nx+1
    x(i) = dx*(i-1);
end
Z(1) = 1;
lamb(1) = -1;

% Iteration for Z
for i = 1:Nx
    K1 = dx * DZ(x(i),Z(i));
    K2 = dx * DZ(x(i) + 0.5*dx, Z(i)+0.5*K1);
    K3 = dx * DZ(x(i) + 0.5*dx, Z(i)+0.5*K2);
    K4 = dx * DZ(x(i) +     dx, Z(i)+    K3);
    Z(i+1) = Z(i) + 1/6*(K1 + 2*K2 + 2*K3 + K4);
    lamb(i+1) = -Z(i+1)/(1+x(i))/(1+x(i));
end

% Plot the curve
plot(x,Z,'r-',x,lamb,'b-.');
title('4-order Runge-Kutta method solution');
xlabel('x');
ylabel('Z or \lambda');
legend('Z(x)','\lambda(x)');
end

function [dz] = DZ(x,z)
lamb = -z/(1+x)/(1+x);
dz = g(lamb)*(1+x) + 2/((1+x)^3)*h(lamb)*z*z;
end

function [G] = g(x)
G = (15120 - 2784*x + 79*x*x + 5/3*x*x*x)/(12-x)/(37 + 25/12 * x);
end
function [F] = h(x)
F = (8 + 5/3*x)/(12-x)/(37 + 25/12 * x);
end
