function LorenzSolver(x0,delx0,ns,fn)
%------------------------------------------------
% Description:
%   Use 4-order Runge-Kutta method to solve the Lorenz equations,
%     which has the following character equation:
%   (C + lambda)[lambda^2 + lambda*(A + 1) + A*(1 - B)] = 0.
%     and several critical Points:
%   B = 0; 1; 1.3456; 13.926; 24.06; 24.7368.
% Definition of parameters:
%   ns: sample number
%   fn: numerical frequency
% Author:         Yunfan Huang (Tsinghua Univ.)
% Last Updated:   Mar 16
%------------------------------------------------

% set the default value
if nargin < 2
  x0 = [0.1 0.1 0.1];
  delx0 = [1e-2 1e-2 1e-2];
end
if nargin < 4
  ns = 2^10;
  fn = 2^8;
end
dt = 1/fn;    % time step
% set global variables
global A B C
A = 10;
B = 22.5;
C = 8/3;
% initialization of the trajectory vector
x = zeros(ns,6);
x1 = x0 + delx0;
x(1,1:3) = x0;
x(1,4:6) = x1;

%% calculate the trajectory
tic
for k=2:ns
  [dx0,dy0,dz0] = dxdt_Lorenz3(x(k-1,1),x(k-1,2),x(k-1,3),dt);
  [dx1,dy1,dz1] = dxdt_Lorenz3(x(k-1,4),x(k-1,5),x(k-1,6),dt);
  x(k,1) = x(k-1,1) + dx0;
  x(k,2) = x(k-1,2) + dy0;
  x(k,3) = x(k-1,3) + dz0;
  x(k,4) = x(k-1,4) + dx1;
  x(k,5) = x(k-1,5) + dy1;
  x(k,6) = x(k-1,6) + dz1;
end
toc

%% save the results to file
filename = ...
  sprintf('%s','data_B',num2str(B),...
    '_x',num2str(x0(1,1)),'_y',num2str(x0(1,2)),'_z',num2str(x0(1,3)),...
    '_dx',num2str(delx0(1,1)),'_dy',num2str(delx0(1,2)),'_dz',num2str(delx0(1,3)),...
    '_spec.mat');
save(filename,'x','x0','delx0','A','B','C','fn');

% figure for test
figure
plot3(x0(1,1),x0(1,2),x0(1,3),'b*',x1(1,1),x1(1,2),x1(1,3),'r*');
legend('x_0','x_1');
hold on
plot3(x(:,1),x(:,2),x(:,3),'b-',x(:,4),x(:,5),x(:,6),'r-');
legend('x_0','x_1')
axis('equal');
grid on;

end

%% 4-order Runge-Kutta method
function [dx,dy,dz] = dxdt_Lorenz3(x,y,z,h)

K1 = f31(x,y,z);
L1 = f32(x,y,z);
M1 = f33(x,y,z);

K2 = f31(x + h*K1/2,y + h*L1/2,z + h*M1/2);
L2 = f32(x + h*K1/2,y + h*L1/2,z + h*M1/2);
M2 = f33(x + h*K1/2,y + h*L1/2,z + h*M1/2);

K3 = f31(x + h*K2/2,y + h*L2/2,z + h*M2/2);
L3 = f32(x + h*K2/2,y + h*L2/2,z + h*M2/2);
M3 = f33(x + h*K2/2,y + h*L2/2,z + h*M2/2);

K4 = f31(x + h*K3,y + h*L3, z + h*M3);
L4 = f32(x + h*K3,y + h*L3, z + h*M3);
M4 = f33(x + h*K3,y + h*L3, z + h*M3);

dx = (K1 + 2*K2 + 2*K3 + K4)*h/6;
dy = (L1 + 2*L2 + 2*L3 + L4)*h/6;
dz = (M1 + 2*M2 + 2*M3 + M4)*h/6;
end

%% definition of the time-invariant system
function g=f31(x,y,z)
global A B C
g = - A* x + A * y;
end
function g=f32(x,y,z)
global A B C
g = B * x - y - x*z;
end
function g=f33(x,y,z)
global A B C
g = -C * z + x*y;
end