function FST01_self_excited_oscillation
clear all;
close all;

%% close parallel mode
%core_number=2;
%parpool('local',core_number);

%% parameters definition
global h
% ��ֵ����
m = 10;
%rng(100);
% ��ֵ��Χ
xoc = 0; % center location
yoc = 0;
roin_init = 0.0;   % inner radius for generation
rout_init = 1.0;  % outer radius for generation
phil_init = 0*pi; % lower bound of angle
phir_init = 2*pi; % upper bound of angle
xohw_bound = 8;  % half width for visualization
yohw_bound = 6;
% ��������
n = 1e4;
% ��������
h = 1e-5;

%% core iteration
tic
x = zeros(n,2);
for i = 1:m
  k = 1;
  ro_init = roin_init+(rout_init-roin_init)*(i-1)/(m-1);
  phio_init = phil_init + rand*(phir_init-phil_init);
  x(k,1) = xoc + ro_init*cos(phio_init);
  x(k,2) = yoc + ro_init*sin(phio_init);
  plot(x(k,1),x(k,2),'g*');
  hold on
  while(abs(x(k,1)-xoc)<= xohw_bound && ...
        abs(x(k,2)-yoc)<=yohw_bound && k<=n)
    k = k+1;
    [dx,dy] = dxdt2(x(k-1,1),x(k-1,2));
    x(k,1) = x(k-1,1) + dx;
    x(k,2) = x(k-1,2) + dy;
  end
  plot(x(1:k-1,1),x(1:k-1,2),'b-');
  hold on
end
toc

%% plot the trajectories
% for regular fixed points
%plot(xoc,yoc,'ro');
%hold on
%titler = sprintf('%s%s%s',num2str(roin_init,4),'<\delta r_{init}<',num2str(rout_init,4));
%titlephi = sprintf('%s%s%s%s',num2str(phil_init/pi,2),'\pi<\phi_{init}<',num2str(phir_init/pi,2),'\pi');
%title(sprintf('%s%s%s%s%s%s','Trajectories around \eta_0 (',titler,', ',titlephi,')'),'FontName','Cambria Math','FontSize',14);
% for asymptote attactors
title(sprintf('%s%s%s','Trajectories around the asymptote attractor (\Delta t = ',num2str(h),')'));
xlabel('\eta','FontName','Cambria Math','FontSize',14);
ylabel('d\eta/dt','FontName','Cambria Math','FontSize',14);
axis equal;
grid on;

%% close parallel mode
%delete(gcp('nocreate'));
end

%% 4-order Runge-Kutta method
function [dx,dy] = dxdt2(x,y)
global h

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

%% definition of the time-invariant system
function g=f21(x,y)
g = y;
end
function g=f22(x,y)
g = -x*(1-x*x*(0.020069-0.000069444*x*x)) ...
    -0.02*y*(1-y*y*(1.1-0.1*y*y));
end