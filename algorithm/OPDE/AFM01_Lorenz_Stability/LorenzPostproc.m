function LorenzPostproc(B,fs,flag,x0,delx0)
%------------------------------------------------
% Description:
%   Process the data obtained from the Lorenz equations,
%     which has the following character equation:
%   (C + lambda)[lambda^2 + lambda*(A + 1) + A*(1 - B)] = 0.
%     and several critical Points:
%   B = 0; 1; 1.3456; 13.926; 24.06; 24.7368.
% Definition of parameters:
%   B: character parameter
%   fs: sampling frequency
%   flag: choose the certain processor
%       1 - 3D plotting
%       2 - projection to 2D
%       3 - projection to 1D
%       4 - power spectrum density
%       5 - phase spectrum
%       6 - Lyapunov index
%       7 - Poincare section
%   x0, delx0: initial values
% Author:         Yunfan Huang (Tsinghua Univ.)
% Last Updated:   Mar 16
%------------------------------------------------
if nargin == 0
  B = 24;
end
if nargin < 2
  fs = 2^3;
end
if nargin < 3
  flag = [3 4 6];
end
if nargin < 8
  x0 = [0.1 0.1 0.1];
  delx0 = [0.01 0.01 0.01];
end
% read data from .mat file
filename = ...
  sprintf('%s','data_B',num2str(B),...
    '_x',num2str(x0(1,1)),'_y',num2str(x0(1,2)),'_z',num2str(x0(1,3)),...
    '_dx',num2str(delx0(1,1)),'_dy',num2str(delx0(1,2)),'_dz',num2str(delx0(1,3)),...
    '_spec.mat');
load(filename,'x','fn');
ns_origin = size(x,1);
dt_origin = 1/fn;
x0 = x(:,1:3); x0_init = x0(1,:);
x1 = x(:,4:6); x1_init = x1(1,:);
% for projection(1-3) and spectrum(4-5)
magn = floor(fn/fs);
ns = floor(ns_origin/magn);
dt_magn = dt_origin*magn;
fshift = floor(-(ns-1)/2:(ns-1)/2)/ns/dt_magn;
time_series = 1 + magn * (0:ns-1);
% for Lyapunov index(6)
cr_dlam = 1e-16;
% for Poincare section(7)
cr_poin = 100*dt_origin;

%% 3D plotting
if ismember(1,flag)

  figure(1)
  plot3(x0(1,1),x0(1,2),x0(1,3),'b*',x1(1,1),x1(1,2),x1(1,3),'r*');
  legend('x_0','x_1');
  hold on
  plot3(x0(time_series,1),x0(time_series,2),x0(time_series,3),'b-',...
        x1(time_series,1),x1(time_series,2),x1(time_series,3),'r-');
  title('Trajectories of x_0 and x_1');
  legend('x_0','x_1')
  axis('equal');
  grid on;

end

%% projection to 2D
if ismember(2,flag)

  figure(2)
  subplot(1,3,1),plot(x0(time_series,1),x0(time_series,2),'b-',...
                      x1(time_series,1),x1(time_series,2),'r-');
  title('Trajectories projected to x-y plane');
  legend('x_0','x_1');
  axis('equal');
  grid on;
  hold on;
  plot(x0_init(1,1),x0_init(1,2),'b*',...
       x1_init(1,1),x1_init(1,2),'r*');
  legend('x_0','x_1');
  hold off;
  subplot(1,3,2),plot(x0(time_series,1),x0(time_series,3),'b-',...
                      x1(time_series,1),x1(time_series,3),'r-');
  title('Trajectories projected to x-z plane')
  legend('x_0','x_1')
  axis('equal');
  grid on;
  subplot(1,3,3),plot(x0(time_series,2),x0(time_series,3),'b-',...
                      x1(time_series,2),x1(time_series,3),'r-');
  title('Trajectories projected to y-z plane')
  legend('x_0','x_1')
  axis('equal');
  grid on;
  
end

%% projection to 1D
if ismember(3,flag)

  figure(3)
  subplot(3,1,1),plot(time_series,x0(time_series,1),'b-',...
                      time_series,x1(time_series,1),'r-');
  title('Trajectories projected to x-axis','FontName','Cambria Math');
  legend('x_0','x_1','FontName','Cambria Math');
  grid on;
  xlabel('t','FontName','Cambria Math');
  ylabel('x','FontName','Cambria Math');
  subplot(3,1,2),plot(time_series,x0(time_series,2),'b-',...
                      time_series,x1(time_series,2),'r-');
  title('Trajectories projected to y-axis','FontName','Cambria Math')
  legend('x_0','x_1','FontName','Cambria Math');
  xlabel('t','FontName','Cambria Math');
  ylabel('y','FontName','Cambria Math');
  grid on;
  subplot(3,1,3),plot(time_series,x0(time_series,3),'b-',...
                      time_series,x1(time_series,3),'r-');
  title('Trajectories projected to z-axis','FontName','Cambria Math')
  legend('x_0','x_1','FontName','Cambria Math');
  xlabel('t','FontName','Cambria Math');
  ylabel('z','FontName','Cambria Math');
  grid on;
  
end

if ismember(4,flag) || ismember(5,flag)

  Y_x0 = fft(x0(1+magn*(0:ns-1),1))/ns;
  if ismember(4,flag)
%% power spectrum density
    Pshift_x0 = abs(fftshift(Y_x0)).^2*ns*dt_magn;
    figure(4);
    plot(fshift,Pshift_x0);
    title('Double power density spectrum along x-axis',...
          'FontName','Cambria Math');
    xlabel('f(Hz)','FontName','Cambria Math');
    ylabel('{|P(f)|}^2','FontName','Cambria Math');
  
  else
%% phase spectrum
    Ashift_x0 = angle(fftshift(Y_x0))/pi;
    figure(5);
    plot(fshift,Ashift_x0);
    title('Phase spectrum along x-axis',...
          'FontName','Cambria Math');
    xlabel('f(Hz)','FontName','Cambria Math');
    ylabel('Phase(\pi rad)','FontName','Cambria Math');
  
  end
end


%% Lyapunov index
if ismember(6,flag)
  
  k = 2;
  dx0 = norm(x1(1,:)-x0(1,:));
  dx_1 = norm(x1(2,:)-x0(2,:));
  dx_new = norm(x1(k+1,:)-x0(k+1,:));
  lam = log(dx_1/dx0);
  lam_new = (log(dx_new/dx0))/k;
  while abs(lam_new-lam) > cr_dlam && k < ns_origin-1
    k = k+1;
    dx_new = norm(x1(k+1,:)-x0(k+1,:));
    lam = lam_new;
    lam_new = (log(dx_new/dx0))/k;
  end
  lam
  
end

%%  Poincare section
if ismember(7,flag)
  
  a_poin = zeros(ns,2);
  n0 = 0;
  for k=1:ns
    if abs(x0(k,1)-8)<cr_poin
      n0 = n0 + 1;
      a_poin(n0,1) = x0(k,2);
      a_poin(n0,2) = x0(k,3);
    end
  end
  figure(6)
  plot(a_poin(1:n0,1),a_poin(1:n0,2),'og');
  title('Poincare section of x_0 with x = 8',...
        'FontName','Cambria Math');
  xlabel('y','FontName','Cambria Math');
  ylabel('z','FontName','Cambria Math');
  
end
end