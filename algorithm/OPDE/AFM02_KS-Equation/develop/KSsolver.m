function KSsolver
%------------------------------------------------
% Description:
%   Use spectual method (i.e. use ut = fft(u) in calc)
%     to solve the K-S equation:
%       u_t + u*u_x + u_{xx} + \sigma*u_{xxxx} = 0.
% Three methods used for time-marching
%   1 for ode45 func
%   2 for 4-order Runge-Kutta method
%   3 for ETDRK4 (Exponential Time-diff & 4-order R-K)
% Author:         Yunfan Huang (Tsinghua Univ.)
% Last Updated:   Apr 16
% Reference: 
%   Four-order Time-stepping for stiff PDEs, SIAM 2005
%    AK Kassam and LN Trefethen (July 2002)
%     doi:10.1137/S1064827502410633
%-----------------------------------------------

clear all;

%% include global parameters
include_flags;

%% prepare the parameters
input_paras;

%% K-S solver
if (tmthd == 1 || tmthd == 2)
  if (corrFlag == 1)
    func = @KS;
  else
    func = @KSorg;
  end
  if (tmthd == 1)
  % 1 for ode45 func
    [tt,ut] = ode45(@(t,ut) func(t,ut),t,ut0);
  else
  % 2 for 4-ordeor Runge-Kutta method
    [tt,ut] = RKdir(@(t,ut) func(t,ut),t,ut0);
  end
elseif (tmthd == 3)
  % 3 for ETDRK4 (Exponential Time-diff & 4-order R-K)
  L = k.^2 - s * k.^4; % Fourier multipliers
  E = exp(dt*L); E2 = exp(dt*L/2);
  M = 32; % number of points for complex means
  r = exp(1i*pi*((1:M)-.5)/M); % roots of unity
  LR = dt*L(:,ones(M,1)) + r(ones(Nx,1),:);
  Q  = dt*real(mean( (exp(LR/2)-1)./LR ,2));
  f1 = dt*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
  f2 = dt*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
  f3 = dt*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
  % Main time loop
  ut = ut0; tt = 0;
  g = -0.5i*k;
  for i = 1:Nt-1
    Nv = g.*fft(real(ifft(ut0)).^2);
    a = E2.*ut0 + Q.*Nv;
    Na = g.*fft(real(ifft(a)).^2);
    b = E2.*ut0 + Q.*Na;
    Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = g.*fft(real(ifft(c)).^2);
    ut0 = E.*ut0 + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    if mod(i,Rfplt)==0
      ut = [ut ut0];
      tt = [tt t(i+1)];
    end
  end
end

%% save data to file
% translate variables to real space
u = real(ifft(ut',[],2));
% write data to .mat file
save(['data/' filename '.mat']);

%% plot the figure for test
figure(1)
subplot(1,2,1);
waterfall(x,tt,u);
view(0,90), colorbar;
xlabel x, ylabel t, zlabel u;
title('Evolution in the x-space');
subplot(1,2,2);
waterfall(fftshift(k),tt,[abs(fftshift(ut',2))]);
view(0,90), colorbar, axis tight;
xlabel k, ylabel t, zlabel fft(u);
title('Evolution in the k-space');

%% postprocessing
postproc(x,tt,ut',pltmode);
  
end


%% Definition of the K-S equation
function dut = KS(t,ut)
include_flags;
  padterm = ifftshift([Zero;fftshift(ut);Zero]);
  crsterm = fftshift(fft(ifft(padterm).^2));
  for i = 1:Nx
    if(abs(ut(i,1)) < 1e-12)
      ut(i,1) = 0;
    end
  end
  dut = - 0.5*(1i*k).*ifftshift(crsterm(Nx/4:Nx*5/4-1,1))...
            + k.^2.*ut - s.*(k.^4).*ut;
end
function dut = KSorg(t,ut)
include_flags;
  for i = 1:Nx
    if(abs(ut(i,1)) < 1e-12)
      ut(i,1) = 0;
    end
  end
  dut = - 0.5*(1i*k).*fft(ifft(ut).^2) + k.^2.*ut - s.*(k.^4).*ut;
end

%% Direct Time-Diff with 4-order R-K method
function [tt,uu] = RKdir(u,t,u0)
include_flags;
  Nt = size(t,1);
  uu = u0;
  tt = 0;
  dt = t(2:Nt) - t(1:Nt-1);
  for i = 1:Nt-1
    du = dudt_rk44(u,t(i),u0,dt(i));
    u0 = u0 + du;
    if (mod(i,Rfplt) == 0)
      uu = [uu u0];
      tt = [tt t(i+1)];
    end
  end
end
function du = dudt_rk44(func,t,u,dt)
  K1 = func(t,u);
  K2 = func(t,u+dt*K1/2);
  K3 = func(t,u+dt*K2/2);
  K4 = func(t,u+dt*K3);
  du = (K1 + 2*K2 + 2*K3 + K4)*dt/6;
end
