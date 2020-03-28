function postproc(x,t,ut,pltmode)
%------------------------------------------------
% Description:
%   Process the data obtained from the K-S equation
%     u_t + u*u_x + u_{xx} + \sigma*u_{xxxx} = 0.
% Parameters:
%   x: grids in real space
%   t: time points for plotting
%   ut: solution points for plotting
%   pltmode: plotting mode
%       1 for plain plotting of x(t) and k(t)
%       2 for phase diagram of (E(t),dE(t))
%       3 for bifurcation diagram of Emin
% Author:         Yunfan Huang (Tsinghua Univ.)
% Last Updated:   Apr 16
%------------------------------------------------

include_flags;

% translate variables to real space
u = real(ifft(ut,[],2));

if ismember(1,pltmode)
  figure(1)
  subplot(1,2,1);
  waterfall(x,t,u);
  view(0,90), colorbar;
  xlabel x, ylabel t, zlabel u;
  title('Evolution in the x-space');
  subplot(1,2,2);
  waterfall(fftshift(k),t,[abs(fftshift(ut,2))]);
  view(0,90), colorbar, axis tight;
  xlabel k, ylabel t, zlabel fft(u);
  title('Evolution in the k-space');
  savefig(['data/' filename '_1.fig']);
end
if ismember(2,pltmode)
  E = sum(u.^2,2);
  Ephs = E(1:end-1);
  dEphs = (E(2:end)-E(1:end-1))/dt;
  figure(2)
  subplot(1,2,1)
  plot(t,E); xlabel t, ylabel E;
  title('Evolution of total energy');
  subplot(1,2,2)
  plot(Ephs,dEphs,Ephs(1),dEphs(1),'g*');
  xlabel E, ylabel dE/dt;
  title('Phase diagram of E-dE/dt');
  savefig(['data/' filename '_2.fig']);
end
if ismember(3,pltmode)
  
end

end

