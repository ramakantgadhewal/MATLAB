include_flags;

%% Prescribed parameters
% sigma in K-S eqn
%   u_t + u*u_x + u_{xx} + \sigma*u_{xxxx} = 0.
s = 0.02;
% total iteration time (defining the range of temporal domain)
T = 1e2;
% time step for iteration and plotting
dt = 1e-2; dtplt = 1e-2;
% plotting mode
%     1 for plain plotting of x(t) and k(t)
%     2 for phase diagram of (E(t),dE(t))
%     3 for bifurcation diagram of Emin
pltmode = [2];
% Lexp = log2(L/(2*pi)) (defining the resolution of spectual space)
Lexp = 2;
% kexp = log2(kmax) (defining the detectable range of spectual space)
kexp = 5;
kexpsignal = 0;
% Nxexp = log2(Nx) (related to the computational efficiency)
%
% Note.
% - For a interval of length L, the max wave number that could
%     be detected from Nx sample points is:
%         kmax = (2*pi/L)*(Nx/2).
% - Proof:
%     According to the sampling theorem, the sampling frequency
%     must be at least twice the signal frequency in order to 
%     detect the certain signal.
%     The signal frequency is: fsignal = 1/lambda = kmax/2/pi,
%     The sampling frequency is: fsample = 1/dx = L/Nx.
%     Therefore,
%         fsignal = 2*fsample
%      -> kmax = (2*pi/L)*(Nx/2)
%      -> kexp = Nxexp - Lexp - 1
%     Here, kmax =: 2^kexp; Nx =: 2^Nxexp; L =: 2*pi*2^Lexp.
%
Nxexp = Lexp + kexp + 1;
% method for time-marching
%   1 for ode45 func
%   2 for 4-order Runge-Kutta method
%   3 for ETDRK4 (Exponential Time-diff & 4-order R-K)
tmthd = 3;
% correction flag 
%   1 for using 3/2 rule
%   0 for direct computation
corrFlag = 0;
% construct the filename
filename = ...
  sprintf('%s','KSdata_sigma',num2str(s),'_dt',num2str(dt),...
               '_Lexp',num2str(Lexp),'_kexp',num2str(kexp),...
               '_Tmthd',num2str(tmthd));

%% Temporal and spatial discretization
% temporal grid number for iteration and plotting
Nt = T/dt;
Ntplt = T/dtplt;
% temporal discretization
t = [0:dt:dt*(Nt-1)].';
% ratio of time points for plotting
Rfplt = dtplt/dt;
% wave length L
L = (2^Lexp)*2*pi;
% spatial grid number
Nx = 2^Nxexp;
% spatial discretization
x = L/Nx*[-Nx/2:1:(Nx/2-1)].';
% discretization in spectual space
k = 2*pi/L*[0:1:(Nx/2-1) 0 (-Nx/2+1):1:-1].';
% zero padding in k-space (3/2 rule)
Zero = zeros(Nx/4,1,'like',1i);

%% initial condition
u0 = -sin(2^kexpsignal.*x);
ut0 = fft(u0);
