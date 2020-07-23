function test_vorStr_FTCS_GS
%% FTCS, Gauss-Seidel
  % domain definition
  CFL = 0.1;     
  Re = 1000;       
  Nx = 128;
  Ny = 128;

  % iteration parameters
  err_uMag_thre = 1e-6;
  err_psif_thre = 1e-6;
  nIter_max = 1e5;
  outputInt = 100;
  outputName = 'test';
  cal_p_flag = 0;

  % derived quantities
  uMag = 1;          % in m/s
  lMag = 1;          % in m
  Lx = lMag;
  Ly = lMag;
  dx = Lx/Nx;
  dy = Ly/Ny;
  dt = CFL * min(dx,dy)/uMag;  % CFL := uMag*dt/dx

  %% Initialization
  x = 0:dx:Lx;
  y = 0:dy:Ly;
  U = zeros(Nx+1,Ny+1);     % x-velocity
  V = zeros(Nx+1,Ny+1);     % y-velocity
  w = zeros(Nx+1,Ny+1);     % vorticity
  psif = zeros(Nx+1,Ny+1);  % stream function
  U(2:Nx,Ny+1) = uMag;
  
  % Calculate (u,v)
  err_uMag = err_uMag_thre + 1;
  nIter    = 0;
  while err_uMag > err_uMag_thre && nIter < nIter_max
    U0 = U; V0 = V;
    nIter = nIter + 1;
    
    w(2:Nx,1)    = 1/dy/dy*2*psif(2:Nx,2);
    w(2:Nx,Ny+1) = 1/dy/dy*2*(psif(2:Nx,Ny) + uMag*dy);
    w(1,2:Ny)    = 1/dx/dx*2*psif(2,2:Ny);
    w(Nx+1,2:Ny) = 1/dx/dx*2*psif(Nx,2:Ny);
    w(1,1)       = (w(1,2) + w(2,1))/2;
    w(1,Ny+1)    = (w(1,Ny) + w(2,Ny+1))/2;
    w(Nx+1,Ny+1) = (w(Nx,Ny+1) + w(Nx+1,Ny))/2;
    w(Nx+1,1)    = (w(Nx,1) + w(Nx+1,2))/2;
    
    w(2:Nx,2:Ny) = w(2:Nx,2:Ny) ...
      - U(2:Nx,2:Ny) .* (w(3:Nx+1,2:Ny) - w(1:Nx-1,2:Ny)) * dt/2/dx ...
      - V(2:Nx,2:Ny) .* (w(2:Nx,3:Ny+1) - w(2:Nx,1:Ny-1)) * dt/2/dy ...
      + 1/Re * dt/dx/dx * (w(3:Nx+1,2:Ny) - 2*w(2:Nx,2:Ny) + w(1:Nx-1,2:Ny)) ...
      + 1/Re * dt/dy/dy * (w(2:Nx,3:Ny+1) - 2*w(2:Nx,2:Ny) + w(2:Nx,1:Ny-1));
    
    err_psif = err_psif_thre + 1;
    while (err_psif > err_psif_thre)
      psif0 = psif;
      psif(2:Nx,2:Ny) = ((psif(3:Nx+1,2:Ny) + psif(1:Nx-1,2:Ny))/dx/dx ...
                  + (psif(2:Nx,3:Ny+1) + psif(2:Nx,1:Ny-1))/dy/dy ...
                  - w(2:Nx,2:Ny)) / (2*(1/dx/dx + 1/dy/dy));
      err_psif = max(max(abs(psif-psif0)));
    end
  
    U(2:Nx,2:Ny) = (psif(2:Nx,3:Ny+1) - psif(2:Nx,1:Ny-1))/2/dy;
    V(2:Nx,2:Ny) = (psif(1:Nx-1,2:Ny) - psif(3:Nx+1,2:Ny))/2/dx;
    
    err_uMag = max(max(max(abs(U-U0))),max(max(abs(V-V0))));
    if (mod(nIter,outputInt) == 0)
      fprintf('%s%i%s%e\n','The ',nIter,'th step gets accuracy = ',err_uMag);
    end
  end

  % Calculate p
  if(cal_p_flag == 1)
    p = cal_p(U,V);
  end

  % Plot the flow field
  [Y,X] = meshgrid(y,x);
  contour(X,Y,psif,50);
  Q = quiver(X,Y,U,V);
  Q.AutoScaleFactor = 1.2;
  
  [X,Y]=meshgrid(y,x);
  U = U';
  V = V';
  save([outputName '.mat'],'X','Y','U','V')

end