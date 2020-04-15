%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poiseuilleFlow.m: Lattice Boltzmann sample in Matlab      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified: Yunfan Huang
%   Original: WangLab
%   Referred: Palabos (cylinder.m)
% Note.
%   Relaxation method: SRT
%   Solid boundary: Fullway bounce-back (link-wise)
%   Open boundary: NEBB/Zou-He B.C. (wet-node)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close('all');
clear;

%% SYSTEM SETUP
% Space and time domain
lx  = 600;
ly  = 400;
[x,y] = meshgrid(1:lx,1:ly);  % get coordinate of matrix indices
obst_x = lx/5+1;   % x position of the cylinder
obst_y = ly/2+3;   % y position of the cylinder
% obst_r = 0;      % radius of the cylinder
obst_r = 11;       % radius of the cylinder
maxT  = 5000;      % total iterations
tPlot = 100;       % cycles for plot
% Physical condition
uMax = 0.2;  % maximum velocity of Poiseuille inflow (non-dim)
Re   = 200;   % Reynolds number
nu   = uMax*2*obst_r/Re; % kinematic viscosity (non-dim)
tau  = 3*nu + 1/2;       % relaxation time (non-dim)
g    = 0.0;  % body force

% Open boundary
in = 1;   % position of inlet
out = lx; % position of outlet
col = [2:ly-1];% position of the fluids
width = ly - 2;  % width of the channel
% Solid boundary
is_solid_node = ...% location of cylinder
    (x-obst_x).^2 + (y-obst_y).^2 < obst_r.^2;
is_solid_node([1,ly],:) = 1;
    % A1,B: is_solid_node([1,ly],:) = 1;
    % A2,C: is_solid_node([1,ly],:) = 0;
[y_solid_node,x_solid_node] = find(is_solid_node);
N_solid_node = length(x_solid_node);

figure(10); imagesc(is_solid_node); axis('equal','tight','off');
xlabel('400'); ylabel('100');

% D2Q9 lattice constant
w   = [ 4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
cx  = [   0,   1,   0,  -1,   0,    1,   -1,   -1,    1];
cy  = [   0,   0,   1,   0,  -1,    1,    1,   -1,   -1];
opp = [   1,   4,   5,   2,   3,    8,    9,    6,    7];

%% INITIAL CONDITION: Poiseuille profile at equilibrium
y_phys = y-1.5;
u_x = 4 * uMax / (width*width) * (y_phys.*width-y_phys.*y_phys);
    % A1,A2: u_x = 4 * uMax / (width*width) * (y_phys.*width-y_phys.*y_phys);
    % B,C: u_x = uMax*ones(ly,lx);
u_y = zeros(ly,lx);
rho = ones(ly,lx);
fIn = zeros(ly,lx,9);
for i=1:9
    cu = 3*(cx(i)*u_x+cy(i)*u_y);
    fIn(:,:,i) = rho .* w(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 3/2*(u_x.^2+u_y.^2) );
end
fEq = fIn;
fOut = fIn;

%% MAIN LOOP
for cycle = 1:maxT
    % MACROSCOPIC VARIABLES (density and velocity)
    rho = sum(fIn,3);
    u_x = reshape ( (reshape(fIn,lx*ly,9))*cx(:), ly,lx) ./rho;
    u_y = reshape ( (reshape(fIn,lx*ly,9))*cy(:), ly,lx) ./rho;
    
    % MACROSCOPIC BOUNDARY CONDITIONS: Dirichlet condition
      % Inlet: Poiseuille profile
    y_phys = col - 1.5;
    u_x(col,in) = 4 * uMax/(width*width)*(y_phys.*width-y_phys.*y_phys);
      % A1,A2: u_x(col,in) = 4 * uMax/(width*width)*(y_phys.*width-y_phys.*y_phys);
      % B,C: u_x(col,in) = uMax*ones(length(col),1);
    u_y(col,in) = 0;
    rho(col,in) = 1 ./ (1-u_x(col,in)) .* ...
      ( sum(fIn(col,in,[1,3,5]),3) + 2*sum(fIn(col,in,[4,7,8]),3) );
      % Outlet: Constant pressure
    rho(col,out) = 1;
    u_x(col,out) = -1 + 1 ./ (rho(col,out)) .* ( ...
        sum(fIn(col,out,[1,3,5]),3) + 2*sum(fIn(col,out,[2,6,9]),3) );
    u_y(col,out) = 0;
    
    % MICROSCOPIC BOUNDARY CONDITIONS: Zou/He B.C.
      % Inlet
    fIn(col,in,2) = fIn(col,in,4) + 2/3*rho(col,in).*u_x(col,in); 
    fIn(col,in,6) = fIn(col,in,8) + 1/2*(fIn(col,in,5)-fIn(col,in,3)) ... 
                                    + 1/2*rho(col,in).*u_y(col,in) ...
                                    + 1/6*rho(col,in).*u_x(col,in); 
    fIn(col,in,9) = fIn(col,in,7) + 1/2*(fIn(col,in,3)-fIn(col,in,5)) ... 
                                    - 1/2*rho(col,in).*u_y(col,in) ...
                                    + 1/6*rho(col,in).*u_x(col,in); 
      % Outlet
    fIn(col,out,4) = fIn(col,out,2) - 2/3*rho(col,out).*u_x(col,out); 
    fIn(col,out,8) = fIn(col,out,6) + 1/2*(fIn(col,out,3)-fIn(col,out,5)) ... 
                                      - 1/2*rho(col,out).*u_y(col,out) ...
                                      - 1/6*rho(col,out).*u_x(col,out); 
    fIn(col,out,7) = fIn(col,out,9) + 1/2*(fIn(col,out,5)-fIn(col,out,3)) ... 
                                      + 1/2*rho(col,out).*u_y(col,out) ...
                                      - 1/6*rho(col,out).*u_x(col,out); 
    
    % COLLISION STEP 
    u_xEq = u_x + tau*g;
    u_yEq = u_y;
    for i = 1:9  % Regular collision
      cu = 3*(cx(i)*u_xEq + cy(i)*u_yEq);
      fEq(:,:,i) = rho .* w(i) .* ...
                  ( 1 + cu + 1/2*(cu.*cu) - 3/2*(u_xEq.^2+u_yEq.^2) );
      fOut(:,:,i) = fIn(:,:,i) - (fIn(:,:,i) - fEq(:,:,i))/tau;
    end
    for i = 1:9  % Fullway bounce-back
      for j = 1:N_solid_node
        fOut(y_solid_node(j),x_solid_node(j),i) = ...
         fIn(y_solid_node(j),x_solid_node(j),opp(i));
      end
    end
        
    % SREAMING STEP
    for i=1:9
      fIn(:,:,i) = circshift(fOut(:,:,i), [cy(i),cx(i),0]);
    end
    
    % VISUALIZATION
    if (mod(cycle,tPlot)==1)
      cycle
      u = reshape(sqrt(u_x.^2+u_y.^2),ly,lx);
      u(is_solid_node) = nan;
      figure(1);
      imagesc(u); colorbar; axis equal tight off; drawnow; hold on;
%       quiver(x,y,10*u_x,10*u_y); axis equal tight; drawnow; hold on;
    end
end

disp('Simulation Over!');
% 
% %% VERIFICATION
% figure(100)
% y_phys = y(:,1)-1.5;
% % Model resulcycle as red circles
% plot(y_phys-width/2,u_x(:,1),'ro') 
% hold on
% % Poiseuille velocity profile as blue line
% nu=1/3*(tau-1/2);
% plot(y_phys-width/2,g/(2*nu)*((width/2)^2-(y_phys-width/2).^2));
