function call_2orderVTVD
include_flags;
input_file;

%% Initialize the conservative and flux vector
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
A_abs = zeros(3,3,Nx);     % Absolute value of coefficient matrix A in C.V.
A_np = zeros(3,3,Nx);      % Coefficient matrix A on nodal points
R_np = zeros(3,3,Nx);
L_np = zeros(3,3,Nx);
Lambda = zeros(3,3,Nx);
WL = zeros(3,1);
WR = zeros(3,1);

%% Loop for each timestep
for j = 1:Nt
% Transition to the conservative vector form
    U(1,:) = rho(:,j);
    U(2,:) = m(:,j);
    U(3,:) = epsi(:,j);
% Iteration for (rho, m, epsi)
    % Update flux variablesv
    for i = 1:Nx+1
        % Update the fundamental variables
        u(i,j) = m(i,j)/rho(i,j);
        p(i,j) = (gam - 1)*(epsi(i,j) - (m(i,j)^2)/rho(i,j)/2);
    end
    % Update Roe mean value
    % U_ctrl(1)  A_np(1)  U_ctrl(2)  A_np(2)  ...  U_ctrl(Nx)  A_np(Nx)  U_ctrl(Nx+1)
    for i = 1:Nx
        % Get local flow info
        rhoL = U(1,i);
        rhoR = U(1,i+1);
        uL = U(2,i)/rhoL;
        uR = U(2,i+1)/rhoR;
        EL = U(3,i)/rhoL;
        ER = U(3,i+1)/rhoR;
        HL = EL + (gam - 1) * (EL - 0.5 * uL^2);
        HR = ER + (gam - 1) * (ER - 0.5 * uR^2);
        u_np = (sqrt(rhoL)*uL + sqrt(rhoR)*uR) / (sqrt(rhoL) + sqrt(rhoR));
        H_np = (sqrt(rhoL)*HL + sqrt(rhoR)*HR) / (sqrt(rhoL) + sqrt(rhoR));
        c_np = sqrt((gam - 1) * (H_np - 0.5 * u_np^2));
        % Determine the viscosity term
        Lambda(1,1,i) = u_np - c_np;
        Lambda(2,2,i) = u_np;
        Lambda(3,3,i) = u_np + c_np;
        A_np(:,:,i) = cal_AH(u_np, H_np, gam);
        R_np(:,:,i) = cal_R(u_np, c_np, H_np);
        L_np(:,:,i) = cal_L(u_np, c_np, gam);
        % Determine the viscosity term
        Lambda_abs = abs(Lambda(:,:,i));
        for k = 1:3
            if Lambda_abs(k,k) < fixEnt
                Lambda_abs(k,k) = 0.5*(Lambda_abs(k,k)^2 + fixEnt^2)/fixEnt;
            end
        end
        A_abs(:,:,i) = R_np(:,:,i) * Lambda_abs * L_np(:,:,i);
    end
    % Internal points
    for i = 3:Nx-1
        for k = 1:3
            DL = minmod(  L_np(k,:,i)*(U(:,i+1) - U(:,i)),     L_np(k,:,i)*(U(:,i) - U(:,i-1)) ) /dx;
            WL(k,1) = L_np(k,:,i) * U(:,i) + 0.5 * DL * (dx - Lambda(k,k,i) * dt);
            DR = minmod(  L_np(k,:,i)*(U(:,i+2) - U(:,i+1)),   L_np(k,:,i)*(U(:,i+1) - U(:,i)) ) /dx;
            WR(k,1) = L_np(k,:,i) * U(:,i+1) - 0.5 * DR * (dx + Lambda(k,k,i) * dt);
        end
        U_right_R = R_np(:,:,i) * WR;
        U_right_L = R_np(:,:,i) * WL;
        
        for k = 1:3
            DL = minmod(  L_np(k,:,i-1)*(U(:,i) - U(:,i-1)),  L_np(k,:,i-1)*(U(:,i-1) - U(:,i-2)) ) /dx;
            WL(k,1) = L_np(k,:,i-1) * U(:,i-1) + 0.5 * DL * (dx - Lambda(k,k,i-1) * dt);
            DR = minmod(  L_np(k,:,i-1)*(U(:,i+1) - U(:,i)),  L_np(k,:,i-1)*(U(:,i) - U(:,i-1))   ) /dx;
            WR(k,1) = L_np(k,:,i-1) * U(:,i) - 0.5 * DR * (dx + Lambda(k,k,i-1) * dt);
        end
        U_left_R = R_np(:,:,i-1) * WR;
        U_left_L = R_np(:,:,i-1) * WL;
        
        F_right = 0.5 * (A_np(:,:,i) * (U_right_R + U_right_L) - A_abs(:,:,i) * (U_right_R - U_right_L));
        F_left = 0.5 * (A_np(:,:,i-1) * (U_left_R + U_left_L) - A_abs(:,:,i-1) * (U_left_R - U_left_L));
        Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    end
    % Boundary points
    Up(:,1) = U(:,1);
    Up(:,Nx+1) = U(:,Nx+1);
    i = 2;
    F_right = 0.5 * (A_np(:,:,i) * (U(:,i+1) + U(:,i)) - A_abs(:,:,i) * (U(:,i+1) - U(:,i)));
    F_left = 0.5 * (A_np(:,:,i-1) * (U(:,i) + U(:,i-1)) - A_abs(:,:,i-1) * (U(:,i) - U(:,i-1)));
    Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    i = Nx;
    F_right = 0.5 * (A_np(:,:,i) * (U(:,i+1) + U(:,i)) - A_abs(:,:,i) * (U(:,i+1) - U(:,i)));
    F_left = 0.5 * (A_np(:,:,i-1) * (U(:,i) + U(:,i-1)) - A_abs(:,:,i-1) * (U(:,i) - U(:,i-1)));
    Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    
%% % Transition back to the form of conservative variables
    rho(:,j+1) = Up(1,:);
    m(:,j+1) = Up(2,:);
    epsi(:,j+1) = Up(3,:);
end

%% calculate the fundamental variables
u(:,Nt+1) = m(:,Nt+1)./rho(:,Nt+1);
p(:,Nt+1) = (gam - 1).*(epsi(:,Nt+1) - (m(:,Nt+1).^2)./rho(:,Nt+1)./2);

% realize the fundamental variables
rho = real(rho);
u = real(u);
p = real(p);
end

function [M] = minmod(a,b)
    if a*b > 0
        M = sign(a) * min(abs(a), abs(b));
    else
        M = 0;
    end
end