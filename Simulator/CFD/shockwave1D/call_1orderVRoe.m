function call_1orderVRoe
include_flags;
input_file;

%% Initialize the conservative and flux vector
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
A_abs = zeros(3,3,Nx);     % Absolute value of coefficient matrix A in C.V.
A_np = zeros(3,3,Nx);      % Coefficient matrix A in C.V.

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
        % NOTICE: E_np != gam * (H_np - 0.5 * u_np^2) + 0.5 * u_np^2
        Lambda(1,1) = u_np - c_np;
        Lambda(2,2) = u_np;
        Lambda(3,3) = u_np + c_np;
        A_np(:,:,i) = cal_AH(u_np, H_np, gam);
        R_np = cal_R(u_np, c_np, H_np);
        L_np = cal_L(u_np, c_np, gam);
        Lambda = abs(Lambda);
        for k = 1:3
            if Lambda(k,k) < fixEnt
                Lambda(k,k) = 0.5*(Lambda(k,k)^2 + fixEnt^2)/fixEnt;
            end
        end
        A_abs(:,:,i) = R_np * Lambda * L_np;
    end
    % Internal points
    for i = 2:Nx
        F_right = 0.5 * (A_np(:,:,i) * (U(:,i+1) + U(:,i)) - A_abs(:,:,i) * (U(:,i+1) - U(:,i)));
        F_left = 0.5 * (A_np(:,:,i-1) * (U(:,i) + U(:,i-1)) - A_abs(:,:,i-1) * (U(:,i) - U(:,i-1)));
        Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    end
    % Boundary points
    Up(:,1) = U(:,1);
    Up(:,Nx+1) = U(:,Nx+1);

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