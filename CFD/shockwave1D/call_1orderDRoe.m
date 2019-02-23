function call_1orderDRoe
include_flags;
input_file;

%% Initialize the conservative and flux vector
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
A_cv = zeros(3,3,Nx);      % Coefficient matrix A in C.V.
A_sign = zeros(3,3,Nx); % sign of characteristic matrix of A in C.V.
F = zeros(3,Nx+1);

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
        F(1,i) = m(i,j);
        F(2,i) = (m(i,j)^2)/rho(i,j) + (gam - 1)*(epsi(i,j) - (m(i,j)^2)/rho(i,j)/2);
        F(3,i) = m(i,j)/rho(i,j)*(epsi(i,j) + (gam - 1)*(epsi(i,j) - (m(i,j)^2)/rho(i,j)/2));
    end
    % Update Roe mean value
    % U(1)  A_ctrl(1)  U(2)  A_ctrl(2)  ...  U(Nx)  A_ctrl(Nx)  U(Nx+1)
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
        u_ctrl = (sqrt(rhoL)*uL + sqrt(rhoR)*uR) / (sqrt(rhoL) + sqrt(rhoR));
        H_ctrl = (sqrt(rhoL)*HL + sqrt(rhoR)*HR) / (sqrt(rhoL) + sqrt(rhoR));
        e_ctrl = gam * (H_ctrl - 0.5 * u_ctrl^2);
        E_ctrl = e_ctrl + 0.5 * u_ctrl^2;
        c_ctrl = sqrt((gam - 1) * (H_ctrl - 0.5 * u_ctrl^2));
        % Determine the viscosity term
        Lambda(1,1) = u_ctrl - c_ctrl;
        Lambda(2,2) = u_ctrl;
        Lambda(3,3) = u_ctrl + c_ctrl;
        A_ctrl = cal_A(u_ctrl, E_ctrl, gam);
        R_ctrl = cal_R(u_ctrl, c_ctrl, H_ctrl);
        L_ctrl = cal_L(u_ctrl, c_ctrl, gam);
        % Determine the viscosity term
        A_cv(:,:,i) = A_ctrl;
        A_sign(:,:,i) = R_ctrl * sign(Lambda) * L_ctrl;
    end
    % Internal points
    for i = 2:Nx
        % without entropy fixing
        % HOWEVER, FR - FL is NOT equal to A*(UR - UL), which is
        %   contradictory with the conservative condition!
        F_right = 0.5 * (F(:,i) + F(:,i+1) - A_sign(:,:,i) * (F(:,i+1) - F(:,i)));
        F_left = 0.5 * (F(:,i-1) + F(:,i) - A_sign(:,:,i-1) * (F(:,i) - F(:,i-1)));
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