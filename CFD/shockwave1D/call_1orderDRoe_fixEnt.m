function call_1orderDRoe_fixEnt
include_flags;
input_file;

%% Initialize the conservative and flux vector
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
A_abs = zeros(3,3,Nx);  % Absolute value of coefficient matrix A in C.V.
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
        eL = U(3,i)/rhoL;
        eR = U(3,i+1)/rhoR;
        % Calculate necessary characteristic variables
        EL = eL + 0.5*uL^2;
        ER = eR + 0.5*uR^2;
        u_ctrl = (sqrt(rhoL)*uL + sqrt(rhoR)*uR) / (sqrt(rhoL) + sqrt(rhoR));
        E_ctrl = (sqrt(rhoL)*EL + sqrt(rhoR)*ER) / (sqrt(rhoL) + sqrt(rhoR));
        A_ctrl = cal_A(u_ctrl, E_ctrl, gam);
        [R, Lambda] = eig(A_ctrl);
        % Determine the viscosity term
        Lambda = abs(Lambda);
        for k = 1:3
            if Lambda(k,k) < fixEnt
                Lambda(k,k) = 0.5*(Lambda(k,k)^2 + fixEnt^2)/fixEnt;
            end
        end
        A_abs(:,:,i) = R * Lambda / R;
    end
    % Internal points
    for i = 2:Nx
        % with entropy fixing
        F_right = 0.5 * (F(:,i) + F(:,i+1) - A_abs(:,:,i) * (U(:,i+1) - U(:,i)));
        F_left = 0.5 * (F(:,i-1) + F(:,i) - A_abs(:,:,i-1) * (U(:,i) - U(:,i-1)));
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