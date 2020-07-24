function call_2orderFVS
include_flags;
input_file;

%% Initialize the conservative and flux vector
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
F_posi = zeros(3,Nx+1);
F_nega = zeros(3,Nx+1);

%% Loop for each timestep
for j = 1:Nt
% Transition to the conservative vector form
    U(1,:) = rho(:,j);
    U(2,:) = m(:,j);
    U(3,:) = epsi(:,j);
% Iteration for (rho, m, epsi)
    % Update flux variables
    for i = 1:Nx+1
        % Update the fundamental variables
        u(i,j) = U(2,i)/U(1,i);
        p(i,j) = (gam - 1)*(U(3,i) - (U(2,i)^2)/U(1,i)/2);
        % Update coefficient matrix per node
        A = cal_A(u(i,j),U(3,i)/U(1,i),gam);
        % Split the flux vector
        [R, Lambda] = eig(A);
%         max_lambda = max(max(Lambda));
%         Lambda_posi = (Lambda + max_lambda*eye(3))/2;
%         Lambda_nega = (Lambda - max_lambda*eye(3))/2;
        Lambda_posi = (Lambda + abs(Lambda))/2;
        Lambda_nega = (Lambda - abs(Lambda))/2;
        A_posi = R * Lambda_posi / R;
        A_nega = R * Lambda_nega / R;
        F_posi(:,i) = A_posi * U(:,i);
        F_nega(:,i) = A_nega * U(:,i);
    end
    % Internal points
    for i = 3:Nx-1
        F_right = 1.5*F_posi(:,i) - 0.5*F_posi(:,i-1) + 1.5*F_nega(:,i+1) - 0.5*F_nega(:,i+2);
        F_left = 1.5*F_posi(:,i-1) - 0.5*F_posi(:,i-2) + 1.5*F_nega(:,i) - 0.5*F_nega(:,i+1);
        Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    end
    % Boundary points
    Up(:,1) = U(:,1);
    Up(:,Nx+1) = U(:,Nx+1);
    
    i = 2;
    F_right = F_posi(:,i) + F_nega(:,i+1);
    F_left = F_posi(:,i-1) + F_nega(:,i);
    Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    
    i = Nx;
    F_right = F_posi(:,i) + F_nega(:,i+1);
    F_left = F_posi(:,i-1) + F_nega(:,i);
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