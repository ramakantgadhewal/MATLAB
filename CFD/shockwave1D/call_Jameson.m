function call_Jameson
include_flags;
input_file;

%% Initialize the flux vector and coefficients
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
F = zeros(3,Nx+1);
lambda = zeros(Nx+1);
e2 = zeros(Nx+1);
e4 = zeros(Nx+1);
nu = zeros(Nx+5);

k2 = 10;
k4 = 1/64;

%% Loop for each timestep
for j = 1:Nt
% Transition to the conservative vector form
    U(1,:) = rho(:,j);
    U(2,:) = m(:,j);
    U(3,:) = epsi(:,j);
% Update flux variables
    for i = 1:Nx+1
        % Update the fundamental variables
        u(i,j) = m(i,j)/rho(i,j);
        p(i,j) = (gam - 1)*(epsi(i,j) - (m(i,j)^2)/rho(i,j)/2);
        F(1,i) = m(i,j);
        F(2,i) = (m(i,j)^2)/rho(i,j) + (gam - 1)*(epsi(i,j) - (m(i,j)^2)/rho(i,j)/2);
        F(3,i) = m(i,j)/rho(i,j)*(epsi(i,j) + (gam - 1)*(epsi(i,j) - (m(i,j)^2)/rho(i,j)/2));
        lambda(i) = abs(u(i,j)) + sqrt(gam*p(i,j)/rho(i,j));   
    end
% Determine e2 and e4 
    % nu: 1  2  3  4  5  ... Nx+2 Nx+3 Nx+4 Nx+5
    % p: -1  0 [1  2  3  ...  Nx  Nx+1]Nx+2 Nx+3
    % v:     0  1  2  3  ... 
    % i:        1  2  3  ...  Nx  Nx+1 Nx+2 Nx+3
    nu(1) = abs((p(2,j)-p(1,j))/(p(2,j)+p(1,j)));
    nu(2) = nu(1);
    nu(3) = nu(1);
    for i = 2:Nx
        num_nu = i + 2;
        nu(num_nu) = abs((p(i+1,j) - 2*p(i,j) + p(i-1,j))/(p(i+1,j) + 2*p(i,j) + p(i-1,j)));
    end
    nu(Nx+3) = nu(Nx+2);
    nu(Nx+4) = nu(Nx+2);
    nu(Nx+5) = nu(Nx+2);
    for i = 1:Nx+1
        num_nu = i + 2;
        e2(i) = k2*max(nu(num_nu-1),max(nu(num_nu),max(nu(num_nu+1),nu(num_nu+2))));
        e4(i) = max(0, k4-e2(i));
    end
% Iteration for (rho, m, epsi)
    for i = 3:Nx-1
        lambda_right = max(lambda(i),lambda(i+1));
        lambda_left = max(lambda(i-1),lambda(i));
%         lambda_right = 0.5 * (lambda(i) + lambda(i+1));
%         lambda_left = 0.5 * (lambda(i-1) + lambda(i));
        F_right = 1/2 * (F(:,i+1) + F(:,i)) - lambda_right * e2(i) * (U(:,i+1) - U(:,i)) ...
            + lambda_right * e4(i) * (F(:,i+2) - 3*F(:,i+1) + 3*F(:,i) - F(:,i-1));
        F_left  = 1/2 * (F(:,i) + F(:,i-1)) - lambda_left * e2(i-1) * (U(:,i) - U(:,i-1)) ...
            + lambda_right * e4(i-1) * (F(:,i+1) - 3*F(:,i) + 3*F(:,i-1) - F(:,i-2));
        Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    end
% Deal with boundary points
    Up(:,1) = U(:,1);
    Up(:,Nx+1) = U(:,Nx+1);
    
    i = 2;    
    lambda_right = max(lambda(i),lambda(i+1));
    lambda_left = max(lambda(i-1),lambda(i));
%         lambda_right = 0.5 * (lambda(i) + lambda(i+1));
%         lambda_left = 0.5 * (lambda(i-1) + lambda(i));
    F_right = 1/2 * (F(:,i+1) + F(:,i)) - lambda_right * e2(i) * (U(:,i+1) - U(:,i));
    F_left  = 1/2 * (F(:,i) + F(:,i-1)) - lambda_left * e2(i-1) * (U(:,i) - U(:,i-1));
    Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);

    i = Nx;
    lambda_right = max(lambda(i),lambda(i+1));
    lambda_left = max(lambda(i-1),lambda(i));
%         lambda_right = 0.5 * (lambda(i) + lambda(i+1));
%         lambda_left = 0.5 * (lambda(i-1) + lambda(i));
    F_right = 1/2 * (F(:,i+1) + F(:,i)) - lambda_right * e2(i) * (U(:,i+1) - U(:,i));
    F_left  = 1/2 * (F(:,i) + F(:,i-1)) - lambda_left * e2(i-1) * (U(:,i) - U(:,i-1));
    Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);

% Transition back to the form of conservative variables
    rho(:,j+1) = Up(1,:);
    m(:,j+1) = Up(2,:);
    epsi(:,j+1) = Up(3,:);
end
%% calculate the fundamental variables
u(:,Nt+1) = m(:,Nt+1)./rho(:,Nt+1);
p(:,Nt+1) = (gam - 1).*(epsi(:,Nt+1) - (rho(:,Nt+1).*u(:,Nt+1).^2)/2);

% realize the fundamental variables
rho = real(rho);
u = real(u);
p = real(p);
end