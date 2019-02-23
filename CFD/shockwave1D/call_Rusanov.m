function call_Rusanov
include_flags;
input_file;

%% Initialize the flux vector and coefficients
U = zeros(3,Nx+1);
Up = zeros(3,Nx+1);
F = zeros(3,Nx+1);
lambda = zeros(Nx+1);
e = 0.5;

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
% Iteration for (rho, m, epsi)
    for i = 2:Nx
        lambda_right = max(lambda(i),lambda(i+1));
        lambda_left = max(lambda(i-1),lambda(i));
%         lambda_right = 0.5 * (lambda(i) + lambda(i+1));
%         lambda_left = 0.5 * (lambda(i-1) + lambda(i));
        F_right = 1/2 * (F(:,i+1) + F(:,i)) - lambda_right * e * (U(:,i+1) - U(:,i));
        F_left  = 1/2 * (F(:,i) + F(:,i-1)) - lambda_left * e * (U(:,i) - U(:,i-1));
        Up(:,i) = U(:,i) - dt/dx*(F_right - F_left);
    end
    Up(:,1) = U(:,1);
    Up(:,Nx+1) = U(:,Nx+1);
% Transition back to the form of conservative variables
    rho(:,j+1) = Up(1,:);
    m(:,j+1) = Up(2,:);
    epsi(:,j+1) = Up(3,:);
end 
%% calculate the fundamental variables
u(:,Nt+1) = m(:,Nt+1)./rho(:,Nt+1);
p(:,Nt+1) = (gam - 1).*(epsi(:,Nt+1) - (rho(:,Nt+1).*u(:,Nt+1).^2)/2);

end