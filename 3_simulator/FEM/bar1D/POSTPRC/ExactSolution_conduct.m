% plots the exact stress  
function ExactSolution_conduct
include_flags; 
 
% divide the problem domain into two regions 
x = 0:0.01:20;

subplot(2,1,1); 

% exact displacement for xa 
u = -10*x.^2 + 400*x; 

% plot displacement 
plot(x,u, '--r' );       
legend('sdf','exact'); 

subplot(2,1,2); 

% exact stress for xa 
du = -20*x + 400;     

% plot stress 
plot(x,du, '--r' );
legend('FE','exact')
