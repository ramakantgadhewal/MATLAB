function [b_origin, b_opt, k0_origin, k0_opt] = CFD02_wave1D
% 1st-order 1D wave equation solver using FDM

%% include the global flags and input the data
 include_flags_CFD02;
 input_file_CFD02;

%% discrete the solution domain
% the discrete variable
u  = zeros(Nx+1,1);
x  = zeros(Nx+1,1);
t  = zeros(Nx+1,1);
for i = 1:Nx+1
    x(i,1) = x1 + dx*(i-1);
    u(i,1) = IC(x(i,1),Nphi,phik,k0,epsil,x1,x2,initialCondition_flag);
end
for i = 1:Nt+1
    t(i,1) = i*dt;
end
% the plot variable
leng_T = length(T);
leng_ind = length(ind);
u_plot = zeros(Nx+1,leng_T,leng_ind);
d_T = zeros(leng_T-1);
d_T(1) = T(1);
for i = 1:leng_T-1
    d_T(i+1) = T(i+1) - T(i);
end
%% iteration
tic
b_origin = 'null';
b_opt = 'null';
k0_origin = 'null';
k0_opt = 'null';
for IND = 1:leng_ind
    if ind(IND) <= 10   % Schemes where l,r is subscribed
        if ind(IND) == 1         % Up1Down0F scheme
            [u_plot(:,:,IND)] = cal_Up1Down0(u_plot(:,:,IND),u,leng_T,d_T);
        elseif ind(IND) == 2     % Lax-Wendroff scheme
            [u_plot(:,:,IND)] = cal_Lax_Wendroff(u_plot(:,:,IND),u,leng_T,d_T);
        elseif ind(IND) == 3     % Warming-Beam scheme
            [u_plot(:,:,IND)] = cal_Warming_Beam(u_plot(:,:,IND),u,leng_T,d_T);
        elseif ind(IND) == 4     % Up2Down1 scheme
            [u_plot(:,:,IND)] = cal_Up2Down1(u_plot(:,:,IND),u,leng_T,d_T);
        elseif ind(IND) == 5     % Up1Down12M scheme
            [u_plot(:,:,IND)] = cal_Up1Down2F_M(u_plot(:,:,IND),u,leng_T,d_T);
        end
    else            % UplDownrF schemes [Requirement: l > r]
        if ind(IND) == 11       % original scheme
            [~,~,b_origin,k0_origin] = calScheme_origin(l,r,e_res);
            [u_plot(:,:,IND)] = cal_UprDownlF(u_plot(:,:,IND),u,leng_T,d_T,l,r,b_origin);
        elseif ind(IND) == 12   % using Integration method (p-norm)
            [b_opt,k0_opt] = calScheme_optInt(l,r,prL,prR,nu,e_res);
            [u_plot(:,:,IND)] = cal_UprDownlF(u_plot(:,:,IND),u,leng_T,d_T,l,r,b_opt);
        elseif ind(IND) == 13   % using Taylor expansion method (l2 norm)
            [b_opt,k0_opt] = calScheme_optTaylor(l,r,prL,prR,nu,e_res);
            [u_plot(:,:,IND)] = cal_UprDownlF(u_plot(:,:,IND),u,leng_T,d_T,l,r,b_opt);
        end
    end
end
toc
%% plot the interesting timepoints
if timepointsPlot_flag == 1
    for IND = 1:leng_ind
        p = zeros(leng_T);
        figure(IND);
        for i = 1:leng_T
            p(i) = plot(x(:,1),u_plot(:,i,IND));
            hold on;
        end
        title(strcat(name_ind(ind(IND)),' scheme (IND=',num2str(ind(IND),' %u'),'), k_0=',num2str(k0,'%u')));
        legend(p(1:leng_T),...
            ['t = ',num2str(TIME(1),'%5.2f')],...
            ['t = ',num2str(TIME(2),'%5.2f')],...
            'Location','northwest');
        hold off;
    end
end
%% calculate the 1-norm error
if calError_flag == 1
    u_accu = zeros(Nx+1,1);
    T0 = T(leng_T);
    for i = 1:Nx+1
        x0 = x1 + dx*(i-1) - a*T0;
        u_accu(i,1) = IC(x0,Nphi,phik,k0,epsil,x1,x2,initialCondition_flag);
    end
    err1 = zeros(leng_ind,1);
    err2 = zeros(leng_ind,1);
    movement = floor(T0*mv/dt);
    for IND = 1:leng_ind
        for i = 1:Nx+1
            err1(IND) = err1(IND) + abs(u_accu(i,1) - u_plot(i,leng_T,IND));
        end
        for i = 1:Nx
            moveStep = mod(i-movement,Nx+1);
            if moveStep ~= 0
                err2(IND) = err2(IND) + abs(u_accu(i,1) - u_plot(moveStep,leng_T,IND));
            else
                err2(IND) = err2(IND) + abs(u_accu(i,1) - u_plot(Nx+1,leng_T,IND));
            end
        end
    end
    err1 = err1/(Nx+1);
    err2 = err2/(Nx+1);
    figure(IND+1);
    plot(ind,err1,'r*',ind,err2,'k+');
    xlabel('scheme index');
    ylabel('1-norm error');
    legend('original','modified')
    hold off;
end
Nx
dx
Nt
dt
l
r
prL
prR
end

function [f] = IC(x,Nphi,phik,k0,epsil,x1,x2,initialCondition_flag)
%% input the initial condition
f = 1;
if initialCondition_flag == 1   
    for k = 1:Nphi
        Ek = (k/k0)^4*exp(-2*(k/k0)^2);
        f = f + epsil*sqrt(Ek)*sin(2*pi*k*(x+phik(k)));
    end
elseif initialCondition_flag == 2
    if x<(x2-x1)/2
        f = 0;
    end
end
end