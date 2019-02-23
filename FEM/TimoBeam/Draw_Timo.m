function Draw_Timo()
 CAL_ERROR = 0;
 CONV_ANAY = 0;
 PATCHTEST = 1;
 v0 = [0 2.3148148148e-02 1.8148148148e-01 3.7037037037e-01];
 e0 = [0 1.8518518519e-01 5.1851851852e-01 7.4074074074e-01];
%% convergence analysis
if CONV_ANAY == 1
    h = [1 2 4 8];
    err_M_0 = [1.371739e+00 3.429348e-01 8.573371e-02 2.143345e-02];
    err_M_1 = [2.930884e-01 1.795422e-01 7.001332e-02 2.030200e-02];
    err_F_0 = [3.740312e-01 1.151161e-01 3.029294e-02 7.829946e-03];
    err_F_1 = [1.043677e-01 6.405428e-02 2.513301e-02 7.454413e-03];
    h = log(h);
    err_M_0 = log(err_M_0);
    err_M_1 = log(err_M_1);
    err_F_0 = log(err_F_0);
    err_F_1 = log(err_F_1);
    p1 = polyfit(h(2:4),err_M_0(2:4),1);
    p2 = polyfit(h(2:4),err_M_1(2:4),1);
    p3 = polyfit(h(2:4),err_F_0(2:4),1);
    p4 = polyfit(h(2:4),err_F_1(2:4),1);
    figure(1)
        subplot(1,2,1)
        plot(h,polyval(p1,h),'b',h,err_M_0,'r+');
        xlabel('log h');
        ylabel('log err');
        legend(['\alpha = ',num2str(p1(1))]);
        title('Concentrated M_z (with SRINT)');
        subplot(1,2,2);
        plot(h,polyval(p2,h),'b',h,err_M_1,'r+');
        xlabel('log h');
        ylabel('log err');
        legend(['\alpha = ', num2str(p2(1))]);
        title('Concentrated M_z (without SRINT)');
    figure(2)
        subplot(1,2,1)
        plot(h,polyval(p3,h),'b',h,err_F_0,'r+');
        xlabel('log h');
        ylabel('log err');
        legend(['\alpha = ',num2str(p3(1))]);
        title('Concentrated F_Q (with SRINT)')
        subplot(1,2,2);
        plot(h,polyval(p4,h),'b',h,err_F_1,'r+');
        xlabel('log h');
        ylabel('log err');
        legend(['\alpha = ',num2str(p4(1))]);
        title('Concentrated F_Q (without SRINT)');
end
if CAL_ERROR == 1
% Number of elements
    N = 1;
% Number of Gauss points
    n = 3;
    [weight,GaussPoints] = GetGauss(n);
    len = 1/N;
% Calculate the err of energy-norm
    INTk = 0;
    INTy = 0;
    for i = 1:N
        J = 0.5 * len;
        xe1 = (i-1)/N;
        xe2 = i/N;
        ve1 = v0(i);
        ve2 = v0(i+1);
        ee1 = e0(i);
        ee2 = e0(i+1);
        FEMdv = (ve2 - ve1)/len;
        FEMde = (ee2 - ee1)/len;
        Intk = 0;
        Inty = 0;
        for j = 1:n
          % Calculate the Gauss points in the physical coordinate
            gaussj = 0.5*(xe1 + xe2) + 0.5*GaussPoints(j)*(xe2 - xe1);
            FEMej = Linear(xe1,ee1,xe2,ee2,gaussj);
            [EXAvj,EXAej] = Fsolution(gaussj);
            [EXAdvj,EXAdej] = Fdsolution(gaussj);
          % Do the Gauss integration
            Intk = Intk + J*weight(j)*(EXAdej - FEMde)^2;
            Inty = Inty + J*weight(j)*((FEMdv - FEMej) - (EXAdvj - EXAej))^2;
        end
        INTk = INTk + 0.5*1000*0.00135*Intk;
        INTy = INTy + 0.5*400*0.18/1.2*Inty;
    end
    err = INTk + INTy;
    fprintf('INTk = %d \n',INTk);
    fprintf('INTy = %d \n',INTy);
    fprintf('Error = %d \n',err);
end
%% patch test
if PATCHTEST == 1
    x0 = [0 0.25 0.7 1];
% Calculate the disp and rota
  % Number of points plotting the accurate solution
    N = 100;
    x = zeros(N+1);
    v = zeros(N+1);
    e = zeros(N+1);
    for i = 1:N+1
     x(i) = (i-1)/N;
     [v(i),e(i)] = Msolution(x(i));
    end
  % Plot the disp and rota
    figure(1)
    ax1 = subplot(1,2,1);
    plot(ax1,x0,v0,'r+',x,v,'b-');
    legend(ax1,'FEM solution','Exact solution','Location','southeast');
    title(ax1,'Displacement')
    ax2 = subplot(1,2,2);
    plot(ax2,x0,e0,'r+',x,e,'b-');
    title(ax2,'Rotation')
    legend(ax2,'FEM solution','Exact solution','Location','southeast');
% % Calculate the bending and shear strain
%     x00 = [0.125, 0.375, 0.625, 0.875];
%     dv0 = zeros(1,4);
%     de0 = zeros(1,4);
%     e00 = zeros(1,4);
%     for i = 1:4
%         dv0(1,i) = (v0(1,i+1) - v0(1,i))/0.25;
%         de0(1,i) = (e0(1,i+1) - e0(1,i))/0.25;
%         e00(1,i) = (e0(1,i+1) + e0(1,i))/2;
%     end
%     k0 = de0;
%     y0 = dv0 - e00;
%     dv = zeros(N+1);
%     de = zeros(N+1);
%     for i = 1:N+1
%      x(i) = (i-1)/N;
%      [dv(i),de(i)] = Mdsolution(x(i));
%     end
%     k = de;
%     y = dv - e;
%   % Plot the bending and shear strain
%     figure(2)
%     ax11 = subplot(1,2,1);
%     plot(ax11,x00,k0,'r+',x,k,'b-');
%     title(ax11,'Bending strain')
%     legend(ax11,'FEM solution','Exact solution','Location','southeast');
%     ax21 = subplot(1,2,2);
%     plot(ax21,x00,y0,'r+',x,y,'b-');
%     title(ax21,'Shear strain')
%     legend(ax21,'FEM solution','Exact solution','Location','southeast');
end
end

function [v,e] = Msolution(x)
 v = 1/2000/0.00135*x*x;
 e = 1/1000/0.00135*x;
end
function [dv,de] = Mdsolution(x)
 dv = 1/1000/0.00135*x;
 de = 1/1000/0.00135;
end

function [v,e] = Fsolution(x)
 v = -1/6000/0.00135*x*x*x + 1/2000/0.00135*x*x + 1.2/400/0.18*x;
 e = -1/2000/0.00135*x*x + 1/1000/0.00135*x;
end
function [dv,de] = Fdsolution(x)
 dv = -1/2000/0.00135*x*x + 1/1000/0.00135*x + 1.2/400/0.18;
 de = -1/1000/0.00135*x + 1/1000/0.00135;
end

function y = Linear(x1,y1,x2,y2,x)
y = y1 + (y2 - y1)/(x2 - x1) * (x - x1);
end
function [weight, GaussPoints] = GetGauss(n)
if n == 3
 weight = [5/9 8/9 5/9];
 GaussPoints = [-sqrt(3/5) 0 sqrt(3/5)];
end
end