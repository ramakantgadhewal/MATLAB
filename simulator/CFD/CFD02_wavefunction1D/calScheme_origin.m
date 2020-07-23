function [K,f,b,k0_plot] = calScheme_origin(l,r,e_res)
%% Unoptimized scheme calculation
% calculate b
Num = l + r + 1;
K = zeros(Num);
f = zeros(Num,1);
f(2,1) = 1;
K(1,:) = 1;
for i = 2:Num
    for j = 1:Num
        K(i,j) = (- l + j - 1)^(i-1);
    end
end
b = K\f;
% Calculate dispersion and dissipation error
N = 100;
PI = floor(pi*N);
k_plot = 0:1/N:PI/N;
Rek_plot = zeros(PI+1,1);
Imk_plot = zeros(PI+1,1);
res_flag = 1;
k0_plot = 0;
for i = 1:PI+1
    ki = (i-1)/N;
    for j = 1:Num
        Rek_plot(i,1) = Rek_plot(i,1) + b(j) * sin((- l + j - 1)*ki);
        Imk_plot(i,1) = Imk_plot(i,1) - b(j) * cos((- l + j - 1)*ki);
    end
    if i > 1 
        if (abs(Rek_plot(i,1)/k_plot(i) - 1) < e_res && res_flag == 1)
            k0_plot = k0_plot + 1/N;
        else
            res_flag = 0;
        end
    end
end
figure(100);
plot(k_plot,k_plot,'k-',k_plot,Rek_plot,'r-',k_plot,Imk_plot,'b-');
title('Original UplDownrF scheme');
legend('Re-accurate','Re-dispersion','Im-dissipation');
xlabel('wavenumber k');
ylabel('Rek^{\prime} or Imk^{\prime}');
hold off;
end
