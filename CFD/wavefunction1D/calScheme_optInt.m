function [b_opt,k0_opt] = calScheme_optInt(l,r,prL,prR,nu,e_res)
%% Optimization using Integration
[K,f,~,~] = calScheme_origin(l,r,e_res);
% treat the last several variables as parameters
Num = l + r + 1;
Num_m = Num - prL - prR;
Num_d = prL + prR;
Kmd(1:Num_m, 1:prL) = K(1:Num_m, 1:prL);
Kmm = K(1:Num_m, prL+1:Num-prR);
Kmd(1:Num_m, prL+1:Num_d) = K(1:Num_m, Num-prR+1:Num);
fm = f(1:Num_m);
bd = sym('optArg_',[Num_d,1]);
bm = Kmm\(fm-Kmd*bd);

% optimize the parameters
syms k Rek_opt Imk_opt;
Rek_opt = 0;
Imk_opt = 0;
for j = 1:Num_m
    Rek_opt = Rek_opt + bm(j) * sin((- l + j - 1)*k);
    Imk_opt = Imk_opt - bm(j) * cos((- l + j - 1)*k);
end
E = exp(-nu*pi)*int(exp(nu*(pi-k))*((Rek_opt-k)^2),k,[0,pi]);
Ex = symvar(E);
dE = sym(zeros(length(Ex),1));
for i = 1:length(Ex)
    dE(i) = diff(E,Ex(i),1);
end
if length(Ex) > 1
    bd_opt = struct2cell(solve(dE == 0, Ex));
else
    bd_opt = solve(dE == 0, Ex);
end
bd_temp = double(subs(bd,bd,bd_opt));
b_opt(1:prL,1) = bd_temp(1:prL);   % [IMPORTANT] convert the symbolic variable to numeric one
b_opt(prL+1:Num-prR,1) = double(subs(bm,bd,bd_opt));
b_opt(Num-prR+1:Num,1) = bd_temp(prL+1:Num_d);

% Calculate dispersion and dissipation error
N = 100;
PI = floor(pi*N);
k_opt = 0:1/N:PI/N;
Rek_opt = zeros(PI+1,1);
Imk_opt = zeros(PI+1,1);
res_flag = 1;
k0_opt = 0;
for i = 1:PI+1
    ki = (i-1)/N;
    for j = 1:Num
        Rek_opt(i,1) = Rek_opt(i,1) + b_opt(j) * sin((- l + j - 1)*ki);
        Imk_opt(i,1) = Imk_opt(i,1) - b_opt(j) * cos((- l + j - 1)*ki);
    end
    if i > 1 
        if (abs(Rek_opt(i,1)/k_opt(i) - 1) < e_res && res_flag == 1)
            k0_opt = k0_opt + 1/N;
        else
            res_flag = 0;
        end
    end
end
figure(200);
plot(k_opt,k_opt,'k-',k_opt,Rek_opt,'r-',k_opt,Imk_opt,'b-');
title('Optimized UplDownrF scheme (using Integration method)');
legend('Re-accurate','Re-dispersion','Im-dissipation');
xlabel('wavenumber k');
ylabel('Rek^{\prime} or Imk^{\prime}');
hold off;
end