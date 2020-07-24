function [K,u,r] = FEM01_directStiff1D(m,n,p,k,F,L0)
%-------变量说明-------%
%   m   单元数         %
%   n   系统节点数     %
%   p   位移约束个数   %
%   k   单元刚度向量   %
%   F   外力变量       %
%   L0  标号转移核(2*m)%
%---------end---------%
syms u(t) r(t); % 定义符号变量：位移、约束力

K0 = [1 -1
      -1 1];              % 单元刚度矩阵核 赋值
L  = zeros(2,n);          % 标号转移阵 初始化
K  = zeros(n);            % 系统刚度阵 初始化

for j = 1:m               % 组装系统刚度阵
    for i = 1:2
        L(i, L0(i,j)) = 1;
    end
    K1 = k(j) * K0;
    K = L'* K1 * L + K;
    L = zeros(2,n);       % 标号转移阵清零
end

p = p+1;                  % 计算系统刚度阵、节点位移、约束力并输出
u = sym(K(p:n,p:n))\sym(F(p:n));
p = p-1;
r = sym(K(1:p,(p+1):n))*sym(u);
end

