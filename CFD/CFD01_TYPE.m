function [L,lam,s] = CFD01_TYPE(A)
% 功能描述：
% -------- 输入变量 -------- %
%       符号矩阵     A 
% -------- 输出变量 -------- %
%       符号         s
%       特征值       lam
%       左特征向量    L
% ---------- end ---------- %
s = symvar(A);

syms lam L;       % 特征值，左特征矩阵
[L,lam] = eig(A.');
L = L.';
end

