function [L,lam,s] = CFD01_TYPE(A)
% ����������
% -------- ������� -------- %
%       ���ž���     A 
% -------- ������� -------- %
%       ����         s
%       ����ֵ       lam
%       ����������    L
% ---------- end ---------- %
s = symvar(A);

syms lam L;       % ����ֵ������������
[L,lam] = eig(A.');
L = L.';
end
