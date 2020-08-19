function [I,I0,err,exp,k] = quad_Romberg(f,a,b,eps)
% 功能：利用Romberg算法计算任意给定函数在给定区间上给定精度的积分，并给出误差与迭代次数
% 输入变量：被积函数f，被积区间[a,b]，计算精度eps(默认为1e-7)
% 输出变量：数值积分值I，精确积分值I0，积分计算误差err，迭代次数 k.[程序运行时间"tic-toc"]
% 变量设置说明：
%   [1]符号变量：
%       函数f(x)的自变量     x
%   [2]数值变量：
%       求积步长             h(:)
%       减半后中点求和       H(:)
%       梯形公式减半加密值   T(:,1)
%       Romberg加速迭代值    T(:,n)[其中，n>1]
%       数值积分值           I
%       精确积分值           I0
%       积分计算误差         err
tic
if nargin == 3
    eps = 1e-7;
end
% 预分配内存
H = zeros(50,1);
T = zeros(50);
h = zeros(50,1);
% 迭代初值计算
T01 = (b-a) * (f(a) + f(b)) / 2;
h0 = b - a;
H0 = f(a + h0/2);
k = 1;
T(k,1) = (T01+h0*H0)/2;
T(k,k+1) = (4*T(k,1)-T01)/3;
error = abs(T01 - T(k,k+1));
% Romberg迭代加速
while (error >= eps)
    % 等距网格减半加密
    h(k) = (b - a)/(2^k);
    H(k) = 0;
    for l = 1:1:2^k
        H(k) = H(k) + f(a + (l-1/2)*h(k));
    end
    T(k+1,1) = (T(k,1)+h(k)*H(k))/2;
    % Richardson外推
    for j = 1:1:(k+1)
        T(k+1,j+1) = (4^j * T(k+1,j) - T(k,j))/(4^j - 1);
    end
    error = abs(T(k,k+1) - T(k+1,k+2));
    k = k + 1;
    if k > 50
        error('迭代可能不收敛！');
    end
end
if k > 1
    k = k - 1;
end
I   = vpa(T(k,k+1));
I0  = int(f,symvar(sym(f)),a,b);
err = vpa(abs(I-I0));
exp = log10(err);
toc
end