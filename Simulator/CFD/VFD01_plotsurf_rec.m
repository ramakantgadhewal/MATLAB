function [dtheta] = VFD01_plotsurf_rec(a,b,n,N)
% -------------输入变量-------------
%   a   接触角
%   b   无量纲量 b:=rho*g*L^2/gamma
%    其中，
%     密度         rho  1000(水) 
%     重力加速度    g    9.8
%     特征长度      L    0.001
%     表面张力系数 gamma 0.0728(水-空气)
%   n   Plot footstep length
%   N   Verification footstep length
% -------------输出变量-------------
%   dtheta  接触角误差
% ---------------end---------------
global const alpha sigma
alpha = a;
sigma = b;
% calculate the initial value of z and the integral constant
z0    = sqrt(2*(1 - sin(alpha))/sigma);
sq0   = sqrt(1-sigma*z0*z0/4);
const = - log(1+sq0) + 2*sq0 + log(z0);
% plot the curve
x = zeros(1,n);
y = zeros(1,n);
d = z0/n;
for i = 1:n
    y(i) = z0 - d*(i-1);
    x(i) = leng(y(i));
end
plot(x,y);
% verificate the curve
y1 = z0;
y2 = z0 - 1/N;
x1 = leng(y1);
x2 = leng(y2);
dtheta = atan((x2-x1)/(y1-y2)) - alpha;
end

function [x] = leng(z)
global const sigma
% define the curve
sq = sqrt(1-sigma*z*z/4);
x  = (log(1+sq) - 2*sq - log(z) + const)/sqrt(sigma);
end

