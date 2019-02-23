function [x,y] = initialization(m,n)
include_flags;

x = zeros(2*m+num+1,n+1);
y = zeros(2*m+num+1,n+1);
xmax = x0(1,1)+len*fac;
ymax = len*fac;

lx = xmax/(exp(mag)-1);
ly = ymax/(exp(mag)-1);
for i = 1:(m+1)
    for j = 1:(n+1)
        x(i,j) = x0(1,1) + (exp(mag*(m+1-i)/m)-1)*lx;
        y(i,j) = (exp(mag*(j-1)/n)-1)*ly;
        iacross = 2*m + num + 2-i;
        x(iacross,j) = x(i,j);
        y(iacross,j) = -y(i,j);
    end
end

for i = (m+2):(m+num)
    % determines the angle distribution of far-field boundary points on the left;
    thetai = cal_theta(i);
    x(i,n+1) = x0(1,1)+0.5*ymax*cos((thetai+0.5)*pi);
    y(i,n+1) = ymax*sin((thetai+0.5)*pi);
    x(i,1) = x0(i-m,1);
    y(i,1) = y0(i-m,1);
    lx = (x(i,n+1)-x(i,1))/(exp(mag)-1);
    ly = (y(i,n+1)-y(i,1))/(exp(mag)-1);
    for j = 2:n
        x(i,j) = x(i,1) + (exp(mag*(j-1)/n)-1)*lx;
        y(i,j) = y(i,1) + (exp(mag*(j-1)/n)-1)*ly;
    end
end
end

function [thetai] = cal_theta(i)
% theta() determines the angle distribution of far-field boundary points on the left;
% - where it is required that theta(m+1) = 0 and theta(m+num) = 1.
% - It is also required that the points are concentrated around the leading
% - edge and trailing edge of the airfoil.
    include_flags;
    angle_i = 2*pi*(i-m-1)/(num-1);
    if angle_i < pi
        thetai = 1/4-1/4*cos(angle_i);
    else
        thetai = 3/4-1/4*cos(angle_i-pi);
    end
end
