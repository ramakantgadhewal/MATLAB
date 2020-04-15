function [a,b] = airfoildata
    include_flags;
    a = zeros(num,1);
    b = zeros(num,1);
    a(1,1) = 1;
    b(1,1) = 0;
    a(num/2+1,1) = 0;
    b(num/2+1,1) = 0;
    for i = 2:num/2
        % determine the value of x along airfoil
        a(i,1) = cal_xi(i);
        % calculate the thickness on x
        b(i,1) = upper(a(i,1));
        a(num-i+2,1) = a(i,1);
        b(num-i+2,1) = - b(i,1);
    end 
end

function [xi] = cal_xi(i)
% xi() determines the distrqibution of x along airfoil;
% - where it is required that xi(1) = 1; xi(num/2+1) = 0.
% - It is also required that the points are concentrated around the leading
% - edge and trailing edge of the airfoil.
    include_flags;
    angle_i = pi * 2*(i-1)/num;
    xi = 1/2+1/2*cos(angle_i);
end

function [y] = upper(x)
% calculate y-coordinate of upper-surface points on the airfoil
    y = 0.6*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4);
end