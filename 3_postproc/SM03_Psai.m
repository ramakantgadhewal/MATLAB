function [d,y] = SM03_Psai(I, b)
b = b/180*pi;
d = asin(2*sqrt(I)/(sin(2*b)*(1+I)));
y = 0.5*atan(tan(2*b)*cos(d))*180/pi;
end