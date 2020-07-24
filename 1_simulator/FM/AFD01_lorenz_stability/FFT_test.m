% sampling frequency
fs = 8;
% footstep in iteration
h = 1/fs;
% time series
t = 0:1/fs:100-1/fs;
% sample number
n = length(t);
% signal initialization and generation
x0 = zeros(1,n);
for i = 1:10
  x0 = x0 + cos(2*pi*i/10*t);
end
% fftshift and power spectrum
Y_x0 = fft(x0);
Pshift_x0 = abs(fftshift(Y_x0)).^2/n;
fshift_x0 = (-n/2:n/2-1)*fs/n;
plot(fshift_x0,Pshift_x0);