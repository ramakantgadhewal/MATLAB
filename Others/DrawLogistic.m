function DrawLogistic()
lambda = 3.1;
n = 1e3;
t = 1:1:n;
X = zeros(n);
X(1) = 1-1/3.1-1e-8;
for i = 2:n
    X(i) = lambda*X(i-1)*(1-X(i-1));
end
plot(t,X,'r-');
end

