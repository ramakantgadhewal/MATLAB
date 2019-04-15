function SP01_draw_disper(b)
% 一维双原子链、力常数交替变化的晶格振动色散关系图
    n = 100;
    for j = 1:length(b)
        w0 = 1;
        q = -2*pi:(pi/n):2*pi;
        w1 = zeros(4*n+1);
        w2 = zeros(4*n+1);
        for i = 1:(4*n+1)
            w1(i) = w0^2*((2+b(j)) + sqrt((2+b(j))^2-4*(1+b(j))*sin(q(i)/2)*sin(q(i)/2)));
            w2(i) = w0^2*((2+b(j)) - sqrt((2+b(j))^2-4*(1+b(j))*sin(q(i)/2)*sin(q(i)/2)));
        end
        plot(q,w1,'r',q,w2,'k');
        hold on;
    end
    
end