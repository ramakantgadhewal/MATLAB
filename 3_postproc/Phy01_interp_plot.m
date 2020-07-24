%    function Phy01_interp_plot(x1, y1, x2, y2, num)
% % 》》主程序Ⅰ《《
% % 输入参数：
% %       x1, y1, x2, y2       插值点坐标           
% % 过程参数：
% %       x_interp, y_interp   插值曲线坐标
% 
% % 判断数组长度是否符合要求
%     m1 = length(x1);
%     n1 = length(y1);
%     m2 = length(x2);
%     n2 = length(y2);
%     if (m1 ~= n1) || (m2 ~= n2)
%         disp '  >> Array length mismatch error! <<'
%         return
%     end
% % 作拟合曲线
%     [x1_interp, y1_interp] = inter_spline(x1, y1, num);
%     [x2_interp, y2_interp] = inter_spline(x2, y2, num);
%     p1 = plot(x1,y1,'b*');
%     hold on;
%     p2 = plot(x2,y2,'bo');
%     p3 = plot(x1_interp, y1_interp,'k-.');
%     p4 = plot(x2_interp, y2_interp,'k');
%     hold off;
% % 给出图例
%     title('\fontsize{14}R-L-C串联电路幅频特性曲线');
%     xlabel('频率f (Hz)');
%     ylabel('电流I (mA)');
%     legend('collected data (C=1\muF,R=20\Omega)','collected data (C=1\muF,R=50\Omega)','I_{1}( f ) [f_{0} = 482.56Hz]','I_{2}( f ) [f_{0} = 482.42Hz]');
%     legend('boxoff');
%   end
    
  function [x_maxima, y_maxima] = Phy01_interp_plot(x, y, num)
% 》》主程序Ⅱ《《
% 输入参数：
%       x, y                 插值点坐标           
% 过程参数：
%       x_interp, y_interp   插值曲线坐标
% 输出参数：
%       x_maxima, y_maxima   分段极值点坐标

% 判断数组长度是否符合要求
    m = length(x);
    n = length(y);
    if (m ~= n)
        disp '  >> Array length mismatch error! <<'
        return
    end
% 作拟合曲线
    [x_interp, y_interp] = Phy01_inter_spline(x,y,num);
%     [x_maxima, y_maxima] = find_maxima(x_interp, y_interp, num, 2);
    p1 = plot(x,y,'b*');
    hold on;
%     p2 = plot(x_maxima, y_maxima,'ro');
    p2 = plot(x_interp, (1+cos(x_interp/90*pi))/2,'g');
    p3 = plot(x_interp, y_interp,'r');
    hold off;
% 给出图例
    title('\fontsize{14}相对透射率p~差角\theta关系曲线');
    xlabel('差角\theta (\circ)');
    ylabel('相对透射率p');
%     legend([p1,p2],'collected data',['(x_{1}, y_{1}) = (',num2str(x_maxima(1),'%5.2f'),', ',num2str(y_maxima(1),'%5.2f'),')'; ...
%         '(x_{2}, y_{2}) = (',num2str(x_maxima(2),'%5.2f'),', ',num2str(y_maxima(2),'%5.2f'),')']);
    legend([p1,p2,p3],'collected data','p = cos^2\theta','interp curve');
    legend('boxoff');
end

function [x_interp, y_interp] = Phy01_inter_linear(x,y,num)
% 分段线性插值
%     range_of_x = max(x) - min(x);
%     x_min = min(x) - range_of_x * 0.1;
%     x_max = max(x) + range_of_x * 0.1;
    x_min = min(x);
    x_max = max(x);
    x_interp = linspace(x_min, x_max, num);
    y_interp = interp1(x, y, x_interp, 'linear');
end

function [x_interp, y_interp] = Phy01_inter_spline(x,y,num)
% 三次样条插值
%     range_of_x = max(x) - min(x);
%     x_min = min(x) - range_of_x * 0.1;
%     x_max = max(x) + range_of_x * 0.1;
    x_min = min(x);
    x_max = max(x);
    x_interp = linspace(x_min, x_max, num);
    y_interp = interp1(x, y, x_interp, 'spline');
end

function [x_interp, y_interp] = Phy01_inter_pchip(x,y,num)
% 共形插值
    range_of_x = max(x) - min(x);
    x_min = min(x) - range_of_x * 0.1;
    x_max = max(x) + range_of_x * 0.1;
    x_interp = linspace(x_min, x_max, num);
    y_interp = interp1(x, y, x_interp, 'pchip');
end

function [x_maxima, y_maxima] = Phy01_find_maxima(x,y,num,n_of_maxima)
% 寻找分段极值点
    x_maxima = zeros(n_of_maxima,1);
    y_maxima = zeros(n_of_maxima,1);
    len = num / n_of_maxima;
    i = 1;
    while(i <= n_of_maxima)
        n1 = int16(len*(i-1)+1);
        n2 = int16(len*i);
        [y_maxima(i),t(i)] = max(y(n1:n2)); 
        x_maxima(i) = x(n1 + t(i) - 1);
        i = i + 1;
    end
end
