function [p] = Phy02_linear_plot(x,y)
    r = min(min(corrcoef(x,y)));
    p = polyfit(x,y,1);
    p1 = plot(x,y,'b*');
    hold on;
    p2 = plot(x,polyval(p,x),'k');
    hold off;
    title('\fontsize{14} lg ( I_e / T^2 ) ~ 1/T 线性拟合曲线');
    xlabel('1/T ( T单位: K )');
    ylabel('lg ( I_e / T^2) ( I单位: \muA, T单位: K )');
%     legend([p1,p2],'collected data',['(x_{1}, y_{1}) = (',num2str(x_maxima(1),'%5.2f'),', ',num2str(y_maxima(1),'%5.2f'),')'; ...
%         '(x_{2}, y_{2}) = (',num2str(x_maxima(2),'%5.2f'),', ',num2str(y_maxima(2),'%5.2f'),')']);
    legend([p1,p2],'Location','northeast',['\fontname{Cambria Math} collected data'],['\fontname{Cambria Math} lg ( I_e / T^2) = b+ k/T :', ...
        sprintf('\n'),'\fontname{Cambria Math} k = ',num2str(p(1)),sprintf('\n'),'\fontname{Cambria Math} r = ',num2str(r,'%0.6g')]);
    legend('boxoff');
end

