function [x,y] = CFD03_gridgene2D
    include_flags_CFD03;
    input_file_CFD03;
    
    [x,y]=initialization;
    
    figure(1);
    plot(x0,y0,'g-');
    
    figure(2);
    plot(x,y,'r-',x',y','r-');

    [x,y]=optimization(x,y);

    figure(3);
    plot(x,y,'r-',x',y','r-');
end