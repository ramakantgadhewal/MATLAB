function [x,y] = main
    include_flags;
    input_file;
    
    [x,y]=initialization;
    
    figure(1);
    plot(x0,y0,'g-');
    
    figure(2);
    plot(x,y,'r-',x',y','r-');

    [x,y]=optimization(x,y);

    figure(3);
    plot(x,y,'r-',x',y','r-');
end