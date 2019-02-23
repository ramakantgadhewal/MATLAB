function [x,y,u,v,w,pf] = CFD05_squaredriven2D
    include_flags;
    input_file;
    
    % Iteration of (u,v)
    err_steady = 1;
    iterator = 0;
    tic
    while err_steady > err_steady_sup
        u0 = u;
        v0 = v;
        [u,v,w,pf] = iter(u,v,w,pf);
        % calculate the error
        UVmax = max(max(max(abs(u0),abs(v0))));
        err_steady = sqrt(sum(sum((u - u0).^2 + (v - v0).^2)))/max(UVmax,1);
        if UVmax_flag == 1 && UVmax > 1
            error('  u or v is too big!');
        end
        % print the err_steady every 10 steps
        iterator = iterator + 1;
        if (mod(iterator,10) == 0)
            fprintf('%s%i%s%f\n','The ',iterator,'th step gets accuracy =',err_steady);
        end
    end
    toc
    
    % Calculation of p
    if(cal_p_flag == 1)
        p = cal_p(u,v);
    end
    
    % Plot the flow field (and p distribution)
    [Y,X] = meshgrid(y,x);
    contour(X,Y,pf,50);
end