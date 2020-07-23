function plotbar; 
include_flags; 

% check if user defined the bar plot 
if strcmpi(plot_bar,'yes')==1;
    XX = zeros(1,nen+1);
    YY = zeros(1,nen+1);
    ZERO = zeros(1,nen+1);
    for i = 1:nel
        for j = 1:nen 
            XX(1,j) = x(IEN(j,i)); 
            YY(1,j) = y(IEN(j,i));  
        end
        XX(1,nen+1) = XX(1,1);
        YY(1,nen+1) = YY(1,1); 
        line(XX,YY);hold on;        
        line(XX,-YY);hold on; 
        plot(XX, ZERO, '+r');

        % check if user defined the plots of the global node numbering  
        if strcmpi(plot_nod,'yes')==1;    
            for j = 1:nen 
                text(XX(j),-1,sprintf('%0.5g',IEN(j,i)));
            end
        end
    end 
title('Bar Plot'); 
end 

% print some mesh parameters 
fprintf(1,'  Bar Params \n'); 
fprintf(1,'No. of Elements  %d \n',nel); 
fprintf(1,'No. of Nodes     %d \n',nnp); 
fprintf(1,'No. of Equations %d \n\n',neq);