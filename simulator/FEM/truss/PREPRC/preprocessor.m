function  [K,f,d0] = preprocessor;
%% Include global flags
include_flags;

%% input file to include all variables 
%  input_file_example2_2;
%  input_file_example2_8;
%  input_file_problem2_4;
 input_file_problem2_7;

%% compute the geometry of the elements
 leng = zeros(1,nel);          % Elements length
 phi  = zeros(nsd,nel);        % the COSINE of the angle
 A(1,:) = x;
 A(2,:) = y;
 if nsd == 3
     A(3,:) = z;
 end
 for e = 1:nel
     for i = 1:nsd
         leng(e) = leng(e) + (A(i,IEN(1,e))-A(i,IEN(2,e)))^2;
     end
     leng(e) = sqrt(leng(e));
     for i = 1:nsd
         phi(i,e) = (A(i,IEN(2,e))-A(i,IEN(1,e))) / leng(e);
     end
 end

%% compute the other normalizing coefficients
 l_ex  = floor(log10(max(leng)));       % Normalizing coefficient for the length of elements
 d_ex  = f_ex + l_ex - E_ex - CArea_ex; % Normalizing coefficient for the node displacements

%% Generate LM array 
for e = 1:nel
    for j = 1:nen
        for m = 1:ndof
            ind = (j-1)*ndof + m;
            LM(ind,e) = ndof*IEN(j,e) - ndof + m;
        end
    end
end
