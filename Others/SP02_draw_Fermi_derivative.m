function SP02_draw_Fermi_derivative
    Nj = 100;
    Ni = 100;
    E0 = 1;
    T = zeros(Ni,1);
    E = zeros(Nj,1);
    f = zeros(Ni,Nj);
    for i = 1:Ni
        T(i,1) = i*1e-4;
    end
    for j = 1:Nj
        E(j,1) = E0 + (j-Nj/2)*1e-5;
    end
    for i = 1:Ni
        for j = 1:Nj
            f(i,j) = T(i,1)*(exp((E(j,1)-E0)/T(i,1))+1)*(exp(-(E(j,1)-E0)/T(i,1))+1);
        end
    end
    [X,Y] = meshgrid(E,T);
    mesh(X,Y,f);
end