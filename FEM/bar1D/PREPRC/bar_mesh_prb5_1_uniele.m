function bar_mesh_prb5_1_uniele
include_flags;

x = 0:20/(nnp-1):20;              % x coordinate
y = 1/sqrt(pi)*ones(1,nnp); % y is used only for the bar plot

% connectivity array
IEN = zeros(nen,nel);
for i = 1:nel
    for j = 1:nen
        IEN(j,i) = (i-1)*(nen-1)+j;
    end
end
plotbar;