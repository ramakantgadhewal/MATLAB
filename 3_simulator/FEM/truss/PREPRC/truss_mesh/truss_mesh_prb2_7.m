function  truss_mesh_prb2_7
% Include global flags
include_flags;

% Node (origin placed at node 1)
for i = 1:6
    x(2*i-1) = 3*(i-1);
    x(2*i)   = 3*(i-1);
    y(2*i-1) = 0;
    y(2*i)   = 4;
end
z = zeros(1,nnp);
% Connectivity array
IEN = [ 1   2      % element 1
        1   3      % element 2
        2   3      % element 3
        2   4      % element 4
        3   4      % element 5
        3   5      % element 6
        4   5      % element 7
        4   6      % element 8
        5   6      % element 9
        5   7      % element 10
        6   7      % element 11
        5   8      % element 12
        6   8      % element 13
        7   8      % element 14
        7   9      % element 15
        7   10     % element 16
        8   10     % element 17
        9   10     % element 18
        9   11     % element 19
        9   12     % element 20
        10  12     % element 21
        11  12]';  % element 22

% plot truss
plottruss;