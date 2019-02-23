function  truss_mesh_ex2_8;
include_flags;


% Node:  1    2    3   (origin placed at node 2) 
%--------------------
x   =  [0.0  1.0  2.0  ];     % X coordinate  
y   =  [0.0  0.0  0.0  ];     % Y coordinate

% connectivity array
IEN =  [1    2          % element 1
        2    3]';       % element 2

% plot truss
plottruss;


