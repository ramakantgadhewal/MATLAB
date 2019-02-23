function  truss_mesh_prb2_4
include_flags;

% Node:   1    2    3    4 (origin placed at node 4) 
%-------------------------
x   =  [-1.0  0.0  1.0  0.0];     % X coordinate  
y   =  [ 1.0  1.0  1.0  0.0];     % Y coordinate
% z   =  [0.0  0.0  0.0  ...];     % Z coordinate

% connectivity array
% Element: 1    2    3
%------------------------
IEN =    [ 1    2    3    
           4    4    4];     

% plot truss
plottruss;