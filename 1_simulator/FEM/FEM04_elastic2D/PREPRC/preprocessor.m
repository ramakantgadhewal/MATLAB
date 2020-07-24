function [K,f,d] = preprocessor; 
include_flags;

% read input file 
% input_file_3T_2ele; 
% input_file_1ele;
% input_file_16ele; 
% input_file_64ele;
input_file_test4Q;

% IEN [identity of element nodes] -- Relates local node numbers to [global node numbers]
% ID  [identity array]            -- The relationship between global node numbers and global equation numbers ...
%                                       as well as nodal boundary condition information
% LM  [location matrix]           -- Gives the [global equation number] that corresponding to the ith DOF of element e

% generate ID and LM arrays
d = setup_ID_LM(d);
