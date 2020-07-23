function CFD04_compressible1D
include_flags_CFD04;
input_file_CFD04;

%% Call iterators
if (scheme_flag == 1)
    call_Rusanov;
elseif (scheme_flag == 2)
    call_Jameson;
elseif (scheme_flag == 3)
    call_1orderFVS;
elseif (scheme_flag == 4)
    call_2orderFVS;
elseif (scheme_flag == 5)
    call_2orderFVS_Runge;
elseif (scheme_flag == 6)
    call_1orderDRoe;
elseif (scheme_flag == 7)
    call_1orderVRoe;
elseif (scheme_flag == 8)
    call_2orderVTVD;    
end

%% plot the shockwave
plot_shockwave;

end