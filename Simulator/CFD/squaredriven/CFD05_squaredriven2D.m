function CFD05_squaredriven2D
  include_flags;
  input_file;
  
  % Test the algorithm index
  if ~ismember(ind_uv2omg,1:3)
    error('  Error in assignment of algorithm index (velocity to vorticity)!')
  end
  if ~ismember(ind_psifun,1)
    error('  Error in assignment of algorithm index (stream function iteration)!')
  end
  if ~ismember(ind_ps2omg,1:2)
    error('  Error in assignment of algorithm index (stream function to vorticity)!')
  end
  if ~ismember(ind_ps2uv,1:2)
    error('  Error in assignment of algorithm index (stream function to velocity)!')
  end
  
  % Calculate (u,v)
  [U,V] = cal_uv(U,V);

  % Calculate p
  pres = cal_pres(U,V);

  % Save the flow field
  [X,Y]=meshgrid(y,x); U = U'; V = V';
  save(['./data/' outputName '.mat'],'X','Y','U','V'); 
  % plot_stream([outputName '.mat'],1);
end