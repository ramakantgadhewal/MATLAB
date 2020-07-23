function CFD05_incompr2D_simple
% Now ONLY valid for STEADY case
% TODO: complete the sub-functions
%    1. update parameters
%    2. estimate and correct
%    3. add pressure field to postproc
  include_flags_CFD05_simple;
  input_file_CFD05_simple;
  
  % Test the algorithm index
  if ~ismember(ind_uiter,1)
    error('  Error in assignment of algorithm index (velocity iteration)!')
  end
  if ~ismember(ind_piter,1)
    error('  Error in assignment of algorithm index (delta pressure iteration)!')
  end
  
  % do the iteration
  nIter = 0;
  while (err_pres > err_pres_thre) && nIter < nIter_max
    update_uiter_paras;
    estimate_uv;
    update_piter_paras;
    correct_uvp;
    nIter = nIter + 1;
    if (mod(nIter,outputInt) == 0)
      fprintf('%s%i%s%e%s%e\n','The ',nIter,"th step gets p accuracy = ",err_pres);
    end
  end

  % Save the flow field
  u_post(1:Nx+1,1:Ny) = (3.*u(1:Nx,:) - u(2:Nx+1,:))./2;
  u_post(1:Nx+1,Ny+1) = (3.*u(2:Nx+1,Ny+1) - u(1:Nx,Ny+1))./2;
  v_post(1:Nx,1:Ny+1) = (3.*v(:,1:Ny) - v(:,2:Ny+1))./2;
  v_post(Nx+1,1:Ny+1) = (3.*v(Nx+1,2:Ny+1) - v(Nx+1,1:Ny))./2;
  [XP,YP]=meshgrid(y_p,x_p);
  [X,Y]=meshgrid(y_u,x_u);
  U = u_post'; V = v_post'; P = p;
  save(['./data/' outputName '.mat'],'X','Y','U','V','XP','YP','P'); 
  % plot_stream([outputName '.mat'],1);
end


function update_uiter_paras
  include_flags;
end

function estimate_uv
  include_flags;
  err_uv = err_uv_thre + 1;
  nIter = 0;
  while err_uv > err_uv_thre && nIter < nIter_max
    nIter = nIter + 1;
    u0 = u; v0 = v;

    err_uv = max(max(max(abs(u - u0))),max(max(abs(v - v0))));
  end
end


function update_piter_paras
  include_flags;
end

function correct_uvp
  include_flags;
  du = zeros(Nx+1,Ny);
  dv = zeros(Nx,Ny+1);
  dp = zeros(Nx,Ny);
  err_pres = err_pres_thre + 1;
  nIter = 0;
  while err_pres > err_pres_thre && nIter < nIter_max
    nIter = nIter + 1;
    p0 = p;

    err_pres = max(max(max(abs(p - p0))));
  end
  p = p + dp;
  u = u + du;
  v = v + dv;
end


function plot_stream(outputName,plotMode)

  % load and record the data
  load(outputName,'X','Y','U','V');
  xx = X(:,:); yy = Y(:,:);
  uu = U(:,:); vv = V(:,:);
  
  % determine the numbers
  N_over = 4;
  vert_streamline = 500;
  arrow_width = 0.5;
  
  % calculate the required info
  xmin = min(min(min(xx)));
  xmax = max(max(max(xx)));
  ymin = min(min(min(yy)));
  ymax = max(max(max(yy)));
  uMag = sqrt(uu.^2+vv.^2);
  uMag_max=max(max(max(uMag)));
  uMag_min=min(min(min(uMag)));
  
  if plotMode == 1
    Nsx = (size(xx,1)-1)/N_over;
    Nsy = (size(xx,2)-1)/N_over;
    Ns = Nsx*Nsy;
    Ns_gid = zeros(Ns);
    for kx = 1:Nsx
      ksx = (kx-1)*N_over + N_over/2;
      for ky = 1:Nsy
        ksy = (ky-1)*N_over + N_over/2;
        Ns_gid((kx-1)*Nsx+ky) = (ksx-1)*size(xx,2) + ksy;
      end
    end
  elseif plotMode == 2
    Nsx = (size(xx,1)-1)/N_over+1;
    Nsy = (size(xx,2)-1)/N_over+1;
    Ns = Nsx*Nsy;
    Ns_gid = zeros(Ns);
    for kx = 1:Nsx
      ksx = (kx-1)*N_over + 1;
      for ky = 1:Nsy
        ksy = (ky-1)*N_over + 1;
        Ns_gid((kx-1)*Nsx+ky) = (ksx-1)*size(xx,2) + ksy;
      end
    end
  else
    error('   Error in assignment of plotMode!')
  end
  
  % generate stream lines
  for k = 1:Ns
    startx = xx(Ns_gid(k));starty = yy(Ns_gid(k));
    sl_i = stream2(xx,yy,uu,vv,startx,starty,[0.1,vert_streamline]);
    sl_sum{k} = sl_i{1};
  end

  % generate color and color index
  N_color = 32;
  mcp = colormap(parula(N_color));
  [~,~,uMag_color] = histcounts(uMag(:),linspace(uMag_min,uMag_max,N_color));

  % draw the stream lines
  figure(1)
  hold on
  xlim([xmin,xmax])
  ylim([ymin,ymax])
  xy_ratio = get(gca, 'DataAspectRatio');
  xy_ratio = xy_ratio(1:2);
  xy_lim = axis;
  for k = 1:Ns
    plot(sl_sum{k}(:,1),sl_sum{k}(:,2),'color',mcp(uMag_color(Ns_gid(k)),:))
    if size(sl_sum{k},1) <= 1
      continue
    end
    arrow_direction = sl_sum{k}(end,:)-sl_sum{k}(end-1,:);
    plot_arrow(xy_lim,xy_ratio,sl_sum{k}(end,:),...
               mcp(uMag_color(Ns_gid(k)),:),arrow_width,arrow_direction)
  end
  hold off
  caxis([uMag_min,uMag_max])
  colorbar
end

function plot_arrow(xy_lim,xy_ratio,xy_arrow,arrow_color,arrow_width,arrow_direction)
  % initialize the arrow's shape
  arrow_0 = [0,0;-1,0.5;-1,-0.5];
  % normalize the direction
  a_dn = arrow_direction(:)./xy_ratio(:);
  a_dn = a_dn/sqrt(sum(a_dn.^2));
  d = (xy_lim(4)-xy_lim(3)+xy_lim(2)-xy_lim(1))/2;
  % zoom the arrow (according to the window)
  arrow_1 = arrow_0 * arrow_width * 0.03 * d;
  % rotate the arrow
  arrow_2 = arrow_1 * [a_dn(1),a_dn(2);-a_dn(2),a_dn(1)];
  % deform the arrow
  xy_ratio_n = xy_ratio/sqrt(sum(xy_ratio.^2));
  arrow_3 = arrow_2.*xy_ratio_n + xy_arrow;
  fill(arrow_3(:,1),arrow_3(:,2),arrow_color,'EdgeColor','none')
end