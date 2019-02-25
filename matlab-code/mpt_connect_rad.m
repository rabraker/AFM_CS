% 
% xtraj0, ytraj0:  cell-arrays of xy-trajectories
% 
% xtrajc, ytrajc: arrays (for now, eventually, need to be somethings like a
% cell array also, so we can embed the meta data for when not to sample.)
% 
% On entry, (xtrajc, ytrajc) shall contain at least 1 mu-path trajectory. THis
% mu path should not be contained in (xtraj0, ytraj0). mpt_connect_rad will
% search over (xtraj0, ytraj0) for the mu-path whose start is closest to
% (xtrajc(end), ytrajc(end)). If one is found, that trajectory will be appended
% to xtrajc, ytrajc and removed from (xtraj0, ytraj0). 
% 
% 
% mpt_connect_rad() will recursively call itself until no connection is made.

function [mptc_s, mptc] = mpt_connect_rad(mptc_s, mptc, opts)
%   rad_min, volts_per_sample, ax)

  
  if length(mptc_s) <2
    return
  end
  x1e = mptc.xt(end);
  y1e = mptc.yt(end);
  [xr, yr] = get_xr_yr(mptc_s);
  
  rdiff = sqrt( (xr - x1e).^2 + (yr - y1e).^2);
  
  [rmin, rmin_idx] = min(rdiff);
  
  if opts.ensure_forward
    optional_cond = (x1e < mptc_s{rmin_idx}.xt(1)); % x2s
  else
    optional_cond = true;
  end
  
  if (rmin < opts.rad_min) && optional_cond
    % Extract the mu-path to which we are connecting, delete it from the
    % (xtraj0, ytraj0) cells, and append it and a connection trajectory 
    % to xtc, ytc, and recursively call back into mpt_connect_rad().
    
    % Extract:
    xt2 = mptc_s{rmin_idx}.xt;
    yt2 = mptc_s{rmin_idx}.yt;
    x2s = xt2(1);
    y2s = yt2(1);
    

    % Build the connecting trajectory. 
    if abs(x1e - x2s) < abs(x1e - x2s) 
      xcon = (x1e:sign(x2s-x1e)*opts.volts_per_sample:x2s)';
      ycon = linspace(y1e, y2s, length(xcon))';
    else
      ycon = (y1e:sign(y2s-y1e)*opts.volts_per_sample:y2s)';
      xcon = linspace(x1e, x2s, length(ycon))';
    end
    % The meta data is negative of scan idx.
    met_idx_con = ones(length(xcon), 1) * (-1*mptc.met_idx(1));
    met_idx_ext = ones(length(xt2), 1) * mptc.met_idx(1);
    
    if ~isempty(opts.ax)
      plot(opts.ax, xcon, ycon, 'r')
      plot(opts.ax, xt2, yt2, '--g')
      plot(opts.ax, mptc.xt, mptc.yt, '--k')
      xlim(opts.ax, [0,1])
      ylim(opts.ax, [0,1])
      drawnow()
    end
    mptc.xt = [mptc.xt(:); xcon(:); xt2];
    mptc.yt = [mptc.yt(:); ycon(:); yt2];
    mptc.met_idx = [mptc.met_idx(:); met_idx_con(:); met_idx_ext(:)];
    
     % Delete
     mptc_s(rmin_idx) = [];
    
    [mptc_s, mptc] = mpt_connect_rad(mptc_s, mptc, opts);
  end
  
  % Otherwise, do nothing.

end

function [xr, yr] = get_xr_yr(mptc_s)
  % Returns the first element of each array contained in the cell arrays xtraj
  % ytraj.
  xr = zeros(length(mptc_s), 1);
  yr = zeros(length(mptc_s), 1);
  
  for k =1:length(mptc_s)
    xr(k) = mptc_s{k}.xt(1);
    yr(k) = mptc_s{k}.yt(1);
  end

end