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

function [xtraj0, ytraj0, xtc, ytc] = mpt_connect_rad(xtraj0, ytraj0, xtc, ytc, opts)
%   rad_min, volts_per_sample, ax)
 
  if nargin < 6
    ax = [];
  end
  
  if length(xtraj0) <2
    return
  end
  x1e = xtc(end);
  y1e = ytc(end);
  [xr, yr] = get_xr_yr(xtraj0, ytraj0);
  
  rdiff = sqrt( (xr - x1e).^2 + (yr - y1e).^2);
  
  [rmin, rmin_idx] = min(rdiff);
  
  if opts.ensure_forward
    optional_cond = (x1e < xtraj0{rmin_idx}(1)); % x2s
  else
    optional_cond = true;
  end
  
  if (rmin < opts.rad_min) && optional_cond
    % Extract the mu-path to which we are connecting, delete it from the
    % (xtraj0, ytraj0) cells, and append it and a connection trajectory 
    % to xtc, ytc, and recursively call back into mpt_connect_rad().
    
    % Extract:
    xt2 = xtraj0{rmin_idx};
    yt2 = ytraj0{rmin_idx};
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
    
    if ~isempty(opts.ax)
      plot(opts.ax, xcon, ycon, 'r')
      plot(opts.ax, xt2, yt2, '--g')
      plot(opts.ax, xtc, ytc, '--k')
      xlim(opts.ax, [0,1])
      ylim(opts.ax, [0,1])
      drawnow()
    end
    xtc = [xtc(:); xcon(:); xt2];
    ytc = [ytc(:); ycon(:); yt2];
     
     % Delete
    xtraj0(rmin_idx) = [];
    ytraj0(rmin_idx) = [];
    
    [xtraj0, ytraj0, xtc, ytc] = mpt_connect_rad(xtraj0, ytraj0, xtc, ytc, opts);
  end
  
  % Otherwise, do nothing.

end

function [xr, yr] = get_xr_yr(xtraj, ytraj)
  % Returns the first element of each array contained in the cell arrays xtraj
  % ytraj.
  xr = zeros(length(xtraj), 1);
  yr = zeros(length(ytraj), 1);
  
  for k =1:length(xtraj)
    xr(k) = xtraj{k}(1);
    yr(k) = ytraj{k}(1);
  end

end