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