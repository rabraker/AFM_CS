
function [ pixmat, pixelifsampled, m_s] = bin_raster_really_slow(xyu_dat,...
    nperiods, samps_per_period, volt2pix, line_detrender)
  if size(xyu_dat, 1) ~= nperiods*samps_per_period
    s = sprintf('Expected length(xyu_dat) == nperiods*samps_per_period')
    s = sprintf('%s \n but have length(xyu_dat)=%d, nperiods*samps_per_period = %d',...
      s, length(xyu_dat), nperiods*samps_per_period)
    error(s)
  end
  if ~exist('line_detrender', 'var') || isempty(line_detrender)
    line_detrender = @(x) line_detrender_local(x);
  end
  xpix = nperiods;
  ypix = nperiods; % TODO: make this work with rectangular image.
  
  
  % detrend the drift from the height data.
  % xyu_dat(:,3) = detrend(xyu_dat(:,3));
  % Get the indeces corresponding to trace data only.
  trace_inds = get_trace_indeces(nperiods, samps_per_period);
  
  
  xdat_trace = xyu_dat(trace_inds, 1);
  ydat_trace = xyu_dat(trace_inds, 2);
  udat_trace = xyu_dat(trace_inds, 3);
  
  % Make the image be at (0,0, --).
  xdat_trace = (xdat_trace - min(xdat_trace))*volt2pix;
  ydat_trace = (ydat_trace - min(ydat_trace))*volt2pix;
  % make
  
  pixelifsampled = zeros(xpix, ypix);
  pixmat = zeros(xpix,ypix)*mean(udat_trace);
  m_s = zeros(ypix,1);
  for j_row = 0:ypix-1
    ind_y = j_row*(samps_per_period/2)+1:(j_row+1)*(samps_per_period/2);
    x_dat_j = xdat_trace(ind_y)';
    U_dat_j_init = udat_trace(ind_y)';
    
    [U_dat_j, mb] = line_detrender(U_dat_j_init);
    
    for i_col = 0:xpix-1
      ind_x = find(x_dat_j >= i_col & x_dat_j < i_col+1);
      
      u_mean_ij = mean(U_dat_j(ind_x));
      if ~isnan(u_mean_ij)
        pixmat(j_row+1, i_col+1) = u_mean_ij;
        pixelifsampled(j_row+1, i_col+1) = 1;
      else
        %             keyboard
      end
      
    end
    
  end

end

function [x, mb] = line_detrender_local(x)
  mb = [];
  return
end
