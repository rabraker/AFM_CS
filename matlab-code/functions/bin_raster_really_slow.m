
function [ pixmat, pixelifsampled, m_s] = bin_raster_really_slow(xyu_dat, nperiods, samps_per_period, volt2pix, Hz)
if size(xyu_dat, 1) ~= nperiods*samps_per_period
   s = sprintf('Expected length(xyu_dat) == nperiods*samps_per_period')
   s = sprintf('%s \n but have length(xyu_dat)=%d, nperiods*samps_per_period = %d',...
       s, length(xyu_dat), nperiods*samps_per_period)
    error(s)
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

% Make the image be at (0,0, 0).
xdat_trace = (xdat_trace - min(xdat_trace))*volt2pix;
ydat_trace = (ydat_trace - min(ydat_trace))*volt2pix;

% make 
% udat_trace = (udat_trace - min(udat_trace));
% udat_trace = udat_trace/max(udat_trace);


pixelifsampled = zeros(xpix, ypix);
pixmat = zeros(xpix,ypix)*mean(udat_trace);
m_s = zeros(ypix,1);
for j_row = 0:ypix-1
%     ind_y = find(ydat_trace >= j_row & ydat_trace < j_row+1);
    ind_y = j_row*(samps_per_period/2)+1:(j_row+1)*(samps_per_period/2);
    x_dat_j = xdat_trace(ind_y)';
    U_dat_j = udat_trace(ind_y)';
    
    % detrend each row
%     keyboard
if exist('Hz', 'var')
  t = (0:length(U_dat_j)-1)'*Hz.Ts;
  U_dat_j = lsim(Hz, detrend(U_dat_j'), t);
  %filtfilt(Hz.num{1}, Hz.den{1}, detrend(U_dat_j'));
%   figure(100); hold on;
%   plot(detrend(U_dat_j'))
%   plot(U_, '--');
%   keyboard
else
    [U_dat_j, mb] = detrend(U_dat_j');
%     U_dat_j = fft_notch(U_dat_j(:), 40e-6, 210, 216);
%     U_dat_j = fft_notch(U_dat_j(:), 40e-6, 28, 31);
end
    U_dat_j = U_dat_j';
%     m_s(j_row+1) = mb(1);
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

