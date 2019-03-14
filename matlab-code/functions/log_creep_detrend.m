
function uz_detrended = log_creep_detrend(uz, idx_state_s)



t = [0:1:length(uz)-1]'*AFM.Ts;


idx_scan = [];
for k=1:length(idx_state_s.scan)
  idx_ = idx_state_s.scan{k};
  idx_scan = [idx_scan; idx_(:)];
end
% idx_scan = [2:1:length(uu)-1]';

uz_scan = uz(idx_scan);
N = length(uz_scan);
tt = t(idx_scan);
alp = 0.1;
% regressor matrix
A = [ones(N,1), tt, log10(tt/alp)];

vvbg = A\uz_scan;

Vo = vvbg(1);
gam = vvbg(3);
b = vvbg(2);

z_fit = (Vo + b*t + gam*log10(t/alp));

uz_detrended = uz - z_fit;
end
