function [y, freqs] = fft_notch(x, Ts, lb, ub )
  % [pow_spec, freqs] = fft_notch(x, Ts, lb, ub )
  

  Fs = 1/Ts;

  Y = fft(x);
  L = length(x);
  P2 = abs(Y/L);
%   pow_spec = P2(1:L/2+1);
%   pow_spec(2:end-1) = 2*pow_spec(2:end-1);
  freqs = Fs*[0:L/2]/L;

  idx_lower = find(freqs <=lb, 1, 'last');
  idx_upper = find(freqs >= ub, 1, 'first');

%   keyboard

  Y(idx_lower:idx_upper) = 0;
  Y(L-idx_upper+2:L-idx_lower+2) = 0;

  y = ifft(Y);


end

