function [y, freqs] = fft_highpass(x, Ts, ub )
  % [pow_spec, freqs] = power_spectrum(x, Ts )
  % For a vector of data x at sample rate Ts,
  % returns the power spectrum and associated frequencies.
  % Roughly this is the modulus of the first half of the FFT of the signal.

  Fs = 1/Ts;

  Y = fft(x);
  L = length(x);
  P2 = abs(Y/L);
  freqs = Fs*[0:L/2]/L;

  idx_lower = 2
  idx_upper = find(freqs >= ub, 1, 'first');

%   keyboard

  Y(idx_lower:idx_upper) = 0;
  Y(L-idx_upper+2:L-idx_lower+2) = 0;

  y = ifft(Y);


end

