function [x_psd, freqs] = power_spectrum(x, Ts, win)
  % [pow_spec, freqs] = power_spectrum(x, Ts, win)
  % For a vector of data x at sample rate Ts,
  % returns the power spectrum and associated frequencies.
  % Roughly this is the modulus of the first half of the FFT of the signal.
  % If win is not included or empty, win is a rectangular window.
  Fs = 1/Ts;

  if ~exist('win', 'var') || isempty(win)
    win = rectwin(length(x));
  end
  
  % For a rectangular window, S = N.
  S = sum(win.^2);
  N = length(x);
  
  % Apply window. Identity operation if win = rectwin.
  x_windowed = x.*win;
  x_dft = fft(x_windowed, N);
  
  x_psd = x_dft.*conj(x_dft); %abs(x_dft).^2; % same as x_dft.*conj(x_dft).
  
  % take only one side. Do not scale DC by 2.
  x_psd = x_psd(1:floor(N/2)+1);
  x_psd = (2*Ts/S) * x_psd;
  x_psd(1) = x_psd(1)/2;
  
  freqs = Fs*[0:N/2]/N;
  
  
%   Y = fft(x);
%   L = length(x);
%   P2 = abs(Y/L);
%   pow_spec = P2(1:floor(L/2)+1);
%   pow_spec(2:end-1) = 2*pow_spec(2:end-1);
%   freqs = Fs*[0:L/2]/L;

end

