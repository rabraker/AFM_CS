function [ x_filt, x_fft] = ideal_lpf(x, f_bw, Ts)
%IDEAL_LPF Summary of this function goes here
%   Detailed explanation goes here
 L = length(x);
 f_s = (1/Ts)*[0:(L/2)]'/L;
 
 
%  keyboard
x_fft = fft(x);

k = find(f_s > f_bw, 1,'first');
kend = find(x_fft == conj(x_fft(k)), 1, 'last');

x_fft(k:kend) = 0;

x_filt = ifft(x_fft);

end

