function [ x_filt, x_fft] = fft_sniper(x, Ts, freqs)
%IDEAL_LPF Summary of this function goes here
%   Detailed explanation goes here
 L = length(x);
 f_s = (1/Ts)*[0:(L/2)]'/L;
 
 
%  keyboard
x_fft = fft(x);

k_s = freqs*0;
k_s_conj = k_s;
for i=1:length(freqs)
    [~,idx] = min(abs(f_s - freqs(i)));
    k_s(i) = idx;
    k_s_conj(i) = find(x_fft == conj(x_fft(idx)), 1, 'last');
end



x_fft(k_s) = 0;
x_fft(k_s_conj) = 0;

x_filt = ifft(x_fft);

end

