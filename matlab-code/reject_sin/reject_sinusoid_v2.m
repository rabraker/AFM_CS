clear, clc
addpath('functions')
Ts = 40e-6;
Ki_s = [0.005, 0.01, 0.05];
LPF = zpk([], [0.7, 0.7], 1, Ts);
LPF = LPF/dcgain(LPF);
models = load(fullfile(PATHS.sysid, 'x-axis_sines_infoFourierCoef_10-21-2018-03.mat'));
G = balreal(-models.modelFit.G_zdir);
sys_og = balreal(-models.modelFit.G_zdir_red);

% sys = ss(zpk([], 0.95, 1-0.95, Ts));
% G = sys;
Ki_s = [0.005, 0.01, 0.05];


root = '/media/labserver/afm-cs/z-scope';


% fnames = {'data-out_KI_p005.csv',...
%   'data-out_KI_p01.csv',...
%   'data-out_KI_p1.csv',...
%   'data-out_KI_p00.csv',...
%   'data-out_KI_p05.csv'};
f2 = 'data-out_KI_p01.csv';
fname = 'data-out_KI_p1.csv';
dat_k = csvread(fullfile(root, f2));
dat_k = csvread(fullfile(root, fname));

ze_k = detrend(dat_k(:,1));
uz_k = detrend(dat_k(:,2));
t_k = (0:length(ze_k)-1)'*Ts;

z = zero(sys_og)
p = pole(sys_og)
p_eject = p(end);
z_eject = z(end);
g_eject = zpk(p_eject, z_eject, (1-z_eject)/(1-p_eject), Ts)

p_eject2 = p(end-2:end-1);
z_eject2 = z(end-2:end-1);
g_eject2 = zpk(p_eject2, z_eject2, 1, Ts);
g_eject2 = g_eject2/dcgain(g_eject2)

sys = minreal(g_eject2*g_eject*sys_og);
G = sys;
%%
[Z_fft, freqs] = power_spectrum(ze_k, Ts);
figure(150)
semilogx(freqs, Z_fft)


%%

dist = timeseries(ze_k, t_k);
n = length(t_k);
% surface = zeros(n,1); 
% surface(n/2:end) = 0.05*0;
w_surf = 5*2*pi;
surface = sign(sin(w_surf*t_k))*0.05;

dist.Data = surface + sin(29.8*2*pi*t_k)*0.01; %dist.Data;
ref = timeseries(surface, t_k);
wo = 29.8*2*pi;

Dz = zpk(0, 1, 0.05, Ts);

trun = t_k(end)

g_filt = tf([1, 2*0.01*wo, wo^2], conv([1, 200], [1, 200]));
g_filtz = c2d(g_filt, Ts);
g_filtz = g_filtz/dcgain(g_filtz);
bode(g_filtz)
grid on
%%
sim('afm_z_dist_sin')

figure(10); clf
subplot(3,2,1)
plot(dist)
hold on
% plot(surf_est, '--')
title('30hz sin estimate')
% legend('truth', 'estimate')
grid
% subplot(3,2,3)
% hold on
% plot(y_hat1)
% title('y-hat')
% grid

subplot(3,2,5)
hold on
plot(uz.Time, surface)
title('surface')
grid

subplot(3,2,2)
plot(ze)
hold on
title('ze')

subplot(3,2,4)
plot(uz)
title('uz')
hold on
% plot(uz_filt)

%%
uz_filt = filtfilt(g_filtz.num{1}, g_filtz.den{1}, uz.Data);
% ze_filt = fft_notch(ze.Data, Ts, 29, 30);
% uz_filt = lsim(Dz, ze_filt, t_k);
figure(12)
plot(uz.Time, uz_filt)
ylim([-0.5, 0.5])
%%
fs = 25e3;             %#sampling rate
f0 = 29.8;                %#notch frequency
fn = fs/2;              %#Nyquist frequency
freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.

notchWidth = 0.1;       %#width of the notch

%Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%#Compute poles
notchPoles = (1-notchWidth) * notchZeros;

figure;
zplane(notchZeros.', notchPoles.');

b = poly( notchZeros ); %# Get moving average filter coefficients
a = poly( notchPoles ); %# Get autoregressive filter coefficients

figure;
freqz(b,a,32000,fs)

%#filter signal x
y = filter(b,a, uz.Data);

figure(6)
plot(y)
hold on
plot(uz.Data, '--')

%%
uzz = uz.Data;
close all
No=ceil(To/Ts);
To = Noo*Ts;
woo = 2*pi/Too;
tt = (0:No-1)'*Ts;
yy = cos(woo*tt);
M = floor(length(uzz)/No);
[y_fft, freqs] = power_spectrum(uzz(1:No*M), Ts);
semilogx(freqs, y_fft);
ze_filt = fft_notch(uzz(1:No*M), Ts, 29, 30);
plot(ze_filt)
%%
H = ones(No, 1);
H(2) = 0;
H(end) = 0;
overlap = No - 1;
N = 4*overlap;
step_size = N - overlap;

pos = 0;
y = uzz*0;
while pos+N < length(uzz)
  yt = ifft( fft( uzz(1+pos : N+pos),N ) .* H, N);
  
  y(1+pos : step_size+pos) = yt(M : N);
  pos = pos + step_size;
  
end


%%
[Z_fft, freqs] = power_spectrum(uz.Data, Ts);
figure(150)
semilogx(freqs, Z_fft)

[Z_fft, freqs] = power_spectrum(uz_filt, Ts);
hold on
semilogx(freqs, Z_fft, '--')




% figure(200), clf
% hold on
% plot(uz.Time, uz.Data, '--')
%%


% %%
% a_t = a - Lz*c;
% c_t = [sys.c*0, Gwz.c];
% H1 = ss(a_t, b, c_t, 0, Ts);
% H2 = ss(a_t, Lz, c_t, 0, Ts);
% 
% H_de = minreal(  (H2 -1)/( 1+ Dz*(G-H1-H2*G)));
% 
% figure(7)
% 
% bode(minreal(H_de*Dz))
% u = lsim(Dz*H_de, surface.Data, surface.Time);
% figure(200), clf
% plot(surface.Time, u)
% %%

function [Sens, H_dc_d, H_ud] = loops(G, Gw, LL, K, Cc, sys_obs)
  ns = size(G.a,1);
  nd = size(Gw.a,1);
  
  Lx = LL(1:ns);
  Ld = LL(ns+1:end);
  Nbar = SSTools.getNbar(G, K);
  [A, B, C] = ssdata(G);
  [Ad, Bd, Cd] = ssdata(Gw);
  
  At = [A - B*K - Lx*C,  -Lx*Cd - B*Nbar*Cd;
       -Ld*C,             Ad - Ld*Cd];
  Lt = [Lx; Ld];
  Kt = [K, Cd*Nbar];
  
  g_1 = ss(At, Lt, Kt, 0, G.Ts);
  Loop = G*g_1;
  
  Sens = minreal(1/(1 + Loop));
  
  H_ud = minreal(-g_1*Sens);
  
  
  g_2 = ss(At, Lt, Cc, 0, G.Ts);
  
  H_dc_d = minreal(g_2*Sens);
  
end





