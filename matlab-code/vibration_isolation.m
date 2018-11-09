% A quick script to check out the vibration isolation chamber I built. The
% first section computes the a rough estimate of the spring constant and
% natural frequency. The second section plots the new and old data with the
% cantilevar just sitting on the surface. THis is a great success, at least
% when the xy-stage is turned off. 
%
% Oddly, when the xy-stage is turned on, we get a very visible sinusoid at
% about 360 hz. I don't know what to make of this. 
% 

W = 47.4; % lbs
l = 11; % inches

L = l * 2.54 / 100;
 
lbs2newtons = 4.4822162825086;
lbs2kilos = 0.45359237; 
M_kilos = W * lbs2kilos;

F_newtons = W * lbs2newtons;

% F = k * x
K = F_newtons / L;

wo = sqrt( K / M)
fo = wo/2/pi
%%
clc
% dat1 = csvread('/media/labserver/afm-cs/z-scope/data-out_KI_p00_nobungee_10-26-2018.csv');
root = fullfile(PATHS.exp, 'z-scope')
dat1 = csvread(fullfile(root, 'data-out_KI_p00_nobungee_10-30-2018.csv'));
dat2 = csvread(fullfile(root, 'data-out_KI_p00_withbungee_10-30-2018.csv'));
dat3 = csvread(fullfile(root, 'data-out_KI_p00_withbungee_xyon_disconnected_10-30-2018.csv'));
dat4 = csvread(fullfile(root, 'data-out_KI_p00_nobungee_xyon_connected_10-30-2018.csv'));

dat1 = dat1(4000:end,:);
dat2 = dat2(4000:end,:);
dat3 = dat3(4000:end,:);
dat4 = dat4(4000:end,:);

N1 = size(dat1,1);
N2 = size(dat2,1);
N3 = size(dat3,1);
N4 = size(dat4,1);

N = min([N1, N2, N3, N4]);

ze1 = detrend(dat1(end-N+1:end,1));
ze2 = detrend(dat2(end-N+1:end,1));
ze3 = detrend(dat3(end-N+1:end,1));
ze4 = detrend(dat4(end-N+1:end,1));



figure(8); clf

h1 = plot(ze1);
h1.DisplayName = 'no-bungee, xy-off';
hold on;
h2 = plot(ze2);
h2.DisplayName = 'with bungee, xy-off';
h3 = plot(ze3);
h3.DisplayName = 'with bungee, xy-on, discon';
h4 = plot(ze4);
h4.DisplayName = 'no bungee, xy-on, con';

legend([h1, h2, h3, h4])

[Z1, freqs] = power_spectrum(ze1, Ts);
[Z2, freqs] = power_spectrum(ze2, Ts);
[Z3, freqs] = power_spectrum(ze3, Ts);
[Z4, freqs] = power_spectrum(ze4, Ts);

figure(9); clf
semilogx(freqs, Z1);
hold on
semilogx(freqs, Z2)
semilogx(freqs, Z3)
% semilogx(freqs, Z4)

legend([h1, h2, h3])

%%
clc
% dat1 = csvread('/media/labserver/afm-cs/z-scope/data-out_KI_p00_nobungee_10-26-2018.csv');
root = fullfile(PATHS.exp, 'z-scope')
dat1 = csvread(fullfile(root, 'data-out_KI_p00_xyon_nondisconnected.csv'));
dat2 = csvread(fullfile(root, 'data-out_KI_p00_xon_ydisconnected.csv'));
dat3 = csvread(fullfile(root, 'data-out_KI_p00_yon_xdisconnected.csv'));
% dat2 = csvread(fullfile(root, 'data-out_KI_p00_xyon_disconnected.csv'));

dat1 = dat1(4000:end,:);
dat2 = dat2(4000:end,:);
dat3 = dat3(4000:end,:);

N1 = size(dat1,1);
N2 = size(dat2,1);
N3 = size(dat3,1);


ze1 = detrend(dat1(:,1));
ze2 = detrend(dat2(:,1));
ze3 = detrend(dat3(:,1));

% N = length(dat1(:,1));
N = min([N1, N2, N3]);

figure(8); clf

h1 = plot(ze1(1:N));
h1.DisplayName = 'all connected';
hold on;
h2 = plot(ze2(1:N));
h2.DisplayName = 'xon, ydisconnected';
h3 = plot(ze3(1:N));
h3.DisplayName = 'yon, x-disconnected';

[Z1, freqs] = power_spectrum(ze1(1:N), Ts);
[Z2, freqs] = power_spectrum(ze2(1:N), Ts);
[Z3, freqs] = power_spectrum(ze3(1:N), Ts);
legend([h1, h2, h3])
%%
clc
figure(9); clf


h11 = semilogx(freqs, log10(Z1));
h11.DisplayName = 'xy-on, both connected';
hold on, grid on
h22 = semilogx(freqs, log10(Z2));
h22.DisplayName = 'xy-on, y-disconnected';

h33 = semilogx(freqs, log10(Z3));
h33.DisplayName = 'xy-on, x-disconnected';
% uistack(h11, 'top')
legend([h11, h22, h33]);

%%
rmpath functions
addpath ~/matlab/afm_mpc_journal/functions/canon/
addpath ~/matlab/afm_mpc_journal/functions

[plants_xy, frf_xy] = CanonPlants.plants_ns14(9, 2);

rmpath ~/matlab/afm_mpc_journal/functions/canon/
rmpath ~/matlab/afm_mpc_journal/functions/
addpath functions

rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.mat');
models_z = load(rand_fname);

yyaxis right
h44 = semilogx(frf_xy.freqs_Hz, log10(abs(frf_xy.G_uz2stage)), '-k', 'LineWidth', 2);
h44.DisplayName = '$G_{z,u_z}$';
h55 = semilogx(models_z.modelFit.frf.freqs_Hz, log10(abs(models_z.modelFit.frf.G_uz2stage)), '-r', 'LineWidth', 2);
h55.DisplayName = '$G_{z,u_z}$';

legend([h11, h22, h33, h44, h55])
%%

dat = csvread(fullfile(root, 'zscope_longseries_ki_0_11-7-2018_mica_01.csv'));


Z_ = dat(:,1);
figure(8)
xlm = xlim;
plot(Z_)
xlim(xlm)
%%
N = 1024*2^3;
M = floor(size(Z_, 1)/N);

Z = Z_(1:N*M);


Z_mat = reshape(Z, [], M);


zdft_mag_sum = zeros(N/2, 1);
zdft_mag2_sum = zeros(N/2, 1);

for k=1:M

  [zdft_k, freqs] = power_spectrum(Z_mat(:, k), Ts);

  zdft_k(1) = [];
  freqs(1) = [];
  zdft_mag_k = abs(zdft_k);

  zdft_mag_sum = zdft_mag_sum + zdft_mag_k;
 
end

zdft_mean = zdft_mag_sum/M;

for k=1:M
   zdft_mag2_sum = zdft_mag2_sum + (zdft_mag_k - zdft_mean).^2;
end


zdft_cov = zdft_mag2_sum/(M-1);
zdft_stdd = sqrt(zdft_cov);


figure(10); %clf
semilogx(freqs, log10(zdft_mean), '-k')
semilogx(freqs, log10(abs(zdft_mean-zdft_stdd)), '--r')
semilogx(freqs, log10(abs(zdft_mean+zdft_stdd)), '--r')


hold on
yyaxis right
h44 = semilogx(frf_xy.freqs_Hz, log10(abs(frf_xy.G_uz2stage)), '-k', 'LineWidth', 2);

%%
clc
L = N;
[pxx, w, pxxc] = pwelch(Z,hamming(L),50,1024,1/Ts, 'ConfidenceLevel',0.97);
pxx(1) = []; w(1) = []; pxxc(1, :) = [];
lower = log10(pxxc(:,1));
upper = log10(pxxc(:,2));

figure(200); clf

h = ciplot(lower, upper, w,'b');
hold on
plot(w, log10(pxx), 'k')
alpha(h, '.5')
grid
ax = gca();
ax.XScale = 'log';
%%


hold on

%%

%%

semilogx(w, log10(pxxc(:,1)), '--r')
semilogx(w, log10(pxxc(:,2)), '--r')
%%
frf_root = [PATHS.sysid, '/multi-axis'];
ls(frf_root)

start_row = 2;
ux_xdir_dat = csvread(fullfile(frf_root, 'ux_to_xdir.csv'), start_row);
ux_ydir_dat = csvread(fullfile(frf_root, 'ux_to_ydir.csv'), start_row);
ux_zdfl_dat = csvread(fullfile(frf_root, 'ux_to_zdfl.csv'), start_row);

uy_ydir_dat = csvread(fullfile(frf_root, 'uy_to_ydir.csv'), start_row);
uy_xdir_dat = csvread(fullfile(frf_root, 'uy_to_xdir.csv'), start_row);
uy_zdfl_dat = csvread(fullfile(frf_root, 'uy_to_zdfl.csv'), start_row);

uz_zdfl_dat = csvread(fullfile(frf_root, 'uz_to_zdfl.csv'), start_row);


%%
clc
[G_ux_x, freqs, C_ux_x] = process_rand_frf(ux_xdir_dat);
[G_ux_y, freqs, C_ux_y] = process_rand_frf(ux_ydir_dat);
[G_ux_z, freqs, C_ux_z] = process_rand_frf(ux_zdfl_dat);

[G_uy_y, freqs, C_uy_y] = process_rand_frf(uy_ydir_dat);
[G_uy_x, freqs, C_uy_x] = process_rand_frf(uy_xdir_dat);
[G_uy_z, freqs, C_uy_z] = process_rand_frf(uy_zdfl_dat);

[G_uz_zdfl, freqs, C_uz_z] = process_rand_frf(uz_zdfl_dat);

F = figure(100); clf

ax1 = subplot(3,3,1);
ax2 = subplot(3,3,2);

ax3 = subplot(3,3,3);
ax4 = subplot(3,3,4);
ax5 = subplot(3,3,5);

ax7 = subplot(3,3,7);
ax8 = subplot(3,3,8);
ax9 = subplot(3,3,9);

yyaxis(ax1, 'right')
semilogx(ax1, freqs, C_ux_x)
%%
yyaxis(ax1, 'left')
h1= frf_bode_mag(G_ux_x, freqs, ax1, 'Hz', '-k');
%%

frf_bode_mag(G_ux_y, freqs, ax4, 'Hz', '-k');
yyaxis(ax4, 'right')
semilogx(ax4, freqs, C_ux_y)

frf_bode_mag(G_ux_z, freqs, ax7, 'Hz', '-k');
yyaxis(ax7, 'right')
semilogx(ax7, freqs, C_ux_z)

frf_bode_mag(G_uy_x, freqs, ax2, 'Hz', '-k');
yyaxis(ax2, 'right')
semilogx(ax2, freqs, C_uy_x)

frf_bode_mag(G_uy_y, freqs, ax5, 'Hz', '-k');
yyaxis(ax5, 'right')
semilogx(ax5, freqs, C_uy_y)

frf_bode_mag(G_uy_z, freqs, ax8, 'Hz', '-k');
yyaxis(ax8, 'right')
semilogx(ax8, freqs, C_uy_z)

frf_bode_mag(G_uz_zdfl, freqs, ax9, 'Hz', '-k');
yyaxis(ax9, 'right')
semilogx(ax9, freqs, C_uz_z);

%%
title(ax1, 'Control: $u_x$')
title(ax2, 'Control: $u_y$')
title(ax3, 'Control: $u_z$')
ax3.XTick = []
ax3.YTick = []
%%
ylabel(ax1, 'to $x$-dir', 'FontSize', 16);
ylabel(ax4, 'to y-dir', 'FontSize', 16);
ylabel(ax7, 'to $z$-dfl', 'FontSize', 16);

% yl = ylabel(ax2, 'to $y$-dir');
% yl = ylabel(ax4, 'to $z$-dfl');


%%









function [G, freqs, Coh] = process_rand_frf(dat)
  
  freqs = dat(:,1);
  mags_db = dat(:,2);
  phase_deg = dat(:,4);
  Coh = dat(:, 6);
  G = (10.^(mags_db/20)).*exp(1j * phase_deg * pi / 180);  
end