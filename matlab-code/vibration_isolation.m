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
%%
yyaxis right
h44 = semilogx(frf_xy.freqs_Hz, log10(abs(frf_xy.G_uz2stage)), '-k', 'LineWidth', 2);
h44.DisplayName = '$G_{z,u_z}$';
h55 = semilogx(models_z.modelFit.frf.freqs_Hz, log10(abs(models_z.modelFit.frf.G_uz2stage)), '-r', 'LineWidth', 2);
h55.DisplayName = '$G_{z,u_z}$';

legend([h11, h22, h33, h44, h55])