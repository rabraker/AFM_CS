
% clear, clc
Ts = 40e-6;


Dz = zpk([0], [1], 0.05, Ts);

models = load('G_zdir.mat')
G2 = models.G2;

G1 = zpk(models.G1*models.G2)
models = load(fullfile(PATHS.sysid, 'x-axis_sines_infoFourierCoef_10-21-2018-03.mat'))
G1 = -models.modelFit.G_zdir;
models.G2 = G2;
%%
figure(1); clf
bode(G1, Dz*G1)


r = 0;

f = 10;
To = 1/f;
n = 10;

N = floor(To*n/Ts)
t = (0:N-1)'*Ts;
d_ = sin(2*pi*f*t);
surface = timeseries(sign(d_), t);
trun = t(end);

% noise
eta_ = mvnrnd(0, 0.001, size(t,1));
eta = timeseries(eta_, t);


sim('afm_z_simple')

figure(2)
subplot(2,1,1)
plot(ze)
title('ze')
subplot(2,1,2)
plot(uz)
title('uz')

H2 = 1 + G1*Dz;
H1 = minreal(H2/Dz);

% ze_lpf = lsim(models.G2*models.G2, ze.Data, ze.Time);
LPF = models.G2*models.G2*models.G2*models.G2;
surf_est = lsim(H2*LPF, ze.Data, ze.Time);

surf_est2 = lsim(-H1*LPF, uz.Data, uz.Time);

figure(3); clf
subplot(2,1,1)
plot(surface)
hold on
plot(t, surf_est, '--')
title('(1+D(z)G(z))*ze(z)')

subplot(2,1,2),
plot(surface)
hold on
plot(t, surf_est2, '--')
title('(1+D(z)G(z)/D(z))*ze(z)')

%%
dat_root = fullfile(PATHS.exp, 'imaging', 'raster')
dat_name = 'raster_scan_512pix_5mic_01Hz_out_10-15-2018-09.csv';
parent_name = 'raster_scan_5mic_01Hz.csv';

sub_dir = '5microns';

parent_name = get_parent_name(dat_name, '_out_')
meta_name = strrep(parent_name, '.csv', '-meta.mat')


dat_path = fullfile(dat_root, sub_dir, dat_name);
parent_path = fullfile(dat_root, sub_dir, parent_name);
meta_in_path = fullfile(dat_root, sub_dir,  meta_name);
meta_out_path = strrep(dat_path, '.csv', '-meta.mat');

load(meta_out_path);

datmat = csvread(dat_path);
parent_dat = csvread(parent_path);

Ts = 40e-6;
samps_per_period = size(parent_dat,1)/2 % twice as many in here for x & y.
samps_per_line = samps_per_period/2

nperiods = 512;
pix = nperiods;
datmat = datmat([1:nperiods*samps_per_period], :);

%%
xy = datmat(:,1:2);
uz = datmat(:, 4);
clear datmat;
%%
t = (0:length(uz)-1)'*Ts;
surf_est2 = lsim(H1*LPF, uz, t);

width = 5;
% volts2micron = 50/10;
micron2pix = pix/width;
volts2pix = volts2microns * micron2pix;
%%
[pixmat2, pixelifsampled] = bin_raster_really_slow([xy, uz], pix, samps_per_period, volts2pix);


thresh = (20/7)*(1/1000)*20;
pixmat2 = pixmat2 - mean(pixmat2(:));
F10 = figure(5); clf
ax1 = gca();
f11 = figure(6); clf
ax2 = gca();
lo = min(min(pixmat2));
hi = max(max(pixmat2));
imshow_dataview(pixmat2, [-thresh, thresh], ax1, ax2)
%%
save('standard_raster.mat', 'pixmat2', 'thresh')