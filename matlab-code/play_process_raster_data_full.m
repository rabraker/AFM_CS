% This script reads in the full, raw data from a raster scan run. This data
% will be produced by the vi play-raster-scan.vi. That vi produces two .csv
% files. One contains the the pre-processed data that labview does in
% realtime for visualization purposed. For slow scans, using that data is
% sufficient. The vi also produces a csv file with all of the raw data. For
% faster scans, we want to use that data, so we can process it better.
% That's what this script does.
clc
clear
close all
addpath('functions')

if ispc
% dat_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\raster\'

dat_root = 'Z:\afm-cs\imaging\raster';
  % data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data'
  % dat_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\raster';
else
  % dat_root = '/home/arnold/gradschool/afm-cs/afm_imaging/data';
  dat_root = fullfile(PATHS.exp, 'imaging', 'raster')
end

dat_name = 'raster_scan_512pix_5mic_01Hz_out_10-4-2018-03.csv'; % Ki=0.01;
dat_name = 'raster_scan_512pix_5mic_01Hz_out_10-15-2018-09.csv';
parent_name = 'raster_scan_5mic_01Hz.csv';

sub_dir = '5microns';

parent_name = get_parent_name(dat_name, '_out_')
meta_name = strrep(parent_name, '.csv', '-meta.mat')


% dat_name = '';
% parent_name = 'raster_scan__256pix_20mic_1.25e-01Hz.csv';

dat_path = fullfile(dat_root, sub_dir, dat_name);
parent_path = fullfile(dat_root, sub_dir, parent_name);
meta_in_path = fullfile(dat_root, sub_dir,  meta_name);
meta_out_path = strrep(dat_path, '.csv', '-meta.mat');

tic

% load(meta_in_path);
load(meta_out_path);

datmat = csvread(dat_path);
parent_dat = csvread(parent_path);


xyref = reshape(parent_dat', 2, [])';
xref = xyref(:,1);
%

% figure(1); plot(datmat(:,1));
% ax1 = gca
% figure(2); plot(datmat(:,3));
% ax3 = gca;
% figure(3); plot(datmat(:,4));
% ax4 = gca;
% linkaxes([ax1, ax3, ax4])
Ts = 40e-6;
samps_per_period = size(parent_dat,1)/2 % twice as many in here for x & y.
samps_per_line = samps_per_period/2

nperiods = 512;
pix = nperiods;
datmat = datmat([1:nperiods*samps_per_period], :);

% visualize tracking error.
if 1
    np = 3;
    xx = datmat(:,1);
    x_np = xx(1:np*length(xref));
    x_np = x_np - min(x_np);
    figure(200); clf; hold on
    t = [0:1:length(xref)*np-1]'*Ts;
    xref_np = repmat(xref, np, 1);
    p1 = plot(t, xref_np*volts2microns);
    p1.DisplayName = '$x_{ref}$';
    p2 = plot(t, x_np*volts2microns);
    p2.DisplayName = '$x(k)$';

    ylm = ylim;
    ylim([0, ylm(2)+.1*ylm(2)])
    ylabel('x-dir [\mu m]')
    xlabel('time [s]')
    leg1 = legend([p1, p2]);
    set(leg1, 'FontSize', 14, 'interpreter', 'latex', 'orientation', 'horizontal')
    leg1.Position=[0.6436    0.8590    0.2611    0.0640];
end
%  Now try by actually using the x-y measured data;
%%

models = load('G_zdir.mat');
% 
LPF = models.G2*models.G2*models.G2*models.G2;
% Dz = zpk([0], [1], 0.05, Ts);
% G1 = zpk(models.G1*models.G2)
% Hz = G1
% clc

% bsFilt = designfilt('bandstopfir','FilterOrder',20, ...
%          'CutoffFrequency1',2000,'CutoffFrequency2',2200, ...
%          'SampleRate',25e3);
% fvtool(bsFilt);
clf
[z, p, k] = zpkdata(models.G1);
g = zpk(p, [.95, .945], 1, Ts);
gg = tf(g)/dcgain(g)
bode(gg)
%%
clf
% % % w = 214*2*pi;
% % % z = 0.001;
% % % gct = tf(w^2, [1, 2*z*w, w^2])
% % % g_z = c2d(gct, Ts, 'matched');
% % % [z, p] = zpkdata(g_z, 'v')
% % % gg2 = zpk(p, [.992, .993], 1, Ts)
% % % gg2 = gg2/dcgain(gg2)
% % w = 214*2*pi;
% % z = 0.0001;
% % gct = tf(w^2, [1, 2*z*w, w^2])
% % g_z = c2d(gct, Ts, 'matched');
% % [z, p] = zpkdata(g_z, 'v')
% % gg2 = zpk(p, [.95, .945], 1, Ts)
% % gg2 = gg2/dcgain(gg2)
% % 
% % bode(gg2)
% % grid
% % 
% % G = gg2; %*gg
models = load(fullfile(PATHS.sysid, 'x-axis_sines_infoFourierCoef_10-21-2018-03.mat'))
G = -models.modelFit.G_zdir;
bode(G), grid on

D = zpk([0], [1], 0.05, Ts);
H = minreal(ss((1+G*D)/D)*LPF);
figure(5)
step(H)

figure
bode(H)
% hold on
% plot(t, yy, '--')
%%
clc
width = 5;
% volts2micron = 50/10;
micron2pix = pix/width;
volts2pix = volts2microns * micron2pix;


% [pixmat2, pixelifsampled] = bin_raster_really_slow(datmat(:,[1,2,4]), pix, samps_per_period, volts2pix, (H));
[pixmat2, pixelifsampled] = bin_raster_really_slow(datmat(:,[1,2,4]), pix, samps_per_period, volts2pix);

% pixmat2=pixmat2(:,1:end);

% pixmat2 = detrend_sampled_plane(pixmat2, pixelifsampled)
thresh = (20/7)*(1/1000)*20;
pixmat2 = pixmat2 - mean(pixmat2(:));
F10 = figure(5+0); clf
ax1 = gca();
f11 = figure(6+0); clf
ax2 = gca();
lo = min(min(pixmat2));
hi = max(max(pixmat2));
imshow_dataview(pixmat2, [-thresh, thresh], ax1, ax2)

%%
axis('on')

Kiz = Cluster.raster_scan_params.PI_params.Ki;
freq = meta.raster_freq;
stit = sprintf('%.2fHz, Kiz = %f', freq, Kiz);

title(stit)
traster_recon = toc
fprintf('Total data processing time:%.3f', traster_recon)

%%
fig_path = strrep(dat_path, '.csv', '-fig.fig');
saveas(F10, fig_path)
