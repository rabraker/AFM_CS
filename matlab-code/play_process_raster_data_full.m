% This script reads in the full, raw data from a raster scan run. This data
% will be produced by the vi play-raster-scan.vi. That vi produces two .csv
% files. One contains the the pre-processed data that labview does in
% realtime for visualization purposed. For slow scans, using that data is
% sufficient. The vi also produces a csv file with all of the raw data. For
% faster scans, we want to use that data, so we can process it better.
% That's what this script does.
clc

init_paths();

clear size
size_dir = '5microns';
dat_root = PATHS.raster_image_data(size_dir, '11-13-2018');
dat_name = 'raster_scan_512pix_5mic_01Hz_out_11-13-2018-01.csv';

raster_paths = get_raster_paths(dat_root, dat_name)
%%

% meta_name = strrep(parent_name, '.csv', '-meta.mat');
% fprintf('Using parent Data file\n%s\n', parent_name)

% dat_path = fullfile(dat_root, sub_dir, dat_name);
% parent_path = fullfile(dat_root, sub_dir, parent_name);
% meta_in_path = fullfile(dat_root, sub_dir,  meta_name);
% meta_out_path = strrep(dat_path, '.csv', '-meta.mat');


Ts = 40e-6;
npix = 512;
rast_exp = RasterExp(raster_paths, npix, width);


% fprintf('Loading Meta file...\n%s\n', raster_paths.meta_path)
% load(raster_paths.meta_path);
% fprintf('Loading Data file...\n%s\n', raster_paths.data_path)
% datmat = csvread(raster_paths.data_path);
% 
% parent_dat = csvread(raster_paths.parent_path);
% xyref = reshape(parent_dat', 2, [])';
% xref = xyref(:,1);
% samps_per_period = size(parent_dat,1)/2 % twice as many in here for x & y.
% samps_per_line = samps_per_period/2
% 
% datmat = datmat([1:npix*samps_per_period], :);

% figure(1); plot(datmat(:,1));
% ax1 = gca
% figure(2); plot(datmat(:,3));
% ax3 = gca;
% figure(3); plot(datmat(:,4));
% ax4 = gca;
% linkaxes([ax1, ax3, ax4])

%%
% visualize tracking error.

if 1
    np = 3;
    
    x_np = rast_exp.x(1:np*length(rast_exp.xref));
    x_np = x_np - min(x_np);
    figure(200); clf; hold on
    t = [0:1:length(rast_exp.xref)*np-1]'*rast_exp.Ts;
    xref_np = repmat(rast_exp.xref, np, 1);
    p1 = plot(t, xref_np*volts2microns);
    p1.DisplayName = '$x_{ref}$';
    p2 = plot(t, x_np*volts2microns);
    p2.DisplayName = '$x(k)$';

    ylm = ylim;
    ylim([0, ylm(2)+.1*ylm(2)])
    ylabel('x-dir [$\mu$ m]')
    xlabel('time [s]')
    leg1 = legend([p1, p2]);
    set(leg1, 'FontSize', 14, 'interpreter', 'latex', 'orientation', 'horizontal')
    leg1.Position=[0.6436    0.8590    0.2611    0.0640];
end

figure(2)
subplot(2,1,1)
ax1 = gca();
N1 = 314;
N2 = N1+4;
rast_exp.plot_n_periods('ze', gca, N1, N2)
% plot(datmat(N1*samps_per_line+1:N2*samps_per_line,3))
subplot(2,1,2)
rast_exp.plot_n_periods('uz', gca, N1, N2)
% plot(datmat(N1*samps_per_line+1:N2*samps_per_line,4))
ax2 = gca();

linkaxes([ax1, ax2], 'x')
%  Now try by actually using the x-y measured data;
%%


w_lpf = 600*2*pi;
z_lpf = exp(-w_lpf*Ts);
LPF1 = zpk([], [z_lpf], 1-z_lpf, Ts);
LPF = LPF1*LPF1;
rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.mat');
models = load(rand_fname);
G1 = models.modelFit.G_zdir;
Dinv1 = models.modelFit.Dinv;
gdrift1 = models.modelFit.gdrift;

KI = -Cluster.raster_scan_params.PI_params.Ki;
D1 = zpk([0], 1, KI, G1.Ts) * Dinv1;
H = -minreal(feedback(D1, G1))*G1*LPF;
figure, step(H)


zz = lsim(gdrift1, datmat(:,4), (0:length(datmat(:,4))-1)'*G1.Ts);
figure(1); clf
% subplot(2,1,1)
%%
plot(detrend2(detrend(zz))); %(1:1:samps_per_line*8))
hold on

% subplot(2,1,2)
% plot(datmat(1:samps_per_line*80,4))

%%
figure(20)
plot(datmat(1:samps_per_line*4, 3))
%%
idd = iddata(datmat(1:samps_per_line, 3), datmat(1:samps_per_line, 4), Ts);

sys = ssest(idd, 12);

zes = lsim(sys, datmat(1:samps_per_line, 3), (0:samps_per_line-1)'*Ts)

plot(zes)

%%
clc

% [pixmat2, pixelifsampled] = bin_raster_really_slow([datmat(:,[1,2]), zz], pix, samps_per_period, volts2pix);
rast_exp.bin_raster_really_slow();

pixmat2 = rast_exp.pix_mat;
% x = [1:512];
% y = x;
% [xx, yy, zz] = prepareSurfaceData(x, y, pixmat2);
% %%
% pixmat3 = pixmat2*0;
% f= fit([xx, yy], zz, 'poly23')
% for k=1:512
%   pixmat3(k, :) = pixmat2(k,:) - f(k, x');
% end
% [xx, yy, zzz] = prepareSurfaceData(x, y, pixmat3);
% figure, plot([xx, yy], zzz)



thresh = (20/7)*(1/1000)*20;
pixmat2 = pixmat2 - mean(pixmat2(:));
F10 = figure(5+2); clf
ax1 = gca();
f11 = figure(6+2); clf
ax2 = gca();
lo = min(min(pixmat2));
hi = max(max(pixmat2));


pixmat3 = detrend_plane(detrend2(pixmat2(10:end-10,10:end-25)));
imshow_dataview(pixmat3 - mean(pixmat3(:)), [-thresh, thresh], ax1, ax2)

figure,
mesh(pixmat3)
% [xx, yy, zzz] = prepareSurfaceData(x, y, pixmat3);
% figure, plot([xx, yy], zzz)
%%
save('tuesday-figs/11-5-2018/derate_retrace.mat', 'pixmat2', 'thresh', 'Cluster')
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
