% Compare raster scanning on the flats to CS-scanning on the flats. How much
% extra oscillation to we actually get?

% initialize paths.
% init_paths();
fig_root = fullfile(PATHS.CS_root, 'matlab-code', 'notes', 'figures');
Ts = 40e-6;

img_size = '5microns';



%move   |  lower  |  settle  | scan   | up 
% 4.035 | 0.971  | 9.658   | 22.607  | 0.428 |
% Total Imaging time: 37.70
cs_exp_data_name_s{1} = 'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-24-2018-03.csv';
data_root = PATHS.cs_image_data(img_size, '11-24-2018');


chan_map = ChannelMap([1:5]);
% -------------------
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);

gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);

cs_paths = get_cs_paths(data_root, cs_exp_data_name_s{1});
hole_depth = (20);
cs_exp = CsExp(cs_paths, chan_map, Ts, hole_depth, gg);
cs_exp.print_state_times();
fprintf('Total Imaging time: %.2f\n', cs_exp.time_total)

% % % cs_exp.process_cs_data(false);
% % % cs_exp.solve_basis_pursuit()
% % % pixmat_cs = cs_exp.Img_bp;
% % % %%
% % % F12 = mkfig(7, 4, 6, true); clf
% % % ax1 = axes('Position', [0.11438 0.4513 0.775 0.51537]);
% % % ax2 = axes('Position', [0.13 0.11 0.775 0.27194]);
% % % ax2.Color = fig_color;
% % % 
% % % imshow_dataview(pixmat_cs, [-thresh, thresh], ax1, ax2);
% % % 
% % % row2_1 = 150;
% % % row2_2 = 175;
% % % hold(ax1, 'on')
% % % plot(ax1, [1,512], [row2_1, row2_1], 'r')
% % % plot(ax1, [1,512], [row2_2, row2_2], 'b')
% % % 
% % % plot(ax2, [1:512], pixmat_cs(512-row2_1, :), 'r')
% % % hold(ax2, 'on')
% % % plot(ax2, [1:512], pixmat_cs(512-row2_2, :), 'b')
% % % ylim(ax2, [-thresh, thresh])
% % % U_cs = 0;
% % % N = 1;

%%
figbase = 20;

[~, axs] = make_traj_figs(figbase);
cs_exp.plot_all_cycles(axs{:});


%% Load Raster and plot
clc

raster_root = PATHS.raster_image_data(img_size, '11-25-2018');
dat_name = 'raster_scan_512pix_5mic_01Hz_out_11-25-2018bungee_extfan-01.csv';

npix = 512;
raster_paths = get_raster_paths(raster_root, dat_name);

rast_exp = RasterExp(raster_paths, npix, width)
%%

% thresh = (20/7)*(1/1000)*20;

clc
rast_exp.bin_raster_really_slow(@detrend)
F11 = mkfig(7, 4, 6, true); clf
ax1 = axes('Position', [0.11438 0.4513 0.775 0.51537]);
ax2 = axes('Position', [0.13 0.11 0.775 0.27194]);
ax2.Color = fig_color;

imshow_dataview(rast_exp.pix_mat, [-thresh, thresh], ax1, ax2);
%%
row2_1 = 155;
row2_2 = 175;
hold(ax1, 'on')
plot(ax1, [1,512], [row2_1, row2_1], 'r')
plot(ax1, [1,512], [row2_2, row2_2], 'b')

[uz2_1, x2_1] = rast_exp.get_row(row2_1);
[uz2_2, x2_2] = rast_exp.get_row(row2_2);

[uz1_1, x1_1] = rast_exp.get_row(row2_1);
[uz1_2, x1_2] = rast_exp.get_row(row2_2);

hold(ax2, 'on');
plot(ax2, x1_1(1:end), uz1_1(1:end), 'r');
plot(ax2, x1_2(1:end), uz1_2(1:end), 'b');


xlim(ax2, [0, 512])
ylim(ax2, [-thresh, thresh])
grid(ax2, 'on')

%%
% Take the raster data and subsample it. Then reconstruct.
sub_samp_rast = rast_exp.pix_mat.*cs_exp.pix_mask;
[n m] = size(sub_samp_rast);
tic
I_vector = PixelMatrixToVector(sub_samp_rast);

pix_mask_vec = PixelMatrixToVector(cs_exp.pix_mask);
I_vector = I_vector(find(pix_mask_vec>0.5));

A = @(x) IDCTfun(x,pix_mask_vec);
At = @(x) DCTfun(x,pix_mask_vec);

Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
Ir_bp = real(Ir_bp);
Img_bp_from_raster = PixelVectorToMatrix(Ir_bp,[n m]);
time_bp = toc;
fprintf('BP Time: %f\n', time_bp);

figure(2)
ax1 = subplot(3,1,1:2)
ax2 = subplot(3,1,3);
imshow_dataview(Img_bp_from_raster - mean(Img_bp_from_raster(:)), [-thresh, thresh], ax1, ax2)

figure(3)
ax3 = subplot(3,1,1:2)
ax4 = subplot(3,1,3);
imshow_dataview(cs_exp.Img_bp, [-thresh, thresh], ax3, ax4)

stit = sprintf('Experimental reconstruction. Total time = %.1f', cs_exp.time_total)
title(ax3, stit)
%%
save('tuesday-figs/11-27-2018/img_data.mat', 'cs_exp', 'pixmat2', 'thresh', 'Img_bp_from_raster',...
  'sub_samp_rast')
%%
figure(4)
ax1 = subplot(3,1,1:2)
ax2 = subplot(3,1,3);
imshow_dataview(sub_samp_rast - mean(sub_samp_rast(:)), [-thresh, thresh]*.5, ax1, ax2)
%%

pixmat3 = pixmat2;
align_y_idx = 220; 
u_floor = mean(detrend(pixmat3(:, align_y_idx)));
for row = 1:512
  u_k_yidx = pixmat3(row, align_y_idx);
  pixmat3(row, :) = pixmat3(row, :) - u_k_yidx;
end
pixmat3 = detrend_plane(pixmat3);
figure(5)
ax1 = subplot(3,1,1:2)
ax2 = subplot(3,1,3);
imshow_dataview(pixmat3 - mean(pixmat3(:)), [-thresh, thresh]*.5, ax1, ax2)

%%
% compare deviations in flat CS cycle to the raster scan over flat
starts = [9.0, 17.0];
ends = [10.63, 18.65];
[U_cs_psd2, freqs2] = cs_exp.psd_from_intervals('uz', 'scan', starts, ends);

U_rast = 0;
N = 1;
win = window(@hann, samps_per_line+1);
for k=[147:155, 197:210] % rows to take the psd from. These *should* be over flat.
  [uz_row] = get_row(datmat2, k, samps_per_line, npix);
  [U_rast_k, freqs_rast] = power_spectrum(uz_row, Ts, win);
%   [U_rast_k, freqs_rast] = pwelch(uz_row, win, [],[],1/Ts);
  U_rast = U_rast + U_rast_k;
  N = N+1;
end
U_rast = U_rast/N;

f12 = figure(12); clf
ax = gca();
% [U_rast, f] = pwelch(uz1_1, [], [],[],1/Ts);

h1 = semilogx(ax, freqs_rast, 10*log10(U_rast));
hold(ax, 'on')
h2 = semilogx(ax, freqs, 10*log10(U_cs_psd));
h3 = semilogx(ax, freqs2, 10*log10(U_cs_psd2));

grid on

xlim([10, 12500])
xlabel('frequency [Hz]')
h1.DisplayName = sprintf('Raster (%d lines)', N);
h2.DisplayName = sprintf('CS (%d scans)', size(uz_mat,2));
leg = legend([h1, h2]);
set(leg, 'location', 'SouthWest')

%%
if 0
  root_psd = '/media/labserver/afm-cs/z-scope/psd-test-fan';
  fname_ext_psd = 'xz_engaged_frf_externalpower.json';
  
  fpath_ext = fullfile(root_psd, fname_ext_psd);
  psd_ext = LvPSD(fpath_ext);
  
  figure(12);
  h_psd = semilogx(psd_ext.freqs, psd_ext.psd(:,3), 'g');
  h_psd.DisplayName = 'No Motion'

end

% saveas(f12, 'tuesday-figs/11-27-2018/psd_scans_01.fig', 'fig')

if 1
  yyaxis right
  % Load frf data
  z_frf = RandNoiseFRF('/media/labserver/afm-cs/sysID/z-axis_frf_zero-offset.json');

  y_frf = RandNoiseFRF('/media/labserver/afm-cs/sysID/multi-axis/ydrive_frf_11-10-2018.csv');
  x_frf = RandNoiseFRF('/media/labserver/afm-cs/sysID/multi-axis/xdrive_frf_11-10-2018.csv');

  hf_x = frf_bode_mag(x_frf.G_frd(1,1), x_frf.freqs, f12, 'Hz', 'k');
  hf_y = frf_bode_mag(y_frf.G_frd(1,1), y_frf.freqs, f12, 'Hz', '--k');
  hf_z = frf_bode_mag(z_frf.G_frd(1,1), z_frf.freqs, f12, 'Hz', ':k');
  
  hf_x.DisplayName = 'xdir FRF';
  hf_y.DisplayName = 'ydir FRF';
  hf_z.DisplayName = 'zdir FRF';
  ax = gca();
  ax.YColor = 'k';
  saveas(f12, 'tuesday-figs/11-27-2018/psd_scans_02.fig', 'fig')
end
% save_fig(f12, fullfile(fig_root, 'CS_raster_flats_PSD_withFRF'), false);
%%
% Now, compare this to a z-bounce data set. Exactly the same control parameters
% and tip.
% This one is before I re-designed the Dinv.
cs_paths_zb1 = get_cs_paths(fullfile(PATHS.exp, 'z-bounce', '11-26-2018'),...
 'cs-traj-z-bounce_500nmequiv_out_11-26-2018-05.csv');
gg2 = zpk([], [], 1, Ts);
cs_exp_zb1 = CsExp(cs_paths_zb1,  chan_map, Ts, hole_depth, gg2);
cs_exp_zb1.process_cs_data(false)
%%
clc
clear CsExp
[U_zb1, freqs_zb1] = cs_exp_zb1.psd_from_intervals('uz', 'scan', 1,7);


[~, axs] = make_traj_figs(figbase+100);
cs_exp_zb1.plot_all_cycles(axs{:});

h_zb1 = semilogx(ax, freqs_zb1, 10*log10(U_zb1));
h_zb1.DisplayName = 'Zbounce, before redesign';
% saveas(f12, 'tuesday-figs/11-27-2018/psd_scans_03.fig', 'fig')
%%
% Now the same thing, but after I re-ID Dinv
cs_paths_zb2 = get_cs_paths(fullfile(PATHS.exp, 'z-bounce', '11-26-2018'),...
 'cs-traj-z-bounce_500nmequiv_out_11-26-2018newDinv-03.csv');
gg2 = zpk([], [], 1, Ts);
cs_exp_zb2 = CsExp(cs_paths_zb2,  chan_map, Ts, hole_depth, gg2);
cs_exp_zb2.process_cs_data(false)
[U_zb2, freqs_zb2] = cs_exp_zb2.psd_from_intervals('uz', 'scan', 1,7);

%%
[~, axs] = make_traj_figs(figbase+100);
cs_exp_zb2.plot_all_cycles(axs{:});

figure(12)
h_zb2 = semilogx(ax, freqs_zb2, 10*log10(U_zb2));
h_zb2.DisplayName = 'Zbounce, after redesign';
% saveas(f12, 'tuesday-figs/11-27-2018/psd_scans_03.fig', 'fig')



cs_paths_zb3 = get_cs_paths(fullfile(PATHS.exp, 'z-bounce', '11-26-2018'),...
 'cs-traj-z-bounce_500nmequiv_out_11-26-2018noDinv-01.csv');
gg2 = zpk([], [], 1, Ts);
cs_exp_zb3 = CsExp(cs_paths_zb3,  chan_map, Ts, hole_depth, gg2);
cs_exp_zb3.process_cs_data(false)
[U_zb3, freqs_zb3] = cs_exp_zb3.psd_from_intervals('uz', 'scan', 1,7);

[~, axs] = make_traj_figs(figbase+100);
cs_exp_zb3.plot_all_cycles(axs{:});

figure(12)
h_zb3 = semilogx(ax, freqs_zb3, 10*log10(U_zb3));
h_zb3.DisplayName = 'Zbounce, NO DINV';
% saveas(f12, 'tuesday-figs/11-27-2018/psd_scans_04.fig', 'fig')


%

% --------------------------------------------------------------------
% Full CS imaging cycle, AFTER Dinv re-fit.
droot2 = PATHS.cs_image_data('5microns', '11-26-2018');
fname = 'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-26-2018newDinv-02.csv';
cs_paths_newDinv = get_cs_paths(droot2, fname);

cs_exp_ND2 = CsExp(cs_paths_newDinv,  chan_map, Ts, hole_depth, gg);
cs_exp_ND2.process_cs_data(false)


[~, axs] = make_traj_figs(figbase+100);
% cs_exp_zbounce.plot_all_cycles(axs{:});

starts = [17.5, 21.22, 24];
ends = [19.5, 22.45, 25.5];
[U_nd2, freqs_nd2, k_nd2] = cs_exp_ND2.psd_from_intervals('uz', 'scan', starts, ends);

figure(12)
h_nd2 = semilogx(freqs_nd2, 10*log10(U_nd2), 'r');
h_nd2.DisplayName = sprintf('CS After Redesign(%d scans)', k_nd2);
% saveas(f12, 'tuesday-figs/11-27-2018/psd_scans_05.fig', 'fig')

 
% fname3 = 'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-26-2018newDinv-03.csv';
% cs_paths_newDinv3 = get_cs_paths(droot2, fname3);
% 
% cs_exp_ND3 = CsExp(cs_paths_newDinv3,  chan_map, Ts, hole_depth, gg);
% cs_exp_ND3.process_cs_data(false)
% 
% 
% [~, axs] = make_traj_figs(figbase+100);
% cs_exp_ND3.plot_all_cycles(axs{:});
% %%%
% starts = [7.247, 12, 26];
% ends = [9.563, 13.5, 27.5]
% [U_ND3, freqs_ND3, k_nd3] =cs_exp_ND3.psd_from_intervals('uz', 'scan', starts, ends);
% 
% figure(12)
% h_nd3 = semilogx(freqs_ND3, 10*log10(U_ND3), 'k');
% %%
% h_nd3.DisplayName = sprintf('CS (%d scans)', k_nd3);

%%
% [uz_mat1, n1] = extract_uz_in_time_range(cs_exp_zbounce, 1, 7);
% n = min(n1, n2);
% 
% 
% 
% U_cs_zbounce_psd = 0;
% win = window(@hann, size(uz_mat1,1));
% for k=1:size(uz_mat1,2)
%   u_k = (uz_mat1(:,k));
%   u_k = u_k - mean(u_k);
% % keyboard
% 
%   [u_psd_k, freqs] = power_spectrum(u_k, Ts, win);
% 
%   U_cs_zbounce_psd = U_cs_zbounce_psd + u_psd_k;
% 
% end
% U_cs_zbounce_psd = U_cs_zbounce_psd/k;
% 
% figure(12)
% h4 = semilogx(freqs, 10*log10(U_cs_zbounce_psd));

%%
% % % %%
% % % % Now, lets solve for the CS image and do the same procedure on the
% % % % reconstructed data. Hopefully, this will tell us where we go wrong.
% % % 

% % % 
% % % % equivalent sampling ratio
% % % ts_ = 1/( 512 * 2); % 512 per 0.5 sec
% % % win = window(@hann, 512);
% % % figure(1); clf
% % % for k=[66:77, 219:232, 277:285, 379:391, 482:494]
% % %   row2 = pixmat_cs(512-k, :);
% % % % xq = linspace(1, 512, samps_per_line);  
% % % %   row = interp1([1:512], row2, xq);
% % % %   row = detrend(row);
% % % row = detrend(row2);
% % % if max(abs(row)) > 0.015
% % %   continue
% % % end
% % % figure(1)
% % % plot(row);
% % % hold on
% % %   [U_k, freqs_cs] = power_spectrum(row(:), ts_, win);
% % %   U_cs = U_cs + U_k;
% % %   N = N + 1;
% % % end
% % % U_cs = U_cs/N;
% % % 
% % % figure(12)
% % % h3 = semilogx(freqs_cs, 10*log10(U_cs));
% % % h3.DisplayName = sprintf('CS %d lines. equiv Fs = %.2f', N, 1/ts_)
% % % grid on
% % % 
% % % 
% % % %%





% 
% 
% function [uz, x] = get_row(datmat, row, samps_per_line, npix)
%   samps_per_period = 2*samps_per_line;
%   start_idx = samps_per_period*(row-1) + 1;
%   end_idx = start_idx + samps_per_line;
% 
%   x = datmat(start_idx:end_idx, 1);
%   uz = detrend(datmat(start_idx:end_idx, 4));
%   x = x - x(1);
%   
%   x = (x)*512/x(end);
% end
% 
% function[datmat, samps_period, samps_line] = get_raster_datmat(raster_paths, nperiods)
%   
%   fprintf('Loading Meta file...\n%s\n', raster_paths.meta_path)
%   load(raster_paths.meta_path);
%   fprintf('Loading Data file...\n%s\n', raster_paths.data_path)
%   datmat = csvread(raster_paths.data_path);
%   
%   parent_dat = csvread(raster_paths.parent_path);
%   
%   
%   samps_period = size(parent_dat,1)/2; % twice as many in here for x & y.
%   samps_line = samps_period/2;
%     
%   datmat = datmat([1:nperiods*samps_period], :);
% end

function [figs, axs] = make_traj_figs(figbase)
  Fig_uz = figure(20+figbase); clf
  
  ax1 = gca();
  Fig_ze = figure(30+figbase); clf
  ax2 = gca();
  
  Fig_x = figure(40+figbase); clf
  ax3 = gca();
  Fig_y = figure(50+figbase); clf
  ax4 = gca();
  
  figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
  axs = {ax1, ax2, ax3, ax4};
end
