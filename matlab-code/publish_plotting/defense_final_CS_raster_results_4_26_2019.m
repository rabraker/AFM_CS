clc
clear
%%
close all
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
% initialize paths.
init_paths();
figbase = 20;


size_dir = '5microns';
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
hole_depth = (20);

chan_map = ChannelMap([1:5]);
exp_date = '4-26-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_01Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_04Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_05Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_08Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_10Hz_out_4-26-2019-01.csv',...
};
cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-04Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-05Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-08Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-01Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-02Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-04Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-05Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-08Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
};


%%
rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
  if rast_exps{k}.time_total == 0
      rast_exps{k}.time_total = rast_exps{k}.samps_per_period*rast_exps{k}.npix*AFM.Ts;
  end
end
%

use_ze = false;
x1s = [63, 60, 57, 54, 52];
x2s = [486, 475, 473, 469, 417];
figbase = 10;
%%
for k=1:length(rast_exps)
    rast_exps{k}.bin_raster_really_slow(@detrend, use_ze);
  
  pixmats_raw{k} = rast_exps{k}.pix_mat(1:end, 1:end);
%   rast_exps{k}.pix_mat_pinned = pixmats_raw{k};
  pixmat_ = pin_along_column(rast_exps{k}.pix_mat, x1s(k), x2s(k));
  rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:));
  rast_exps{k}.pin_idx_s = [x1s(k), x2s(k)];
  
      
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit)

  % stit = sprintf('(raw) scan %d', k);
  % plot_raster_data(pixmats_raw{k}, (figbase-5)*k, stit)
end
%%
for k=1:length(rast_exps)
  rast_exps{k}.save()
end
%%
data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg,...
      'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end
%%
[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps{1}.plot_all_cycles(axs{1:4});


%%
bp = true;
use_dct2 = true;
thresh = (20/7)*(1/1000)*20;
opts = l1qc_dct_opts('l1_tol', 0.001);
for k=1:length(cs_exps)
  cs_exps{k}.process_cs_data(false, []);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
    cs_exps{k}.solve_bp(false, use_dct2, opts);

    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, \\%% %.2f sampling\nTotal time: %.2f',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{k}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    f10=mkfig(1001 + 2*k, 6, 7.5); clf
    ax4 = subplot(3,1,[1,2]);
    ax4_2 =subplot(3,1,3);
    
    ImshowDataView.setup(f10);
    cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
    ImshowDataView.imshow(cs_exps{k}.Img_bp-mean(cs_exps{k}.Img_bp(:)),...
      [-thresh, thresh], ax4, ax4_2, cb_exp)
    title(ax4, stit)
  end

%   cs_exps{k}.save();
end
%%
% ------------------- Load Data from ACC ---------------------------------
acc_cs_root = '/media/labserver/acc2018-data/cs-data/5microns/9-22-2017';

cs_name_s = {...
    'cs-traj-512pix-5perc-500nm-5mic-01Hz_out_9-23-2017-04_img-data.mat',...
    'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_9-23-2017-04_img-data.mat',...
    'cs-traj-512pix-15perc-500nm-5mic-01Hz_out_9-23-2017-03_img-data.mat'};

acc_rast_root = '/media/labserver/acc2018-data/raster/5microns/9-22-2017';

% ---------------------- Raster Data ----------------------------------- %
rast_name_s={...
    'raster_scan_512pix_5mic_2.50e-01Hz_out_9-23-2017-02.mat',...
    'raster_scan_512pix_5mic_01Hz_out_9-23-2017-04.mat',...
    'raster_scan_512pix_5mic_05Hz_out_9-23-2017-03.mat',...
    'raster_scan_512pix_5mic_10Hz_out_9-23-2017-01.mat'};

acc_cs_s = {};
for k=1:length(cs_name_s)
    dat = load(fullfile(acc_cs_root, cs_name_s{k}));
    acc_cs_s{k} = dat.img_data;
    acc_cs_s{k}.sub_sample_frac = sum(acc_cs_s{k}.pixelifsampled(:))/512^2;

end

acc_rast_s = {};
x1s = [29, 26, 67, 67];
x2s = [491, 488, 481, 481];
for k=1:length(rast_name_s)
    dat = load(fullfile(acc_rast_root, rast_name_s{k}));
    acc_rast_s{k} = dat.rdat;
    if k < 5
    pix_mat_pinned= pin_along_column(dat.rdat.pixmat, x1s(k), x2s(k));
    else
        pix_mat_pinned= dat.rdat.pixmat;
    end
    pix_mat_pinned = pix_mat_pinned - mean(pix_mat_pinned(:));
    acc_rast_s{k}.pixmat_pinned = pix_mat_pinned;
    stit = sprintf('raster (%.2f Hz)', acc_rast_s{k}.freq);
    plot_raster_data(pix_mat_pinned, 5000+k, stit)
%     figure(5000+k)
%     imagesc(pixmat_pinned, [-thresh, thresh])
%     colormap('gray')
end

%%
Fig_rast = mkfig(2999, 7, 1.8); clf;
ha_rast = tight_subplot(1, 4, [0.025, 0.015], [.01, .1], [.02, .02], true);
Fig_cs = mkfig(3000, 7, 1.8); clf;
ha_cs = tight_subplot(1, 4, [0.025, 0.015], [.01, .1], [.02, .02], true);

Fig_err = mkfig(3001, 7, 9); clf
ha_err = tight_subplot(4, 3, [0.025, 0.015], [.01, .02], [.05, .02], true);

% Fig_rows_acc = mkfig(3003, 7, 2); clf
% ha_row_acc = tight_subplot(1, 1, [0.1, 0.015], [.1, .05], [.085, .02], false);
% ylabel(ha_row_acc(1), 'height [nm]', 'FontSize', 14)
% title(ha_row_acc(1), 'raster (before)', 'FontSize', 14)


Fig_rows = mkfig(3002, 7, 5.5); clf
ha_row = tight_subplot(3, 1, [0.1, 0.015], [.1, .05], [.085, .02], false);

title(ha_row(1), 'raster (before)', 'FontSize', 14)
ylabel(ha_row(1), 'height [nm]', 'FontSize', 14)

title(ha_row(2), 'raster (now)', 'FontSize', 14)
ylabel(ha_row(2), 'height [nm]', 'FontSize', 14)

title(ha_row(3), 'CS', 'FontSize', 14)
ylabel(ha_row(3), 'height [nm]', 'FontSize', 14)
xlabel(ha_row(3), 'x-direction pixel', 'FontSize', 14)

grid(ha_row(1), 'on')
grid(ha_row(2), 'on')
grid(ha_row(3), 'on')

mu = 100;
% mu = Inf;
Img_filts = {};
mxs = [];
thresh = (20/7)*(1/1000)*20;
DRng = 2*thresh;

for k=1:length(cs_exps)
  Img_filts{k} = cs_exps{k}.Img_bp - min(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    Img_filts{k} = SplitBregmanROFAn(Img_filts{k}, mu, 0.001);
  end
%   figure(101+k)
  mx = max(Img_filts{k}(:)) - min(Img_filts{k}(:));
  mxs = [mxs; mx];
%   mesh(Img_filts{k}), colormap('gray')
end

slice = 30:512-30;
master_idx = 1;
im_master = rast_exps{master_idx}.pix_mat_pinned - mean(rast_exps{master_idx}.pix_mat_pinned(:));
if ~isinf(mu)
  %im_master = SplitBregmanROF(im_master, mu, 0.001);
end
fprintf('---------------------------------------------------\n');
scan_metrics = {};
row_idx = 221;
fprintf('\n')
j=1;
excludes = [];
do_rows = [1, 3, 4, 5];
% do_rows = [1, 3, 5];
npix = rast_exps{1}.npix;
for k=do_rows

  imk = rast_exps{k}.pix_mat_pinned - mean(rast_exps{k}.pix_mat_pinned(:));
  if ~isinf(mu)
    %imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
%   [im1_ontok_fit, imk_slice] = align_by_metric(im_master, imk, [], 'psnr');
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);
  damage = rast_exps{k}.damage_metric();
  quality = rast_exps{k}.quality_metric();
  
  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', rast_exps{k}.meta_in.raster_freq,...
    'time', rast_exps{k}.time_total,...
    'coverage', 100,...
    'type', 'raster');
    
  scan_metrics{end+1} = csm;
  
  if ~isempty(intersect(k, excludes))
    continue;
  end

  im_err = im1_ontok_fit - imk_slice;
  imagesc(ha_err(j), im_err, [-thresh, thresh]);
  colormap(ha_err(j), 'gray');
  
  
  imagesc(ha_rast(j), imk, [-thresh, thresh]);
  colormap(ha_rast(j), 'gray');
  stit = sprintf('raster (%.1f Hz, %.1f sec.)',  csm.rate, npix/csm.rate);
  if k==2 
    stit_err = sprintf('(master) raster (%.1f Hz, %.1f sec)', csm.rate, npix/csm.rate);
  else
    stit_err = stit;
  end

  if ~isempty(intersect(k, do_rows)) && csm.rate ~= 8
    hold(ha_rast(j), 'on')
    hold(ha_row(2), 'on')
    
    plot(ha_rast(j), [1, 512], [row_idx, row_idx], 'r');
    hl_row = plot(ha_row(2), imk(row_idx, :)*AFM.volts2nm_z());
    hl_row.DisplayName = sprintf('%.1f Hz', csm.rate);
  end
  title(ha_err(j), stit_err);
  title(ha_rast(j), stit);
  
  fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);

  set(ha_rast(j), 'XTickLabel', [], 'YTickLabel', []);
  set(ha_err(j), 'XTickLabel', [], 'YTickLabel', []);
  
  j=j+1;
end


set(ha_row(2), 'XLim', [1, 512])
axes(ha_row(2))
leg2 = legend();
set(leg2, 'NumColumns', 3, 'FontSize', 11, 'Position', [0.3178 0.6005 0.4223 0.0409]);

% ----------- plot acc rows --------------------------------------------- %

hands_row_acc = gobjects(length(acc_rast_s)-1, 1);
hold(ha_row(1), 'on')
% cla(ha_row_acc)
for k=2:length(acc_rast_s)
%    hands_row_acc(k) = plot(ha_row_acc,...
%        acc_rast_s{k}.pixmat_pinned(row_idx, :)*AFM.volts2nm_z());
   hands_row_acc(k-1) = plot(ha_row(1),...
       acc_rast_s{k}.pixmat_pinned(row_idx, :)*AFM.volts2nm_z());
   hands_row_acc(k-1).DisplayName = sprintf('%.2f Hz', acc_rast_s{k}.freq);
end


xlim(ha_row(1), [1, npix])
leg1 = legend(hands_row_acc);
set(leg1, 'NumColumns', 3, 'FontSize', 11, 'Position', [0.3059 0.9168 0.4550 0.0409])
%
% ------------------------ Plot Current CS ------------------------------ %
j = length(rast_exps);
fprintf('---------------------------------------------------\n');
cs_stats = {};
do_rows = [1, 4, 6, 10];

j = 1;
for k=do_rows

  imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);

  damage = cs_exps{k}.damage_metric();
  quality = cs_exps{k}.quality_metric();

  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
    'time', length(cs_exps{k}.x)*cs_exps{k}.Ts,...
    'coverage', cs_exps{k}.meta_in.actual_sub_samble_perc,...
    'type', 'CS');
  
  fprintf('(CS %.1f Hz, %% %.1f, %.1f) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    csm.rate, csm.coverage, cs_exps{k}.sub_sample_frac*100, psn_1k, ssm_1k, damage, csm.time);

    scan_metrics{end+1} = csm;
  [t_cycle_avg, t_connect] = cs_exps{k}.estimate_mpt_connect_savings();
  cs_stats{k, 1} = csm;
  cs_stats{k, 2} = t_cycle_avg;
  cs_stats{k, 3} = t_connect;
  
  if k==3 || k==7
    continue
  end
  im_err = im1_ontok_fit - imk_slice;
  imagesc(ha_err(j), im_err, [-thresh, thresh]);
  colormap(ha_err(j), 'gray');

  imagesc(ha_cs(j), imk, [-thresh, thresh]);
  colormap(ha_cs(j), 'gray');

  title(ha_cs(j), sprintf('CS (%.1f \\%%, %.1f Hz, %1.f sec.)', csm.coverage, csm.rate, csm.time));
  title(ha_err(j), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));
  
  set(ha_cs(j), 'XTickLabel', [], 'YTickLabel', []);
  set(ha_err(j), 'XTickLabel', [], 'YTickLabel', []);

  if ~isempty(intersect(k, do_rows))
    hold(ha_cs(j), 'on')
    hold(ha_row(3), 'on')
    
    plot(ha_cs(j), [1, 512], [row_idx, row_idx], 'r');
    hl_row = plot(ha_row(3), imk(row_idx, :)*AFM.volts2nm_z());
    hl_row.DisplayName = sprintf('%.1f Hz, %.1f \\%%', csm.rate, csm.coverage);
  end
  

  
  j=j+1;
end

set(ha_row(3), 'XLim', [1, 512])
axes(ha_row(3))
leg3 = legend();
set(leg3, 'NumColumns', 4, 'FontSize', 11, 'Position', [0.1133 0.2728 0.8589 0.0409]);

set(ha_row(1), 'YLim', [-20, 15]);
set(ha_row(2), 'YLim', [-20, 15]);
set(ha_row(3), 'YLim', [-20, 15]);
%%


save_fig(Fig_rast, fullfile(PATHS.defense_fig(), 'raster_images_3-20-2019'), true)
save_fig(Fig_cs, fullfile(PATHS.defense_fig(), 'cs_images_3-20-2019'), true)
% save_fig(Fig_err, 'notes/figures/cs_raster_images_err_3-20-2019')
save_fig(Fig_rows, fullfile(PATHS.defense_fig(), 'cs_raster_pixel_rows_3-20-2019'), true)



%%
Fig_acc = mkfig(3999, 5.5, 1.8); clf;
ha_acc = tight_subplot(1, 3, [0.025, 0.015], [.01, .1], [.02, .02], true);
Fig_now = mkfig(4000, 5.5, 1.8); clf;
ha_now = tight_subplot(1, 3, [0.025, 0.015], [.01, .1], [.02, .02], true);

Fig_now_tv = mkfig(4001, 5.5, 1.8); clf;
ha_now_tv = tight_subplot(1, 3, [0.025, 0.015], [.01, .1], [.02, .02], true);

Fig_row_tv = mkfig(4002, 5.5, 1.); clf;
ha_row_tv = tight_subplot(1, 1, [0.025, 0.015], [.01, .1], [.1, .02], true);
Fig_row = mkfig(4003, 5.5, 1.); clf;
ha_row = tight_subplot(1, 1, [0.025, 0.015], [.01, .1], [.1, .02], true);

row_idx = 220;

imagesc(ha_acc(1), acc_rast_s{2}.pixmat_pinned, [-thresh, thresh])
colormap(ha_acc(1), 'gray')
title(ha_acc(1), 'raster (1.0 Hz, 512 sec)')

imagesc(ha_acc(2), acc_cs_s{2}.bp_im, [-thresh, thresh])
colormap(ha_acc(2), 'gray')
ttot = sum(acc_cs_s{2}.meta.state_counts) * AFM.Ts;
stit = sprintf('CS (1.0 Hz, \\%% %.2f, %.1f sec)',...
    acc_cs_s{2}.sub_sample_frac*100, ttot);
title(ha_acc(2), stit)

imagesc(ha_acc(3), acc_cs_s{3}.bp_im, [-thresh, thresh])
colormap(ha_acc(3), 'gray')
ttot = sum(acc_cs_s{3}.meta.state_counts) * AFM.Ts;
stit = sprintf('CS (1.0 Hz, \\%% %.2f, %.1f sec)',...
    acc_cs_s{3}.sub_sample_frac*100, ttot);
title(ha_acc(3), stit)


row_acc_s{1} = acc_rast_s{2}.pixmat_pinned(row_idx, :);
row_acc_s{2} = acc_cs_s{2}.bp_im(row_idx, :);

hold(ha_row, 'on')
plot(ha_row, row_acc_s{1})
plot(ha_row, row_acc_s{2})


% -------------------- Current Results -----------------------------------
% Plot both with and without TV-denoising.
mu = 100; 

imk = rast_exps{1}.pix_mat_pinned - mean(rast_exps{1}.pix_mat_pinned(:));
imagesc(ha_now(1), imk, [-thresh, thresh])
colormap(ha_now(1), 'gray')
title(ha_now(1), 'raster (1.0 Hz, 512 sec)')

imagesc(ha_now_tv(1), imk, [-thresh, thresh])
colormap(ha_now_tv(1), 'gray')
title(ha_now_tv(1), 'raster (1.0 Hz, 512 sec)')

k = 1;
imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
if ~isinf(mu)
    imk_tv = SplitBregmanROF(imk, mu, 0.001);
end
stit = sprintf('CS (%.1f Hz, \\%% %.2f, %.1f sec)',cs_exps{k}.meta_in.tip_velocity/10,...
    cs_exps{k}.sub_sample_frac*100, cs_exps{k}.time_total);

imagesc(ha_now(2), imk, [-thresh, thresh])
colormap(ha_now(2), 'gray')
title(ha_now(2), stit)

imagesc(ha_now_tv(2), imk_tv, [-thresh, thresh])
colormap(ha_now_tv(2), 'gray')
title(ha_now_tv(2), stit)

plot(ha_row, imk(row_idx, :))

% 1Hz and 15%
k = 6;
imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
if ~isinf(mu)
    imk_tv = SplitBregmanROF(imk, mu, 0.001);
end

stit = sprintf('CS (%.1f Hz, \\%% %.2f, %.1f sec)',cs_exps{k}.meta_in.tip_velocity/10,...
    cs_exps{k}.sub_sample_frac*100, cs_exps{k}.time_total);

imagesc(ha_now(3), imk, [-thresh, thresh])
colormap(ha_now(3), 'gray')
title(ha_now(3), stit)


imagesc(ha_now_tv(3), imk_tv, [-thresh, thresh])
colormap(ha_now_tv(3), 'gray')
title(ha_now_tv(3), stit)
%%
for j=1:3
  set(ha_acc(j), 'XTickLabel', [], 'YTickLabel', []);
  set(ha_now(j), 'XTickLabel', [], 'YTickLabel', []);
  set(ha_now_tv(j), 'XTickLabel', [], 'YTickLabel', []);
end

%%
save_fig(Fig_acc, fullfile(PATHS.defense_fig(), 'final_cs_rast_acc_for_cp'), true)
save_fig(Fig_now, fullfile(PATHS.defense_fig(), 'final_cs_rast_now_for_cp'), true)
save_fig(Fig_now_tv, fullfile(PATHS.defense_fig(), 'final_csTV_rast_now_for_cp'), true)




%%
Fig_rows = mkfig(3002, 7, 4.5); clf
ha_row = tight_subplot(2, 1, [0.1, 0.015], [.1, .05], [.085, .02], false);
xlabel(ha_row(2), 'x-direction pixel', 'FontSize', 14)
ylabel(ha_row(1), 'height [nm]', 'FontSize', 14)
title(ha_row(1), 'raster', 'FontSize', 14)
title(ha_row(2), 'CS', 'FontSize', 14)
ylabel(ha_row(2), 'height [nm]', 'FontSize', 14)

grid(ha_row(1), 'on')
grid(ha_row(2), 'on')


%%

S1 = scan_metrics_table(scan_metrics)

S2 = state_times_table(cs_exps)



fid = fopen('notes/tables/cs_raster_table_3-20-2019_muInf_dct2.tex', 'w+');
fprintf(fid, '%s', S1);
fclose(fid);

fid = fopen('notes/tables/cs_state_times_table_3-20-2019_muInf_dct2.tex', 'w+');
fprintf(fid, '%s', S2);
fclose(fid);


function S = scan_metrics_table(csm_s)
  S = 'type &  rate (Hz) & PSNR & SSIM & time [s] & damage\\';
  S = sprintf('%s\n\\toprule\n', S);
  
  for k=1:length(csm_s)
    csm = csm_s{k};
    if strcmp(csm.type, 'CS')
      description = sprintf('CS (%.2f~\\%%)', csm.coverage);
    else
      description = 'raster';
    end
    S = sprintf('%s%s & %.2f & %.2f & %.2f & %.1f & %.2f\\\\\n', S,...
      description, csm.rate, csm.psnr, csm.ssim, csm.time, csm.damage);
  end
end

function S = state_times_table(cs_exps)
  S = 'description & move & engage & tip-settle & scan & connect & tip-up & total\\';
  S = sprintf('%s\n\\toprule\n', S);
  for k=1:length(cs_exps)
    cst = cs_exps{k}.get_state_times();
    description = sprintf('%.1f Hz, %.2f~\\%%',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc);
  
    S = sprintf('%s%s &%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\\n', S,...
      description, cst.move, cst.tdown, cst.tsettle, cst.scan,...
      cst.connect, cst.tup, cst.total);
  end
end




%%
function plot_raster_data(pixmat2, figbase, stit, plot_mesh)
  if nargin < 4
    plot_mesh = false;
  end
  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
%   F10 = figure(figbase+2); clf
  F10 = mkfig(figbase+2, 6, 7.5, false);
  ax1 = subplot(3,1,[1,2]);
  ax2 =subplot(3,1,3);
    

  lo = min(min(pixmat2));
  hi = max(max(pixmat2));
  
  
  imshow_dataview(flipud(pixmat2 - mean(pixmat2(:))), [-thresh, thresh], ax1, ax2)
  try
    grid(ax1, 'on')
  catch
    keyboard
  end
  
  colormap(ax1, 'gray')
  grid(ax2, 'on')
  ax1.GridAlpha = 1;
  ax2.GridAlpha = 1;
  title(ax1, stit)
  title(ax2, stit)
  if plot_mesh
    figure(figbase+4)
    mesh(pixmat2)
    xlabel('x')
    ylabel('y')
    title(stit)
  end
end
