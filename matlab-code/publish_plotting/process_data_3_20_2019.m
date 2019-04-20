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
exp_date = '3-20-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_04Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_05Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_08Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_10Hz_out_3-20-2019-01.csv',...
};
cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-04Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-05Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-02Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-04Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-05Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
};


%%
rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths, 'load_full', true);
  if rast_exps{k}.time_total == 0
      rast_exps{k}.time_total = rast_exps{k}.samps_per_period*rast_exps{k}.npix*AFM.Ts;
  end
end
%

x1s = [30, 47, 53, 51, 59, 60];
x2s = [496, 465, 464, 468, 472, 424];
figbase = 10;
for k=1:length(rast_exps)
  rast_exps{k}.bin_raster_really_slow(@detrend);
  
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
for k=1:length(cs_exps)
  cs_exps{k}.save()
end
%%
[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps{1}.plot_all_cycles(axs{1:4});


%%
bp = true;
recalc = false;
use_dct2 = true;
thresh = (20/7)*(1/1000)*20;
opts = l1qc_dct_opts('l1_tol', 0.001);
for k=1:length(cs_exps)
  cs_exps{k}.process_cs_data(false, []);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
    cs_exps{k}.solve_bp(recalc, use_dct2, opts);

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
Fig1 = mkfig(3000, 7, 9); clf;
ha1 = tight_subplot(4, 3, [0.025, 0.015], [.01, .02], [.02, .02], true);

Fig_err = mkfig(3001, 7, 9); clf
ha_err = tight_subplot(4, 3, [0.025, 0.015], [.01, .02], [.05, .02], true);
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

mu = 300;
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
master_idx = 2;
im_master = rast_exps{master_idx}.pix_mat_pinned - mean(rast_exps{master_idx}.pix_mat_pinned(:));
if ~isinf(mu)
  %im_master = SplitBregmanROF(im_master, mu, 0.001);
end
fprintf('---------------------------------------------------\n');
scan_metrics = {};
row_idx = 220;
fprintf('\n')
j=1;
excludes = [1];
do_rows = [2, 4, 6];
for k=1:length(rast_exps)
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
  
  
  imagesc(ha1(j), imk, [-thresh, thresh]);
  colormap(ha1(j), 'gray');
  stit = sprintf('raster (%.1f Hz)', csm.rate);
  if k==2 
    stit_err = sprintf('(master) raster (%.1f Hz)', csm.rate);
  else
    stit_err = stit;
  end

  if ~isempty(intersect(k, do_rows))
    hold(ha1(j), 'on')
    hold(ha_row(1), 'on')
    
    plot(ha1(j), [1, 512], [row_idx, row_idx], 'r');
    hl_row = plot(ha_row(1), imk(row_idx, :)*AFM.volts2nm_z());
    hl_row.DisplayName = sprintf('%.1f Hz', csm.rate);
  end
  title(ha_err(j), stit_err);
  title(ha1(j), stit);
  
  fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);

  set(ha1(j), 'XTickLabel', [], 'YTickLabel', []);
  set(ha_err(j), 'XTickLabel', [], 'YTickLabel', []);
  
  j=j+1;
end


set(ha_row(1), 'XLim', [1, 512])
axes(ha_row(1))
leg1 = legend();
set(leg1, 'NumColumns', 3, 'FontSize', 11, 'Position', [0.5574 0.8968 0.4223 0.0500])


j = length(rast_exps);
fprintf('---------------------------------------------------\n');
cs_stats = {};
do_rows = [1, 4, 5, 8];
for k=1:length(cs_exps)

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

  imagesc(ha1(j), imk, [-thresh, thresh]);
  colormap(ha1(j), 'gray');

  title(ha1(j), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));
  title(ha_err(j), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));
  
  set(ha1(j), 'XTickLabel', [], 'YTickLabel', []);
  set(ha_err(j), 'XTickLabel', [], 'YTickLabel', []);

  if ~isempty(intersect(k, do_rows))
    hold(ha1(j), 'on')
    hold(ha_row(2), 'on')
    
    plot(ha1(j), [1, 512], [row_idx, row_idx], 'r');
    hl_row = plot(ha_row(2), imk(row_idx, :)*AFM.volts2nm_z());
    hl_row.DisplayName = sprintf('%.1f Hz, %.1f \\%%', csm.rate, csm.coverage);
  end
  

  
  j=j+1;
end

set(ha_row(2), 'XLim', [1, 512])
axes(ha_row(2))
leg2 = legend();
set(leg2, 'NumColumns', 4, 'FontSize', 11, 'Position', [0.1032 0.4253 0.8589 0.0500])
%%

set(ha1(end), 'Visible', 'off')
set(ha_err(end), 'Visible', 'off')


% save_fig(Fig1, 'notes/figures/cs_raster_images_3-20-2019')
% save_fig(Fig_err, 'notes/figures/cs_raster_images_err_3-20-2019')

save_fig(Fig1, fullfile(PATHS.cs_final_fig(), 'cs_raster_images_3-20-2019'), false)
save_fig(Fig_err, fullfile(PATHS.cs_final_fig(), 'cs_raster_images_err_3-20-2019'),false)
save_fig(Fig_rows, fullfile(PATHS.cs_final_fig(), 'cs_raster_pixel_rows_3-20-2019'), false)





%%

% skip the 0.5Hz scan
sm = scan_metrics(2:end)
S1 = scan_metrics_table(sm)

S2 = state_times_table(cs_exps)



fid = fopen(fullfile(PATHS.cs_final_table(), 'cs_raster_table_3-20-2019_muInf_dct2.tex'), 'w+');
fprintf(fid, '%s', S1);
fclose(fid);

fid = fopen(fullfile(PATHS.cs_final_table(), 'cs_state_times_table_3-20-2019_muInf_dct2.tex'), 'w+');
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
